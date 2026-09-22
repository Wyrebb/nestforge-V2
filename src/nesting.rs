use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};
use std::sync::atomic::{AtomicBool, Ordering};
use rand::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use crate::geometry::*;

#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct Piece { pub id:u32, pub cutting_hull:Vec<Point>, pub qi:usize }

#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct Placement {
    pub piece_id:u32, pub orig_qi:usize, pub deg:f64,
    pub tx:f64, pub ty:f64, pub hull:Vec<Point>, pub xf:Xform,
}

#[derive(Debug,Clone,Serialize,Deserialize)]
pub struct NestConfig {
    pub sheet_w:f64, pub sheet_h:f64, pub margin:f64, pub spacing:f64,
    pub rotations:Vec<f64>, pub ga_pop:usize, pub ga_gens:usize,
    pub continuous:bool, pub max_sheets:usize,
    #[serde(default="default_strategy")] pub strategy:String,
}
fn default_strategy()->String{"left".to_string()}

#[derive(Debug,Serialize,Deserialize)]
pub struct SheetResult {
    pub sheet_num:usize, pub placements:Vec<Placement>,
    pub fitness:f64, pub util_pct:u32,
}
pub type Cb=Arc<dyn Fn(usize,Vec<Placement>,f64,u32)+Send+Sync>;

// ── Cache NFP global persistant entre tous les runs ───────────
lazy_static::lazy_static! {
    static ref GLOBAL_NFP_CACHE: RwLock<HashMap<String,Vec<Point>>> =
        RwLock::new(HashMap::new());
}

// ── Douglas-Peucker ───────────────────────────────────────────
fn dp_dist(p:Point, a:Point, b:Point) -> f64 {
    let dx=b[0]-a[0]; let dy=b[1]-a[1]; let len2=dx*dx+dy*dy;
    if len2<1e-12 { return ((p[0]-a[0]).powi(2)+(p[1]-a[1]).powi(2)).sqrt(); }
    let t=((p[0]-a[0])*dx+(p[1]-a[1])*dy)/len2;
    let t=t.clamp(0.0,1.0);
    ((p[0]-a[0]-t*dx).powi(2)+(p[1]-a[1]-t*dy).powi(2)).sqrt()
}

fn dp_recurse(pts:&[Point], eps:f64, i0:usize, i1:usize) -> Vec<usize> {
    if i1<=i0+1 { return vec![]; }
    let mut max_d=0.0f64; let mut max_i=i0+1;
    for i in i0+1..i1 { let d=dp_dist(pts[i],pts[i0],pts[i1]); if d>max_d{max_d=d;max_i=i;} }
    if max_d>eps {
        let mut r=dp_recurse(pts,eps,i0,max_i);
        r.push(max_i);
        r.extend(dp_recurse(pts,eps,max_i,i1));
        r
    } else { vec![] }
}

pub fn simplify_poly(pts:&[Point], eps:f64) -> Vec<Point> {
    if pts.len()<=3 { return pts.to_vec(); }
    let n=pts.len();
    let mut open=pts.to_vec(); open.push(pts[0]);
    let mut idx=vec![0usize];
    idx.extend(dp_recurse(&open,eps,0,n));
    idx.push(n);
    idx.sort(); idx.dedup();
    let result:Vec<Point>=idx.iter().filter(|&&i|i<n).map(|&i|pts[i]).collect();
    if result.len()>=3 { result } else { pts.to_vec() }
}

// ── NFP orbital (Burke et al. 2006) ──────────────────────────
fn nfp_orbital(a:&[Point], b:&[Point]) -> Option<Vec<Point>> {
    let pa=ensure_cw(a); let pb=ensure_cw(b);
    let mut sa=0usize;
    for i in 1..pa.len() { if pa[i][1]<pa[sa][1]||(almost_equal(pa[i][1],pa[sa][1])&&pa[i][0]<pa[sa][0]){sa=i;} }
    let mut sb=0usize;
    for i in 1..pb.len() { if pb[i][1]>pb[sb][1]||(almost_equal(pb[i][1],pb[sb][1])&&pb[i][0]>pb[sb][0]){sb=i;} }
    let mut tx=pa[sa][0]-pb[sb][0]; let mut ty=pa[sa][1]-pb[sb][1];
    let maxiter=(pa.len()+pb.len()+4)*6;
    let mut nfp:Vec<Point>=Vec::with_capacity(maxiter);
    let mut prev:Option<[f64;2]>=None;
    for _ in 0..maxiter {
        nfp.push([tx,ty]);
        struct E{x:f64,y:f64,nx:f64,ny:f64}
        let mut edges:Vec<E>=Vec::new();
        for i in 0..pa.len() { let j=(i+1)%pa.len(); let ex=pa[j][0]-pa[i][0]; let ey=pa[j][1]-pa[i][1]; let l=(ex*ex+ey*ey).sqrt(); if l>1e-12{edges.push(E{x:ex,y:ey,nx:ex/l,ny:ey/l});} }
        for i in 0..pb.len() { let j=(i+1)%pb.len(); let ex=-(pb[j][0]-pb[i][0]); let ey=-(pb[j][1]-pb[i][1]); let l=(ex*ex+ey*ey).sqrt(); if l>1e-12{edges.push(E{x:ex,y:ey,nx:ex/l,ny:ey/l});} }
        let best=match prev {
            None => edges.iter().min_by(|a,b|a.nx.partial_cmp(&b.nx).unwrap().then(a.ny.partial_cmp(&b.ny).unwrap())),
            Some(pv) => edges.iter().min_by(|a,b|{
                let ang=|e:&&E|{let cr=pv[0]*e.ny-pv[1]*e.nx; let dt=pv[0]*e.nx+pv[1]*e.ny; cr.atan2(dt)};
                ang(a).partial_cmp(&ang(b)).unwrap()
            }),
        };
        let e=match best{Some(e)=>e, None=>break};
        prev=Some([e.nx,e.ny]); tx+=e.x; ty+=e.y;
        if nfp.len()>3&&almost_equal(tx,nfp[0][0])&&almost_equal(ty,nfp[0][1]){break;}
    }
    if nfp.len()>=3{Some(nfp)}else{None}
}

fn nfp_minkowski(a:&[Point], b:&[Point]) -> Vec<Point> {
    let ha=convex_hull(a); let bn:Vec<Point>=convex_hull(b).iter().map(|p|[-p[0],-p[1]]).collect();
    let pa=ensure_ccw(&ha); let pb=ensure_ccw(&bn);
    let mut ia=0usize; let mut ib=0usize;
    for i in 1..pa.len(){if pa[i][1]<pa[ia][1]||(almost_equal(pa[i][1],pa[ia][1])&&pa[i][0]<pa[ia][0]){ia=i;}}
    for i in 1..pb.len(){if pb[i][1]<pb[ib][1]||(almost_equal(pb[i][1],pb[ib][1])&&pb[i][0]<pb[ib][0]){ib=i;}}
    let na=pa.len(); let nb=pb.len();
    let mut res:Vec<Point>=Vec::new(); let mut i=0usize; let mut j=0usize;
    while i<na||j<nb {
        res.push([pa[(ia+i)%na][0]+pb[(ib+j)%nb][0], pa[(ia+i)%na][1]+pb[(ib+j)%nb][1]]);
        let ea=[pa[(ia+i+1)%na][0]-pa[(ia+i)%na][0], pa[(ia+i+1)%na][1]-pa[(ia+i)%na][1]];
        let eb=[pb[(ib+j+1)%nb][0]-pb[(ib+j)%nb][0], pb[(ib+j+1)%nb][1]-pb[(ib+j)%nb][1]];
        let cr=ea[0]*eb[1]-ea[1]*eb[0];
        if i>=na{j+=1;}else if j>=nb{i+=1;}else if cr>0.0{i+=1;j+=1;}
        else if cr<0.0{if ea[1]<eb[1]||(almost_equal(ea[1],eb[1])&&ea[0]<eb[0]){i+=1;}else{j+=1;}}
        else{i+=1;j+=1;}
    }
    res
}

fn nfp_key(ia:u32, da:f64, ib:u32, db:f64, sp:f64) -> String {
    format!("{}_{}_{:04}_{:04}_{:03}", ia, ib, (da*10.0)as i32, (db*10.0)as i32, (sp*10.0)as i32)
}

// ── Précalcul NFP parallèle avec cache global ─────────────────
pub fn precalc_nfps(queue:&[Piece], cfg:&NestConfig) {
    let mut seen=HashSet::new();
    let uniq:Vec<&Piece>=queue.iter().filter(|p|seen.insert(p.id)).collect();

    // Trouver les paires manquantes dans le cache global
    let pairs:Vec<(u32,Vec<Point>,f64,u32,Vec<Point>,f64,f64)> = {
        let cache=GLOBAL_NFP_CACHE.read().unwrap();
        let mut v=Vec::new();
        for pa in &uniq { for pb in &uniq {
            for &da in &cfg.rotations { for &db in &cfg.rotations {
                let key=nfp_key(pa.id,da,pb.id,db,cfg.spacing);
                if !cache.contains_key(&key){
                    v.push((pa.id, pa.cutting_hull.clone(), da,
                            pb.id, pb.cutting_hull.clone(), db, cfg.spacing));
                }
            }}
        }}
        v
    };

    if pairs.is_empty() { return; }

    // Calcul PARALLÈLE de toutes les paires manquantes
    let computed:Vec<(String,Vec<Point>)>=pairs.par_iter().map(|(ia,ha,da,ib,hb,db,sp)|{
        let key=nfp_key(*ia,*da,*ib,*db,*sp);
        let norm_a=simplify_poly(&norm(&rotate(&center(ha),*da)), 0.1);
        let norm_b_raw=if *sp>0.0{offset_poly(&norm(&rotate(&center(hb),*db)),*sp/2.0)}else{norm(&rotate(&center(hb),*db))};
        let norm_b=simplify_poly(&norm_b_raw, 0.1);
        let nfp=nfp_orbital(&norm_a,&norm_b).unwrap_or_else(||nfp_minkowski(&norm_a,&norm_b));
        (key, nfp)
    }).collect();

    // Écriture dans le cache global (lock unique)
    let mut cache=GLOBAL_NFP_CACHE.write().unwrap();
    for (k,v) in computed { cache.insert(k,v); }
}

// ── Placement ─────────────────────────────────────────────────
fn place_piece(piece:&Piece, placed:&[Placement], cfg:&NestConfig, deg:f64) -> Option<(f64,f64,Vec<Point>)> {
    let rh=norm(&rotate(&center(&piece.cutting_hull),deg));
    let b=bbox(&rh); let pw=b.w(); let ph=b.h();
    let x0=cfg.margin; let y0=cfg.margin;
    let x1=cfg.sheet_w-cfg.margin-pw; let y1=cfg.sheet_h-cfg.margin-ph;
    if x1<x0||y1<y0{return None;}

    let pbuf:Vec<(Vec<Point>,BBox)>=placed.iter().map(|pl|(
        if cfg.spacing>0.0{offset_poly(&pl.hull,cfg.spacing/2.0)}else{pl.hull.clone()},
        bbox(&pl.hull)
    )).collect();

    let mut cands:Vec<[f64;2]>=Vec::new();
    if placed.is_empty(){
        cands.push([x0,y0]);
    } else {
        let bcx=pw/2.0; let bcy=ph/2.0;
        let cache=GLOBAL_NFP_CACHE.read().unwrap();
        for pl in placed {
            let key=nfp_key(pl.piece_id,pl.deg,piece.id,deg,cfg.spacing);
            if let Some(nfp)=cache.get(&key){
                let pb2=bbox(&pl.hull);
                let plcx=pl.tx+(pb2.max_x-pb2.min_x)/2.0;
                let plcy=pl.ty+(pb2.max_y-pb2.min_y)/2.0;
                for &[nx,ny] in nfp {
                    let tx=nx+plcx-bcx; let ty=ny+plcy-bcy;
                    cands.push([tx,ty]); cands.push([tx+0.001,ty]); cands.push([tx,ty+0.001]);
                }
            }
        }
        drop(cache);
        let step=(pw.min(ph)*0.08).max(1.5);
        let mut ty=y0; while ty<=y1{let mut tx=x0;while tx<=x1{cands.push([tx,ty]);tx+=step*3.0;}ty+=step*3.0;}
    }

    let mut seen_map:HashMap<String,bool>=HashMap::new();
    let mut valid:Vec<[f64;2]>=Vec::new();
    for [tx,ty] in cands {
        let cx=tx.max(x0).min(x1); let cy=ty.max(y0).min(y1);
        let k=format!("{:.2},{:.2}",cx,cy);
        if seen_map.insert(k,true).is_none(){valid.push([cx,cy]);}
    }

    match cfg.strategy.as_str(){
        "left"  =>valid.sort_by(|a,b|a[0].partial_cmp(&b[0]).unwrap().then(a[1].partial_cmp(&b[1]).unwrap())),
        "corner"=>valid.sort_by(|a,b|(a[0]+a[1]).partial_cmp(&(b[0]+b[1])).unwrap()),
        _       =>valid.sort_by(|a,b|a[1].partial_cmp(&b[1]).unwrap().then(a[0].partial_cmp(&b[0]).unwrap())),
    }

    let mut best:Option<(f64,f64,Vec<Point>)>=None; let mut best_score=f64::MAX;
    for [tx,ty] in valid {
        let cbb=BBox{min_x:tx,min_y:ty,max_x:tx+pw,max_y:ty+ph};
        let needs=pbuf.iter().any(|(_,bb2)|cbb.overlaps(bb2));
        let hull=translate(&rh,tx,ty);
        let ok=if needs{
            let bufh=if cfg.spacing>0.0{offset_poly(&hull,cfg.spacing/2.0)}else{hull.clone()};
            !pbuf.iter().any(|(bp,_)|polys_overlap(&bufh,bp))
        }else{true};
        if ok{
            let gx=placed.iter().fold(tx+pw,|acc,pl|acc.max(bbox(&pl.hull).max_x));
            let gy=placed.iter().fold(ty+ph,|acc,pl|acc.max(bbox(&pl.hull).max_y));
            let score=match cfg.strategy.as_str(){
                "left"  =>gx*cfg.sheet_h+ty,
                "corner"=>gx*cfg.sheet_h+gy*cfg.sheet_w,
                _       =>gy*cfg.sheet_w+tx,
            };
            if score<best_score{best_score=score;best=Some((tx,ty,hull));}
        }
    }
    best
}

// ── Évaluation individu ───────────────────────────────────────
fn eval(perm:&[usize], rots:&[f64], queue:&[Piece], cfg:&NestConfig) -> (Vec<Placement>,f64) {
    let mut placed:Vec<Placement>=Vec::new();
    for (k,&idx) in perm.iter().enumerate(){
        let piece=&queue[idx]; let deg=rots[k];
        if let Some((tx,ty,hull))=place_piece(piece,&placed,cfg,deg){
            let xf=build_xform(&piece.cutting_hull,deg,tx,ty);
            placed.push(Placement{piece_id:piece.id,orig_qi:piece.qi,deg,tx,ty,hull,xf});
        }
    }
    if placed.is_empty(){return(placed,0.0);}
    let all_pts:Vec<Point>=placed.iter().flat_map(|p|p.hull.iter().cloned()).collect();
    let gbb=bbox(&all_pts);
    let ua:f64=placed.iter().map(|p|poly_area(&p.hull)).sum();
    let fitness=(ua/gbb.area().max(1.0))*(placed.len()as f64/perm.len()as f64);
    (placed,fitness)
}

// ── GA ────────────────────────────────────────────────────────
#[derive(Clone)] struct Ind{perm:Vec<usize>,rots:Vec<f64>}

fn rand_perm(n:usize, rng:&mut impl Rng) -> Vec<usize> {
    let mut v:Vec<usize>=(0..n).collect(); v.shuffle(rng); v
}

fn pmx(p1:&[usize], p2:&[usize], rng:&mut impl Rng) -> Vec<usize> {
    let n=p1.len(); let a=rng.gen_range(0..n); let b=rng.gen_range(0..n);
    let lo=a.min(b); let hi=a.max(b);
    let mut c=p1.to_vec(); let mut map:HashMap<usize,usize>=HashMap::new();
    for i in lo..=hi{map.insert(p2[i],p1[i]);c[i]=p2[i];}
    for i in 0..n{
        if i>=lo&&i<=hi{continue;}
        let mut v=c[i]; let mut s=0;
        while let Some(&nv)=map.get(&v){v=nv;s+=1;if s>=n{break;}}
        c[i]=v;
    }
    c
}

pub fn nest_ga(queue:&[Piece], cfg:&NestConfig, sheet:usize, stop:&AtomicBool, cb:&Cb) -> (Vec<Placement>,f64) {
    let n=queue.len(); let mut rng=rand::thread_rng();
    let mut pop:Vec<Ind>=vec![Ind{perm:(0..n).collect(),rots:vec![cfg.rotations[0];n]}];
    for _ in 1..cfg.ga_pop {
        pop.push(Ind{
            perm:rand_perm(n,&mut rng),
            rots:(0..n).map(|_|cfg.rotations[rng.gen_range(0..cfg.rotations.len())]).collect(),
        });
    }
    let mut best_pl:Vec<Placement>=Vec::new(); let mut best_fit=-1.0f64; let mut gen=0usize;

    loop {
        gen+=1;
        // Évaluation PARALLÈLE — le cache global est partagé en lecture seule
        let mut scored:Vec<(Ind,Vec<Placement>,f64)>=pop.par_iter().map(|ind|{
            let (pl,fit)=eval(&ind.perm,&ind.rots,queue,cfg);
            (ind.clone(),pl,fit)
        }).collect();
        scored.sort_by(|a,b|b.2.partial_cmp(&a.2).unwrap());

        if scored[0].2>best_fit {
            best_fit=scored[0].2; best_pl=scored[0].1.clone();
            let util=(best_pl.iter().map(|p|poly_area(&p.hull)).sum::<f64>()/(cfg.sheet_w*cfg.sheet_h)*100.0)as u32;
            cb(sheet,best_pl.clone(),best_fit,util);
        }

        if stop.load(Ordering::Relaxed){break;}
        if !cfg.continuous&&gen>=cfg.ga_gens{break;}

        const ELITE:usize=2; const MUT:f64=0.15; const RMUT:f64=0.25;
        let mut np:Vec<Ind>=scored.iter().take(ELITE).map(|(ind,_,_)|ind.clone()).collect();
        while np.len()<cfg.ga_pop {
            let pick=|rng:&mut ThreadRng|->Ind{
                let i=rng.gen_range(0..scored.len()); let j=rng.gen_range(0..scored.len());
                if scored[i].2>=scored[j].2{scored[i].0.clone()}else{scored[j].0.clone()}
            };
            let p1=pick(&mut rng); let p2=pick(&mut rng);
            let mut cp=pmx(&p1.perm,&p2.perm,&mut rng);
            let mut cr:Vec<f64>=p1.rots.iter().zip(p2.rots.iter()).map(|(&a,&b)|if rng.gen::<f64>()<0.5{a}else{b}).collect();
            for i in 0..cp.len(){if rng.gen::<f64>()<MUT{let j=rng.gen_range(0..cp.len());cp.swap(i,j);}}
            for r in cr.iter_mut(){if rng.gen::<f64>()<RMUT{*r=cfg.rotations[rng.gen_range(0..cfg.rotations.len())];}}
            np.push(Ind{perm:cp,rots:cr});
        }
        pop=np;
    }
    (best_pl,best_fit)
}

// ── Orchestrateur ─────────────────────────────────────────────
pub fn run_nesting(mut queue:Vec<Piece>, cfg:NestConfig, stop:Arc<AtomicBool>, cb:Cb) -> Vec<SheetResult> {
    queue.sort_by(|a,b|poly_area(&b.cutting_hull).partial_cmp(&poly_area(&a.cutting_hull)).unwrap());
    // Précalcul NFP parallèle — saute les entrées déjà dans le cache global
    precalc_nfps(&queue,&cfg);
    let mut results:Vec<SheetResult>=Vec::new();
    let mut remaining:Vec<Piece>=queue; let mut sheet=0;
    while !remaining.is_empty()&&sheet<cfg.max_sheets&&!stop.load(Ordering::Relaxed){
        sheet+=1;
        let (placed,fit)=nest_ga(&remaining,&cfg,sheet,&stop,&cb);
        if placed.is_empty(){break;}
        let util=(placed.iter().map(|p|poly_area(&p.hull)).sum::<f64>()/(cfg.sheet_w*cfg.sheet_h)*100.0)as u32;
        let used:HashSet<usize>=placed.iter().map(|p|p.orig_qi).collect();
        remaining.retain(|p|!used.contains(&p.qi));
        results.push(SheetResult{sheet_num:sheet,placements:placed,fitness:fit,util_pct:util});
    }
    results
}
