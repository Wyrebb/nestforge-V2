use serde::{Deserialize, Serialize};
pub type Point = [f64; 2];
pub type Poly = Vec<Point>;
const TOL: f64 = 1e-9;

pub fn almost_equal(a: f64, b: f64) -> bool { (a - b).abs() < TOL }

pub fn signed_area(pts: &[Point]) -> f64 {
    let n = pts.len(); let mut a = 0.0f64;
    for i in 0..n { let j=(i+1)%n; a += pts[i][0]*pts[j][1] - pts[j][0]*pts[i][1]; }
    a / 2.0
}
pub fn poly_area(pts: &[Point]) -> f64 { signed_area(pts).abs() }

#[derive(Debug, Clone, Copy)]
pub struct BBox { pub min_x:f64, pub min_y:f64, pub max_x:f64, pub max_y:f64 }
impl BBox {
    pub fn w(&self) -> f64 { self.max_x - self.min_x }
    pub fn h(&self) -> f64 { self.max_y - self.min_y }
    pub fn area(&self) -> f64 { self.w() * self.h() }
    pub fn overlaps(&self, o: &BBox) -> bool {
        !(self.max_x < o.min_x || o.max_x < self.min_x
        || self.max_y < o.min_y || o.max_y < self.min_y)
    }
}

pub fn bbox(pts: &[Point]) -> BBox {
    let mut ax=f64::MAX; let mut ay=f64::MAX;
    let mut bx=f64::MIN; let mut by=f64::MIN;
    for p in pts {
        if p[0]<ax{ax=p[0];} if p[1]<ay{ay=p[1];}
        if p[0]>bx{bx=p[0];} if p[1]>by{by=p[1];}
    }
    BBox{min_x:ax, min_y:ay, max_x:bx, max_y:by}
}

pub fn norm(pts: &[Point]) -> Poly {
    let b=bbox(pts); pts.iter().map(|p|[p[0]-b.min_x, p[1]-b.min_y]).collect()
}
pub fn center(pts: &[Point]) -> Poly {
    let b=bbox(pts); let cx=(b.min_x+b.max_x)/2.0; let cy=(b.min_y+b.max_y)/2.0;
    pts.iter().map(|p|[p[0]-cx, p[1]-cy]).collect()
}
pub fn rotate(pts: &[Point], deg: f64) -> Poly {
    let r=deg.to_radians(); let c=r.cos(); let s=r.sin();
    pts.iter().map(|p|[p[0]*c-p[1]*s, p[0]*s+p[1]*c]).collect()
}
pub fn translate(pts: &[Point], tx: f64, ty: f64) -> Poly {
    pts.iter().map(|p|[p[0]+tx, p[1]+ty]).collect()
}
pub fn ensure_cw(pts: &[Point]) -> Poly {
    if signed_area(pts)>0.0 { pts.iter().rev().cloned().collect() } else { pts.to_vec() }
}
pub fn ensure_ccw(pts: &[Point]) -> Poly {
    if signed_area(pts)<0.0 { pts.iter().rev().cloned().collect() } else { pts.to_vec() }
}

pub fn seg_intersect(p1:Point, p2:Point, p3:Point, p4:Point) -> bool {
    let d1=[p2[0]-p1[0], p2[1]-p1[1]]; let d2=[p4[0]-p3[0], p4[1]-p3[1]];
    let x=d1[0]*d2[1]-d1[1]*d2[0]; if x.abs()<TOL{return false;}
    let t=((p3[0]-p1[0])*d2[1]-(p3[1]-p1[1])*d2[0])/x;
    let u=((p3[0]-p1[0])*d1[1]-(p3[1]-p1[1])*d1[0])/x;
    t>TOL && t<1.0-TOL && u>TOL && u<1.0-TOL
}

pub fn pt_in_poly(px: f64, py: f64, pts: &[Point]) -> bool {
    let n=pts.len(); let mut inside=false; let mut j=n-1;
    for i in 0..n {
        let xi=pts[i][0]; let yi=pts[i][1];
        let xj=pts[j][0]; let yj=pts[j][1];
        if ((yi>py)!=(yj>py)) && (px<(xj-xi)*(py-yi)/(yj-yi)+xi) { inside=!inside; }
        j=i;
    }
    inside
}

pub fn polys_overlap(a: &[Point], b: &[Point]) -> bool {
    if a.is_empty()||b.is_empty(){return false;}
    let ba=bbox(a); let bb2=bbox(b);
    if !ba.overlaps(&bb2){return false;}
    for p in a { if pt_in_poly(p[0],p[1],b){return true;} }
    for p in b { if pt_in_poly(p[0],p[1],a){return true;} }
    let na=a.len(); let nb=b.len();
    for i in 0..na {
        for j in 0..nb {
            if seg_intersect(a[i],a[(i+1)%na],b[j],b[(j+1)%nb]){return true;}
        }
    }
    false
}

pub fn offset_poly(pts: &[Point], d: f64) -> Poly {
    if d<=0.0||pts.len()<3{return pts.to_vec();}
    let cx:f64=pts.iter().map(|p|p[0]).sum::<f64>()/pts.len()as f64;
    let cy:f64=pts.iter().map(|p|p[1]).sum::<f64>()/pts.len()as f64;
    pts.iter().map(|p|{
        let dx=p[0]-cx; let dy=p[1]-cy;
        let l=(dx*dx+dy*dy).sqrt().max(1e-9);
        [p[0]+(dx/l)*d, p[1]+(dy/l)*d]
    }).collect()
}

pub fn convex_hull(pts: &[Point]) -> Poly {
    let n=pts.len(); if n<3{return pts.to_vec();}
    let mut s=pts.to_vec();
    s.sort_by(|a,b|a[0].partial_cmp(&b[0]).unwrap().then(a[1].partial_cmp(&b[1]).unwrap()));
    let cross=|o:&Point,a:&Point,b:&Point|(a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);
    let mut lo:Vec<Point>=Vec::new(); let mut hi:Vec<Point>=Vec::new();
    for p in &s { while lo.len()>=2&&cross(&lo[lo.len()-2],&lo[lo.len()-1],p)<=0.0{lo.pop();} lo.push(*p); }
    for p in s.iter().rev() { while hi.len()>=2&&cross(&hi[hi.len()-2],&hi[hi.len()-1],p)<=0.0{hi.pop();} hi.push(*p); }
    hi.pop(); lo.pop(); lo.extend(hi); lo
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Xform { pub cx:f64, pub cy:f64, pub cos_r:f64, pub sin_r:f64, pub off_x:f64, pub off_y:f64, pub deg:f64 }

pub fn build_xform(hull: &[Point], deg: f64, tx: f64, ty: f64) -> Xform {
    let b=bbox(hull); let cx=(b.min_x+b.max_x)/2.0; let cy=(b.min_y+b.max_y)/2.0;
    let r=deg.to_radians(); let cos_r=r.cos(); let sin_r=r.sin();
    let rc:Vec<Point>=hull.iter().map(|p|{let dx=p[0]-cx;let dy=p[1]-cy;[dx*cos_r-dy*sin_r,dx*sin_r+dy*cos_r]}).collect();
    let rb=bbox(&rc);
    Xform{cx, cy, cos_r, sin_r, off_x:tx-rb.min_x, off_y:ty-rb.min_y, deg}
}

pub fn apply_xform(p: Point, xf: &Xform) -> Point {
    let dx=p[0]-xf.cx; let dy=p[1]-xf.cy;
    [dx*xf.cos_r-dy*xf.sin_r+xf.off_x, dx*xf.sin_r+dy*xf.cos_r+xf.off_y]
}
