mod geometry;
mod nesting;

use actix_web::{web, App, HttpServer, HttpResponse};
use actix_cors::Cors;
use actix_files::Files;
use serde::{Deserialize, Serialize};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicBool, Ordering};
use tokio::sync::broadcast;
use async_stream::stream;
use nesting::{Piece, NestConfig, Placement, run_nesting, Cb};

#[derive(Deserialize)]
struct NestReq {
    queue:Vec<PieceReq>, sheet_w:f64, sheet_h:f64, margin:f64, spacing:f64,
    rotations:Vec<f64>, ga_pop:usize, ga_gens:usize, continuous:bool, max_sheets:usize,
    #[serde(default="default_strategy")] strategy:String,
}
fn default_strategy()->String{"left".to_string()}

#[derive(Deserialize)]
struct PieceReq{id:u32, cutting_hull:Vec<[f64;2]>, qi:usize}

#[derive(Serialize,Clone)]
#[serde(tag="type")]
enum Evt {
    Log      {msg:String, level:String},
    Progress {pct:f64},
    Best     {sheet_num:usize, placements:Vec<Placement>, fitness:f64, util:u32},
    SheetDone{sheet_num:usize, placements:Vec<Placement>, fitness:f64, util:u32},
    Done     {total_placed:usize, total:usize, sheet_count:usize},
}

struct State{stop:Arc<AtomicBool>, running:Mutex<bool>}

async fn api_nest(body:web::Json<NestReq>, st:web::Data<State>) -> HttpResponse {
    {
        let mut r=st.running.lock().unwrap();
        if *r{st.stop.store(true,Ordering::Relaxed);tokio::time::sleep(std::time::Duration::from_millis(200)).await;}
        *r=true;
    }
    st.stop.store(false,Ordering::Relaxed);
    let req=body.into_inner(); let total=req.queue.len();
    let queue:Vec<Piece>=req.queue.iter().map(|p|Piece{id:p.id,cutting_hull:p.cutting_hull.clone(),qi:p.qi}).collect();
    let cfg=NestConfig{
        sheet_w:req.sheet_w, sheet_h:req.sheet_h, margin:req.margin, spacing:req.spacing,
        rotations:req.rotations.clone(), ga_pop:req.ga_pop.max(4), ga_gens:req.ga_gens.max(1),
        continuous:req.continuous, max_sheets:req.max_sheets, strategy:req.strategy.clone(),
    };
    let stop=Arc::clone(&st.stop);
    let (tx,_)=broadcast::channel::<String>(512);
    let tx=Arc::new(tx); let tx2=Arc::clone(&tx); let st2=web::Data::clone(&st);

    tokio::task::spawn_blocking(move||{
        let send=|e:Evt|{let j=serde_json::to_string(&e).unwrap_or_default();let _=tx2.send(format!("data:{}\n\n",j));};
        send(Evt::Log{msg:format!("Rust — {} pièces — pop:{} gens:{} stratégie:{}",total,cfg.ga_pop,cfg.ga_gens,cfg.strategy),level:"info".into()});
        let tx3=Arc::clone(&tx2);
        let cb:Cb=Arc::new(move|sheet_num,placements,fitness,util|{
            let e=Evt::Best{sheet_num,placements,fitness,util};
            let j=serde_json::to_string(&e).unwrap_or_default();
            let _=tx3.send(format!("data:{}\n\n",j));
        });
        let results=run_nesting(queue,cfg,stop,cb);
        let tp:usize=results.iter().map(|r|r.placements.len()).sum();
        for r in &results{
            let e=Evt::SheetDone{sheet_num:r.sheet_num,placements:r.placements.clone(),fitness:r.fitness,util:r.util_pct};
            let j=serde_json::to_string(&e).unwrap_or_default();
            let _=tx2.send(format!("data:{}\n\n",j));
        }
        let _=tx2.send(format!("data:{}\n\n",serde_json::to_string(&Evt::Done{total_placed:tp,total,sheet_count:results.len()}).unwrap_or_default()));
        let mut r=st2.running.lock().unwrap(); *r=false;
    });

    let mut rx=tx.subscribe();
    let s=stream!{
        yield Ok::<_,actix_web::Error>(web::Bytes::from("retry:1000\n\n"));
        loop{match rx.recv().await{
            Ok(m)=>{let done=m.contains("\"Done\"");yield Ok(web::Bytes::from(m));if done{break;}}
            Err(_)=>break,
        }}
    };
    HttpResponse::Ok().content_type("text/event-stream")
        .insert_header(("Cache-Control","no-cache"))
        .insert_header(("X-Accel-Buffering","no"))
        .streaming(s)
}

async fn api_stop(st:web::Data<State>)->HttpResponse{
    st.stop.store(true,Ordering::Relaxed);
    HttpResponse::Ok().json(serde_json::json!({"ok":true}))
}

async fn api_status(st:web::Data<State>)->HttpResponse{
    let r=*st.running.lock().unwrap();
    HttpResponse::Ok().json(serde_json::json!({"running":r,"cpus":rayon::current_num_threads(),"engine":"NestForge Rust 1.0"}))
}

#[actix_web::main]
async fn main()->std::io::Result<()>{
    let port=std::env::var("PORT").unwrap_or("5000".into());
    let addr=format!("127.0.0.1:{}",port);
    let state=web::Data::new(State{stop:Arc::new(AtomicBool::new(false)),running:Mutex::new(false)});
    println!("╔══════════════════════════════════════╗");
    println!("║  NestForge Engine v1.0  (Rust)       ║");
    println!("║  http://localhost:{}               ║",port);
    println!("║  CPU threads: {}                    ║",rayon::current_num_threads());
    println!("╚══════════════════════════════════════╝");
    #[cfg(windows)] std::process::Command::new("cmd").args(["/c","start",&format!("http://localhost:{}",port)]).spawn().ok();
    #[cfg(target_os="macos")] std::process::Command::new("open").arg(format!("http://localhost:{}",port)).spawn().ok();
    #[cfg(target_os="linux")] std::process::Command::new("xdg-open").arg(format!("http://localhost:{}",port)).spawn().ok();
    HttpServer::new(move||{
        App::new()
            .app_data(state.clone())
            .app_data(web::JsonConfig::default().limit(64*1024*1024))
            .wrap(Cors::default().allow_any_origin().allow_any_method().allow_any_header())
            .route("/api/nest",   web::post().to(api_nest))
            .route("/api/stop",   web::get().to(api_stop))
            .route("/api/status", web::get().to(api_status))
            .service(Files::new("/",".").index_file("nestforge.html"))
    }).bind(&addr)?.run().await
}
