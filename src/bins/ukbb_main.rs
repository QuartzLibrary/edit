use gloo_worker::Registrable;
use leptos::task::Executor;

use edit::ukbb_worker::WorkerStruct;

fn main() {
    console_error_panic_hook::set_once();

    console_log::init().unwrap();

    _ = Executor::init_wasm_bindgen();

    WorkerStruct::registrar().register();

    log::info!("[Worker] Init");
}
