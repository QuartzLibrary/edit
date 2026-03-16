use gloo_worker::Registrable;
use any_spawner::Executor;

use edit::pgs_worker::PgsWorkerStruct;

fn main() {
    console_error_panic_hook::set_once();

    console_log::init().unwrap();

    _ = Executor::init_wasm_bindgen();

    PgsWorkerStruct::registrar().register();

    log::info!("[Worker] Init");
}
