use futures::stream::AbortHandle;
use gloo_worker::HandlerId;
use std::{cell::RefCell, sync::OnceLock};

use analysis::{pgs_scores::PgsCatalogScores, util::Stats};

use crate::util::{wait_for, HandlerIdOrd};

thread_local! {
    pub static SCORES_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
    pub static EDIT_ANALYSIS_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
}

pub static SCORES: OnceLock<PgsCatalogScores> = OnceLock::new();

pub async fn scores(id: HandlerId) -> &'static PgsCatalogScores {
    wait_for(id, || SCORES.get()).await
}

pub async fn stats(id: HandlerId) -> Option<Stats> {
    static STATS: OnceLock<Option<Stats>> = OnceLock::new();

    if let Some(s) = STATS.get() {
        return *s;
    }

    let scores = scores(id).await;
    let s = scores.stats();
    STATS.set(s).unwrap();
    s
}
