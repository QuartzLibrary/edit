use futures::stream::AbortHandle;
use gloo_worker::HandlerId;
use std::{cell::RefCell, sync::OnceLock};

use analysis::{
    pvalues::{PhenotypeTopPValues, TopPValueVariant},
    scores::Scores,
    util::Stats,
};

use crate::util::{wait_for, HandlerIdOrd};

thread_local! {
    pub static SCORES_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
    pub static EDIT_ANALYSIS_TASK: RefCell<Option<(HandlerIdOrd, AbortHandle)>> = const { RefCell::new(None) };
}

pub static SCORES: OnceLock<Scores> = OnceLock::new();
pub static PVALUES: OnceLock<PhenotypeTopPValues> = OnceLock::new();

pub async fn scores(id: HandlerId) -> &'static Scores {
    wait_for(id, || SCORES.get()).await
}
pub async fn pvalues(id: HandlerId) -> &'static PhenotypeTopPValues {
    wait_for(id, || PVALUES.get()).await
}

pub async fn stats(id: HandlerId, use_hq: bool) -> Option<Stats> {
    static STATS: OnceLock<Option<Stats>> = OnceLock::new();
    static STATS_HQ: OnceLock<Option<Stats>> = OnceLock::new();

    let store = if use_hq { &STATS_HQ } else { &STATS };

    if let Some(s) = store.get() {
        return *s;
    }

    let scores = scores(id).await;
    let s = if use_hq {
        scores.stats_hq()
    } else {
        scores.stats()
    };
    store.set(s).unwrap();
    s
}
/// The variants with the highest p-values.
pub async fn top_pvalue_variants(
    id: HandlerId,
    use_hq: bool,
) -> &'static [&'static TopPValueVariant] {
    fn get_top_variants_hq(pvalues: &PhenotypeTopPValues) -> Vec<&TopPValueVariant> {
        let mut top_variants: Vec<_> = pvalues
            .top_variants
            .values()
            .filter(|v| v.stat.neglog10_pval_meta_hq.is_some() && v.stat.beta_meta_hq.is_some())
            .collect();
        top_variants.sort_unstable_by_key(|v| -v.stat.neglog10_pval_meta_hq.unwrap());
        top_variants
    }
    fn get_top_variants(pvalues: &PhenotypeTopPValues) -> Vec<&TopPValueVariant> {
        let mut top_variants: Vec<_> = pvalues
            .top_variants
            .values()
            .filter(|v| v.stat.neglog10_pval_meta.is_some() && v.stat.beta_meta.is_some())
            .collect();
        top_variants.sort_unstable_by_key(|v| -v.stat.neglog10_pval_meta.unwrap());
        top_variants
    }

    static TOP: OnceLock<Vec<&'static TopPValueVariant>> = OnceLock::new();
    static TOP_HQ: OnceLock<Vec<&'static TopPValueVariant>> = OnceLock::new();

    let store = if use_hq { &TOP_HQ } else { &TOP };

    if let Some(top) = store.get() {
        return top.as_slice();
    }

    log::info!("[Worker][{id:?}] Getting top p-values (use_hq: {use_hq}).");
    let start = web_time::Instant::now();

    let pvalues = pvalues(id).await;
    let top = if use_hq {
        get_top_variants_hq(pvalues)
    } else {
        get_top_variants(pvalues)
    };
    store.set(top).unwrap();

    log::info!(
        "[Worker][{id:?}] Got top p-values in {:?}ms",
        start.elapsed().as_millis()
    );

    store.get().unwrap().as_slice()
}
