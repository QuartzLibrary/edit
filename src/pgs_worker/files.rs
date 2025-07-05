use gloo_worker::HandlerId;
use pgs_catalog::{metadata::Metadata, PgsId};
use std::io;

use analysis::pgs_scores::PgsCatalogScores;

use crate::util::load_large_json;

pub async fn fetch_scores(
    id: HandlerId,
    origin: &str,
    pgs_id: PgsId,
) -> io::Result<PgsCatalogScores> {
    let url = format!("{origin}/data/pgs_catalog/scores/{pgs_id}.json.br");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut scores: PgsCatalogScores = load_large_json(&url).await?;

    // "Based on these thresholds, two individuals, NA21310 and HG02300,
    // were listed as males, but had genotypes consistent with females"
    // https://www.biorxiv.org/content/10.1101/078600v1.full.pdf
    scores.scores.remove("HG02300");
    scores.scores.remove("NA21310");

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(scores)
}
pub async fn fetch_all_metadata(origin: &str, pgs_id: Option<PgsId>) -> io::Result<Metadata> {
    let url = if let Some(pgs_id) = pgs_id {
        format!("{origin}/data/pgs_catalog/metadata/{pgs_id}.json.br")
    } else {
        format!("{origin}/data/pgs_catalog/metadata/all.json.br")
    };

    log::info!("[App] Fetching metadata from: {url}");

    let metadata: Metadata = load_large_json(&url).await?;

    log::info!("[App] Metadata loaded.");

    Ok(metadata)
}
