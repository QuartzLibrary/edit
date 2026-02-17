use gloo_worker::HandlerId;
use std::io;

use pan_ukbb::PhenotypeManifestEntry;
use resource::UrlResource;

use analysis::{pvalues::PhenotypeTopPValues, scores::Scores};

use crate::util::load_large_json;

pub async fn fetch_scores(id: HandlerId, origin: &str, file: String) -> io::Result<Scores> {
    let url = format!("{origin}/data/pan_ukbb/scores/{file}.json.br");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut scores: Scores = load_large_json(&url).await?;

    // "Based on these thresholds, two individuals, NA21310 and HG02300,
    // were listed as males, but had genotypes consistent with females"
    // https://www.biorxiv.org/content/10.1101/078600v1.full.pdf
    scores.scores.remove("HG02300");
    scores.scores.remove("NA21310");

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(scores)
}
pub async fn fetch_pvalues(
    id: HandlerId,
    origin: &str,
    file: String,
) -> io::Result<PhenotypeTopPValues> {
    let url = format!("{origin}/data/pan_ukbb/pvalues/{file}.json.br");

    log::info!("[Worker][{id:?}] Loading initial data from URL: {url}");

    let mut pvalues: PhenotypeTopPValues = load_large_json(&url).await?;

    // "Based on these thresholds, two individuals, NA21310 and HG02300,
    // were listed as males, but had genotypes consistent with females"
    // https://www.biorxiv.org/content/10.1101/078600v1.full.pdf
    pvalues.samples.remove("HG02300");
    pvalues.samples.remove("NA21310");

    for variant in pvalues.top_variants.values_mut() {
        variant.genotypes.remove("HG02300");
        variant.genotypes.remove("NA21310");
    }

    log::info!("[Worker][{id:?}] Data loaded from {url}");

    Ok(pvalues)
}
pub async fn fetch_manifest(origin: String) -> io::Result<Vec<PhenotypeManifestEntry>> {
    let url = format!("{origin}/data/pan_ukbb/phenotype_manifest.tsv");

    log::info!("[App] Fetching manifest from: {url}");

    let entries = PhenotypeManifestEntry::load_async(UrlResource::new(&url)?).await?;

    log::info!("[App] Manifest loaded.");

    Ok(entries)
}
