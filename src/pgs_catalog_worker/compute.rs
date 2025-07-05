use gloo_worker::HandlerId;
use ordered_float::NotNan;
use std::collections::BTreeMap;

use analysis::util::Stats;

use crate::util::yield_now;

use super::{
    globals::scores, OutputGetEditAnalysis, OutputGetScores, VariantInfo, VariantSampleInfo,
};

/// `async` to give it a chance to be interrupted.
pub async fn compute_scores(id: HandlerId, edit_count: isize) -> Option<OutputGetScores> {
    log::info!("[Worker][{id:?}] Computing scores.");

    let values = compute_scores_(id, edit_count).await?;

    let iter = || values.iter().copied();
    let mean = analysis::util::mean(iter()).unwrap();
    let std_dev = analysis::util::std_dev(mean, iter()).unwrap();
    let iter = || values.iter().copied().map(|v| NotNan::new(v).unwrap());
    let min = *iter().min().unwrap();
    let max = *iter().max().unwrap();

    Some(OutputGetScores {
        scores: values,
        stats: Stats {
            mean,
            std_dev,
            min,
            max,
        },
    })
}
/// `async` to give it a chance to be interrupted.
async fn compute_scores_(id: HandlerId, edit_count: isize) -> Option<Vec<f64>> {
    let scores = scores(id).await;

    let mut data = Vec::new();
    for (i, s) in scores.scores().enumerate() {
        data.push(s.edited_score(edit_count)?);

        if i % 400 == 0 {
            yield_now(id, Some(i)).await;
        }
    }
    Some(data)
}

/// `async` to give it a chance to be interrupted.
pub async fn compute_edit_analysis(
    id: HandlerId,
    edit_count: isize,
) -> Option<OutputGetEditAnalysis> {
    log::info!("[Worker][{id:?}] Computing edit analysis.");

    let scores = scores(id).await;

    let decreasing = edit_count < 0;
    let edit_count = edit_count.unsigned_abs();

    let mut edits = BTreeMap::new();
    for (i, score) in scores.scores().enumerate() {
        let all_edits = match decreasing {
            true => &score.worst,
            false => &score.best,
        };
        let all_edits = all_edits.get(..edit_count)?;
        for variant in all_edits.iter() {
            let key = variant.key;
            let edit = edits
                .entry(key)
                .or_insert_with(|| VariantInfo::new(key, scores));

            edit.count += 1;
            edit.ploidy_dosage.increment(VariantSampleInfo {
                dosage: variant.dosage,
                ploidy: variant.ploidy,
            });
        }

        if i % 400 == 0 {
            yield_now(id, Some(i)).await;
        }
    }

    Some(OutputGetEditAnalysis {
        edits: edits.into_values().collect(),
    })
}
