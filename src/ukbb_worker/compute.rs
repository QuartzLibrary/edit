use futures::StreamExt;
use gloo_worker::HandlerId;
use ordered_float::NotNan;
use std::collections::BTreeMap;

use hail::contig::GRCh38Contig;
use pan_ukbb::SummaryStats;
use utile::collections::counting_set::CountingBTreeSet;

use analysis::{
    pvalues::{SampleGenotype, TopPValueVariant},
    util::{Stats, SummaryStatKey},
};

use crate::util::{yield_now, AsJson};

use super::{
    globals::{scores, top_pvalue_variants},
    OutputGetEditAnalysis, OutputGetScores, VariantInfo, VariantSampleInfo,
};

/// `async` to give it a chance to be interrupted.
pub async fn compute_scores(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<OutputGetScores> {
    log::info!("[Worker][{id:?}] Computing scores.");

    let values = compute_scores_(id, edit_count, use_hq, top_pvalues).await?;

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
async fn compute_scores_(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<Vec<f64>> {
    let scores = scores(id).await;

    if top_pvalues == 0 || edit_count == 0 {
        let mut data = Vec::new();
        for (i, s) in scores.scores().enumerate() {
            data.push(if use_hq {
                s.edited_score_hq(edit_count)?
            } else {
                s.edited_score(edit_count)?
            });

            if i % 400 == 0 {
                yield_now(id, Some(i)).await;
            }
        }
        Some(data)
    } else {
        let mut top_variants = top_pvalue_variants(id, use_hq).await;
        if top_variants.len() < top_pvalues {
            return None;
        }
        top_variants = &top_variants[0..top_pvalues];
        if top_variants.len() < edit_count.unsigned_abs() {
            return None;
        }

        yield_now(id, None).await;

        let edits = scores.scores.keys().map(|id| {
            get_top_edits(top_variants, id, use_hq, edit_count)
                .into_iter()
                .map(|(s, g)| actual_edit(s, g, use_hq, edit_count).unwrap())
                .sum::<NotNan<f64>>()
        });

        let base_scores = Box::pin(compute_scores_(id, 0, use_hq, 0)).await?;

        Some(
            futures::stream::iter(base_scores.into_iter().zip(edits))
                .enumerate()
                .then(|(i, (base, edit))| async move {
                    if i % 400 == 0 {
                        yield_now(id, Some(i)).await;
                    }
                    base + *edit
                })
                .collect()
                .await,
        )
    }
}

/// `async` to give it a chance to be interrupted.
pub async fn compute_edit_analysis(
    id: HandlerId,
    edit_count: isize,
    use_hq: bool,
    top_pvalues: usize,
) -> Option<OutputGetEditAnalysis> {
    log::info!("[Worker][{id:?}] Computing edit analysis.");

    let scores = scores(id).await;

    if top_pvalues == 0 || edit_count == 0 {
        let decreasing = edit_count < 0;
        let edit_count = edit_count.unsigned_abs();

        let mut edits = BTreeMap::new();
        for (i, score) in scores.scores().enumerate() {
            let all_edits = match (use_hq, decreasing) {
                (true, true) => &score.worst_hq,
                (true, false) => &score.best_hq,
                (false, true) => &score.worst,
                (false, false) => &score.best,
            };
            let all_edits = all_edits.get(..edit_count)?;
            for variant in all_edits.iter() {
                let key = variant.key.clone().key;
                let edit = edits
                    .entry(key.clone())
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
            use_hq,
            edits: edits.into_values().collect(),
        })
    } else {
        let mut top_variants = top_pvalue_variants(id, use_hq).await;
        if top_variants.len() < top_pvalues {
            return None;
        }
        top_variants = &top_variants[0..top_pvalues];
        if top_variants.len() < edit_count.unsigned_abs() {
            return None;
        }

        yield_now(id, None).await;

        let mut edits = BTreeMap::new();
        for (i, sample) in scores.scores.keys().enumerate() {
            for (s, g) in get_top_edits(top_variants, sample, use_hq, edit_count) {
                let key = SummaryStatKey::new(s);
                let edit = edits.entry(key).or_insert_with(|| VariantInfo {
                    count: 0,
                    stat: AsJson(s.clone()),
                    ploidy_dosage: CountingBTreeSet::new(),
                });
                edit.count += 1;
                edit.ploidy_dosage.increment(VariantSampleInfo {
                    dosage: g.dosage,
                    ploidy: g.ploidy,
                });
            }
            if i % 400 == 0 {
                yield_now(id, Some(i)).await;
            }
        }

        Some(OutputGetEditAnalysis {
            use_hq,
            edits: edits.into_values().collect(),
        })
    }
}

/// The variants with the largest effects in the slice.
fn get_top_edits<'a>(
    top_variants: &'a [&TopPValueVariant],
    id: &str,
    use_hq: bool,
    edit_count: isize,
) -> Vec<(&'a SummaryStats<GRCh38Contig>, SampleGenotype)> {
    let positive_edits = edit_count >= 0;

    let mut edits: Vec<_> = top_variants
        .iter()
        .map(|v| (&v.stat, v.genotypes[id]))
        .collect();

    edits.sort_unstable_by_key(|(s, g)| {
        let edit = actual_edit(s, *g, use_hq, edit_count).unwrap();
        if positive_edits {
            -edit
        } else {
            edit
        }
    });

    edits.truncate(edit_count.unsigned_abs());

    edits
}

fn actual_edit(
    s: &SummaryStats<GRCh38Contig>,
    g: SampleGenotype,
    use_hq: bool,
    edit_count: isize,
) -> Option<NotNan<f64>> {
    match (use_hq, edit_count >= 0) {
        (true, true) => s.max_edit_hq(g.dosage, g.ploidy),
        (true, false) => s.min_edit_hq(g.dosage, g.ploidy),
        (false, true) => s.max_edit(g.dosage, g.ploidy),
        (false, false) => s.min_edit(g.dosage, g.ploidy),
    }
}
