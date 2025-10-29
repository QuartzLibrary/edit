mod compute;
mod globals;

pub mod files;

use gloo_worker::{HandlerId, Worker, WorkerScope};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, ops::Range, time::Duration};

use hail::contig::GRCh38Contig;
use pan_ukbb::{PhenotypeManifestEntry, SummaryStats};
use utile::collections::counting_set::CountingBTreeSet;

use analysis::{
    scores::Scores,
    util::{Stats, SummaryStatKey},
};

use crate::util::{AsJson, spawn_task};

use self::globals::*;

pub struct WorkerStruct {}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub enum Input {
    Set {
        origin: String,
        file: String,
    },
    GetScores {
        edit_count: isize,
        normalise: bool,
        use_hq: bool,
        top_pvalues: usize,
    },
    GetEditAnalysis {
        edit_count: isize,
        normalise: bool,
        use_hq: bool,
        top_pvalues: usize,
    },
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
#[expect(clippy::large_enum_variant)]
pub enum Output {
    Init(OutputInit),
    GetScores(Option<OutputGetScores>),
    GetEditAnalysis(Option<OutputGetEditAnalysis>),
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputInit {
    // Horrible hack (AsJson)
    pub phenotype: AsJson<PhenotypeManifestEntry>,

    pub sample_count: usize,

    pub stats: Option<Stats>,
    pub stats_hq: Option<Stats>,

    pub full_score_range: Option<Range<f64>>,
    pub full_score_range_hq: Option<Range<f64>>,

    pub norm_stats: Option<Stats>,
    pub norm_stats_hq: Option<Stats>,

    pub norm_full_score_range: Option<Range<f64>>,
    pub norm_full_score_range_hq: Option<Range<f64>>,
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputGetScores {
    pub scores: Vec<f64>,

    pub stats: Stats,
}
impl OutputGetScores {
    pub fn normalised(self, mean: f64, std_dev: f64) -> Self {
        Self {
            scores: self
                .scores
                .into_iter()
                .map(|v| (v - mean) / std_dev)
                .collect(),
            stats: self.stats.normalised(mean, std_dev),
        }
    }
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputGetEditAnalysis {
    pub use_hq: bool,
    /// Pre-sorted by p-value.
    pub edits: Vec<VariantInfo>,
}
impl OutputGetEditAnalysis {
    pub fn normalised(self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            use_hq: self.use_hq,
            edits: self
                .edits
                .into_iter()
                .map(|v| v.normalised(std_dev, std_dev_hq))
                .collect(),
        }
    }
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct VariantInfo {
    pub count: usize,
    pub stat: AsJson<SummaryStats<GRCh38Contig>>,
    pub ploidy_dosage: CountingBTreeSet<VariantSampleInfo>,
}
impl VariantInfo {
    pub fn new(key: SummaryStatKey, scores: &Scores) -> Self {
        let SummaryStatKey {
            chr,
            pos,
            ref_allele,
            alt,
            ..
        } = key;

        let stat = scores
            .summary_stats
            .iter()
            .find(|s| s.chr == chr && s.pos == pos && s.ref_allele == ref_allele && s.alt == alt)
            .unwrap()
            .clone();

        VariantInfo {
            count: 0,
            stat: AsJson(stat),
            ploidy_dosage: CountingBTreeSet::new(),
        }
    }
    pub fn normalised(self, std_dev: NotNan<f64>, std_dev_hq: NotNan<f64>) -> Self {
        Self {
            count: self.count,
            stat: AsJson(self.stat.0.normalised(std_dev, std_dev_hq)),
            ploidy_dosage: self.ploidy_dosage,
        }
    }
}
#[derive(Debug, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct VariantSampleInfo {
    pub dosage: u8,
    pub ploidy: u8,
}

impl Worker for WorkerStruct {
    type Input = Input;
    type Message = ();
    type Output = Output;

    fn create(_scope: &WorkerScope<Self>) -> Self {
        Self {}
    }

    fn update(&mut self, _scope: &WorkerScope<Self>, _msg: Self::Message) {
        unreachable!()
    }

    fn received(&mut self, scope: &WorkerScope<Self>, msg: Self::Input, id: HandlerId) {
        log::info!("[Worker][{id:?}] Received: {msg:?}");
        match msg {
            Input::Set { file, origin } => {
                let scope = scope.clone();
                let respond = move |scores: &Scores| {
                    let stats = scores.stats();
                    let stats_hq = scores.stats_hq();

                    let full_score_range = scores.full_score_range();
                    let full_score_range_hq = scores.full_score_range_hq();

                    let norm_stats = stats.map(|s| s.normalised_self());
                    let norm_stats_hq = stats_hq.map(|s| s.normalised_self());

                    let norm_full_score_range = stats.map(|s| {
                        let full_score_range = full_score_range.clone().unwrap();
                        s.normalise_value(full_score_range.start)
                            ..s.normalise_value(full_score_range.end)
                    });
                    let norm_full_score_range_hq = stats_hq.map(|s| {
                        let full_score_range_hq = full_score_range_hq.clone().unwrap();
                        s.normalise_value(full_score_range_hq.start)
                            ..s.normalise_value(full_score_range_hq.end)
                    });

                    let response = Output::Init(OutputInit {
                        phenotype: AsJson(scores.phenotype.clone()),

                        sample_count: scores.scores.len(),

                        stats,
                        stats_hq,

                        full_score_range,
                        full_score_range_hq,

                        norm_stats,
                        norm_stats_hq,

                        norm_full_score_range,
                        norm_full_score_range_hq,
                    });

                    scope.respond(id, response);
                };
                wasm_bindgen_futures::spawn_local(async move {
                    let scores = files::fetch_scores(id, &origin, file.clone())
                        .await
                        .unwrap();
                    SCORES.set(scores.clone()).unwrap();

                    respond(SCORES.get().unwrap());

                    // Give it a second for the initial scores to be sent.
                    utile::time::sleep(Duration::from_millis(1000)).await;

                    let pvalues = files::fetch_pvalues(id, &origin, file.clone())
                        .await
                        .unwrap();
                    PVALUES.set(pvalues).unwrap();
                });
            }
            Input::GetScores {
                edit_count,
                normalise,
                use_hq,
                top_pvalues,
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data =
                        compute::compute_scores(id, edit_count, use_hq, top_pvalues).await;
                    if normalise && let Some(s) = stats(id, use_hq).await {
                        data = data.map(|d| d.normalised(s.mean, s.std_dev));
                    }
                    log::info!("[Worker][{id:?}] Sending scores.");
                    scope.respond(id, Output::GetScores(data))
                };
                spawn_task(&SCORES_TASK, id, f);
            }
            Input::GetEditAnalysis {
                edit_count,
                normalise,
                use_hq,
                top_pvalues,
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data =
                        compute::compute_edit_analysis(id, edit_count, use_hq, top_pvalues).await;
                    if let Some(data) = &mut data {
                        data.edits.sort_unstable_by_key(|v| {
                            if use_hq {
                                -v.stat.neglog10_pval_meta_hq.unwrap()
                            } else {
                                -v.stat.neglog10_pval_meta.unwrap()
                            }
                        });
                    }
                    if normalise
                        && let Some(s) = stats(id, false).await
                        && let Some(s_hq) = stats(id, true).await
                    {
                        let std_dev = NotNan::new(s.std_dev).unwrap();
                        let std_dev_hq = NotNan::new(s_hq.std_dev).unwrap();
                        data = data.map(|d| d.normalised(std_dev, std_dev_hq));
                    }
                    log::info!("[Worker][{id:?}] Sending edit analysis.");
                    scope.respond(id, Output::GetEditAnalysis(data))
                };
                spawn_task(&EDIT_ANALYSIS_TASK, id, f);
            }
        }
    }
}
