mod compute;
mod globals;

pub mod files;

use gloo_worker::{HandlerId, Worker, WorkerScope};
use ordered_float::NotNan;
use pgs_catalog::{metadata::Metadata, PgsId};
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, ops::Range};

use utile::collections::counting_set::CountingBTreeSet;

use analysis::{
    pgs_scores::{Association, PgsCatalogAssociationKey, PgsCatalogScores},
    util::Stats,
};

use crate::util::{spawn_task, AsJson};

use self::globals::*;

pub struct PgsWorkerStruct {}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub enum Input {
    Set { origin: String, pgs_id: PgsId },
    GetScores { edit_count: isize, normalise: bool },
    GetEditAnalysis { edit_count: isize, normalise: bool },
}

#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)]
pub enum Output {
    Init(OutputInit),
    GetScores(Option<OutputGetScores>),
    GetEditAnalysis(Option<OutputGetEditAnalysis>),
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct OutputInit {
    pub metadata: AsJson<Metadata>,

    pub sample_count: usize,

    pub stats: Option<Stats>,
    pub full_score_range: Option<Range<f64>>,
    pub norm_stats: Option<Stats>,
    pub norm_full_score_range: Option<Range<f64>>,
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
    pub edits: Vec<VariantInfo>,
}
impl OutputGetEditAnalysis {
    pub fn normalised(self, std_dev: NotNan<f64>) -> Self {
        Self {
            edits: self
                .edits
                .into_iter()
                .map(|v| v.normalised(std_dev))
                .collect(),
        }
    }
}
#[derive(Debug, Clone)]
#[derive(Serialize, Deserialize)]
pub struct VariantInfo {
    pub count: usize,
    pub association: AsJson<Association>,
    pub ploidy_dosage: CountingBTreeSet<VariantSampleInfo>,
}
impl VariantInfo {
    pub fn new(key: PgsCatalogAssociationKey, scores: &PgsCatalogScores) -> Self {
        let association = scores.associations.get(&key).unwrap().clone();

        VariantInfo {
            count: 0,
            association: AsJson(association),
            ploidy_dosage: CountingBTreeSet::new(),
        }
    }
    pub fn normalised(self, std_dev: NotNan<f64>) -> Self {
        Self {
            count: self.count,
            association: AsJson(self.association.0.normalised(std_dev)),
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

impl Worker for PgsWorkerStruct {
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
            Input::Set { origin, pgs_id } => {
                let scope = scope.clone();
                let respond = move |metadata: Metadata, scores: &PgsCatalogScores| {
                    let stats = scores.stats();

                    let full_score_range = scores.full_score_range();

                    let norm_stats = stats.map(|s| s.normalised_self());

                    let norm_full_score_range = stats.map(|s| {
                        let full_score_range = full_score_range.clone().unwrap();
                        s.normalise_value(full_score_range.start)
                            ..s.normalise_value(full_score_range.end)
                    });

                    let response = Output::Init(OutputInit {
                        metadata: AsJson(metadata),

                        sample_count: scores.scores.len(),

                        stats,

                        full_score_range,

                        norm_stats,

                        norm_full_score_range,
                    });

                    scope.respond(id, response);
                };
                wasm_bindgen_futures::spawn_local(async move {
                    let metadata = files::fetch_all_metadata(&origin, Some(pgs_id))
                        .await
                        .unwrap();
                    let scores = files::fetch_scores(id, &origin, pgs_id).await.unwrap();
                    SCORES.set(scores.clone()).unwrap();

                    respond(metadata, SCORES.get().unwrap());
                });
            }
            Input::GetScores {
                edit_count,
                normalise,
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data = compute::compute_scores(id, edit_count).await;
                    if normalise && let Some(s) = stats(id).await {
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
            } => {
                let scope = scope.clone();
                let f = async move {
                    let mut data = compute::compute_edit_analysis(id, edit_count).await;
                    if normalise && let Some(s) = stats(id).await {
                        let std_dev = NotNan::new(s.std_dev).unwrap();
                        data = data.map(|d| d.normalised(std_dev));
                    }
                    log::info!("[Worker][{id:?}] Sending edit analysis.");
                    scope.respond(id, Output::GetEditAnalysis(data))
                };
                spawn_task(&EDIT_ANALYSIS_TASK, id, f);
            }
        }
    }
}
