use gloo_worker::{Spawnable, WorkerBridge};
use leptos::{IntoView, ev, html, prelude::*, view};
use leptos_ext::signal::{Load, ReadSignalExt, WriteSignalExt};
use ordered_float::NotNan;
use pgs_catalog::{
    PgsId,
    metadata::{Metadata, PerformanceMetric, Publication, ScoreDevelopmentSample},
    simplified::SimplifiedHarmonizedStudyAssociation,
};
use plotly::layout::{Shape, ShapeLine};
use std::{
    sync::{Arc, Mutex},
    time::Duration,
};
use utile::task::Task;

use analysis::{
    pgs_scores::{Association, parse_contig},
    util::Stats,
};

use edit::{
    pgs_worker::{
        Input, Output, OutputGetEditAnalysis, OutputGetScores, OutputInit, PgsWorkerStruct,
        VariantInfo, VariantSampleInfo,
    },
    send_option::SendOption,
    util::{AsJson, PLOTLY_THEME, render_edit_frequency_histogram, render_plotly},
};

const BINS: usize = 100;

#[derive(Debug, Clone)]
struct PageState {
    pgs_id: PgsId,

    bridge: Arc<SendOption<WorkerBridge<PgsWorkerStruct>>>,

    init: ArcRwSignal<Option<OutputInit>>,

    /// The number of edits to apply to the scores.
    /// Negative means reducing the scores.
    edit_count: RwSignal<isize>,
    /// Normalise the scores and edits.
    normalise: RwSignal<bool>,

    show_stats: RwSignal<bool>,
    full_range: RwSignal<bool>,
}

pub fn score_page(pgs_id: PgsId) -> impl IntoView {
    let state = PageState::new(pgs_id);

    let scores = state.scores_signal();
    let edits = state.edits_signal();

    let loading = scores.map(|o| !o.is_ready());
    let scores = scores.keep_if(|o| o.is_ready());

    move || {
        let PageState {
            pgs_id,
            bridge: _,
            init,
            edit_count,
            normalise,
            show_stats: _,
            full_range: _,
        } = state.clone();

        let init = init.get();

        html::div().class("page-content score-viewer").child((
            html::a().href("/").class("back-link").child("< Index"),
            html::h1().child(
                init.as_ref()
                    .map(|init| init.metadata.scores[0].name.clone())
                    .unwrap_or("Loading…".to_owned()),
            ),
            html::div()
                .class("card card-chart")
                .class(("loading", loading))
                .child({
                    let state = state.clone();
                    move || match scores.get() {
                        Load::Ready(Some(scores)) => {
                            render_main_chart(&state, scores.into()).into_any()
                        }
                        Load::Ready(None) => html::div()
                            .class("loading-text")
                            .child("Too many edits, no scores available.")
                            .into_any(),
                        Load::Loading => html::div()
                            .class("loading-text")
                            .child("Loading…")
                            .into_any(),
                    }
                }),
            render_controls(&state),
            init.as_ref().map(|init| {
                let init = init.clone();
                move || {
                    render_full_phenotype_info(
                        pgs_id,
                        &init,
                        scores.get(),
                        loading,
                        normalise.get(),
                    )
                }
            }),
            init.as_ref()
                .map(|init| render_edit_analysis_panel(init.clone(), edits, edit_count.into())),
            html::p()
                .class("footnote")
                .child("All genomic positions are on GRCh38/hg38 and 1-based."),
        ))
    }
}

fn render_main_chart(state: &PageState, data: Signal<OutputGetScores>) -> impl IntoView + use<> {
    const DEEP_SKY_BLUE: &str = "#00BFFF";
    const HOT_PINK: &str = "#FF69B4";
    const LIME_GREEN: &str = "#32CD32";
    const GOLD: &str = "#FFD700";

    const MEAN: &str = "mean";
    const MEAN_PLUS_STDDEV: &str = "+sd";
    const MEAN_MINUS_STDDEV: &str = "-sd";
    const MIN: &str = "min";
    const MAX: &str = "max";

    fn stat_line(x: f64, color: &'static str, name: &'static str) -> Shape {
        Shape::new()
            .x0(x)
            .x1(x)
            .y0(0.0)
            .y1(1.0)
            .y_ref("paper")
            .line(ShapeLine::new().color(color).width(2.0))
            .name(name)
    }
    fn unedited_stat_line(x: f64, color: &'static str, name: &'static str) -> Shape {
        Shape::new()
            .x0(x)
            .x1(x)
            .y0(0.0)
            .y1(1.0)
            .y_ref("paper")
            .line(ShapeLine::new().color(format!("{color}30")).width(1.0))
            .name(name)
    }

    fn create_stat_annotation(
        x: f64,
        text: &str,
        color: &'static str,
        opacity: f64,
    ) -> plotly::layout::Annotation {
        plotly::layout::Annotation::new()
            .x(x)
            .y(1.1)
            .y_ref("paper")
            .text(text)
            .show_arrow(false)
            .font(plotly::common::Font::new().color(color))
            .opacity(opacity)
    }

    let PageState {
        init,
        edit_count,
        normalise,
        show_stats,
        full_range,
        ..
    } = state.clone();

    let init = init
        .get()
        .expect("initial state should be ready before scores");

    let zero_edits = edit_count.is(0).dedup();

    let chart = Signal::derive_local(move || {
        // Untracked to avoid re-rendering before the data is ready. (Always trigger a change in the data.)
        let normalise = normalise.get_untracked();
        let zero_edits = zero_edits.get_untracked();

        let OutputGetScores { scores, stats } = data.get();
        let Stats {
            mean,
            std_dev,
            min,
            max,
        } = stats;

        let mut plot = plotly::Plot::new();
        plot.add_trace(plotly::Histogram::new(scores).n_bins_x(BINS).name("Scores"));

        let current_shapes = [
            stat_line(mean, DEEP_SKY_BLUE, MEAN),
            stat_line(mean + std_dev, HOT_PINK, MEAN_PLUS_STDDEV),
            stat_line(mean - std_dev, HOT_PINK, MEAN_MINUS_STDDEV),
            stat_line(min, LIME_GREEN, MIN),
            stat_line(max, GOLD, MAX),
        ];

        let unedited_shapes = match normalise {
            true => init.norm_stats,
            false => init.stats,
        }
        .map(|stats| {
            let Stats {
                mean,
                std_dev,
                min,
                max,
            } = stats;
            [
                unedited_stat_line(mean, DEEP_SKY_BLUE, MEAN),
                unedited_stat_line(mean + std_dev, HOT_PINK, MEAN_PLUS_STDDEV),
                unedited_stat_line(mean - std_dev, HOT_PINK, MEAN_MINUS_STDDEV),
                unedited_stat_line(min, LIME_GREEN, MIN),
                unedited_stat_line(max, GOLD, MAX),
            ]
        });

        let shapes = if show_stats.get() {
            let unedited_shapes = match unedited_shapes {
                Some(_) if zero_edits => None,
                Some(shapes) => Some(shapes),
                None => None,
            };

            current_shapes
                .into_iter()
                .chain(unedited_shapes.into_iter().flatten())
                .collect()
        } else {
            vec![]
        };

        // Create annotations for statistical lines
        let annotations = if show_stats.get() {
            vec![
                create_stat_annotation(mean, MEAN, DEEP_SKY_BLUE, 1.0),
                create_stat_annotation(mean + std_dev, MEAN_PLUS_STDDEV, HOT_PINK, 1.0),
                create_stat_annotation(mean - std_dev, MEAN_MINUS_STDDEV, HOT_PINK, 1.0),
                create_stat_annotation(min, MIN, LIME_GREEN, 1.0),
                create_stat_annotation(max, MAX, GOLD, 1.0),
            ]
        } else {
            vec![]
        };

        let mut layout = plotly::Layout::new()
            .template(PLOTLY_THEME.get())
            .shapes(shapes)
            .annotations(annotations)
            .margin(
                plotly::layout::Margin::new()
                    .top(60)
                    .right(60)
                    .bottom(60)
                    .left(60),
            );

        if full_range.get() {
            let range = match normalise {
                true => init.norm_full_score_range.clone(),
                false => init.full_score_range.clone(),
            };

            if let Some(range) = range {
                layout = layout.x_axis(
                    plotly::layout::Axis::new()
                        .range(vec![range.start, range.end])
                        .auto_range(false),
                );
            }
        }

        plot.set_layout(layout);

        plot
    });

    render_plotly(chart)
}

fn render_controls(state: &PageState) -> impl IntoView + use<> {
    const EDIT_COUNT_TIP: &str = "The number of edits to apply to the samples.";

    let PageState {
        init,
        edit_count,
        normalise,
        show_stats,
        full_range,
        ..
    } = state.clone();

    let buttons = [
        (
            normalise,
            "Normalise",
            "Normalise the scores to the range of the original data.",
        ),
        (
            show_stats,
            "Show stats",
            "Show some summary statistics in the chart above.",
        ),
        (
            full_range,
            "Full range",
            "Force the chart to show the full range of edits (makes it easert to eyeball the shift in the distribution).",
        ),
    ];

    let edit_limit = init
        .get()
        .map(|init| init.metadata.scores[0].number_of_variants)
        .unwrap_or(500);

    html::div().class("card controls").child((
        slider_with_controls(
            "Edits: ",
            edit_count.double_bind(|i| *i as f64, move |edit| *edit as isize),
            -(edit_limit as f64),
            edit_limit as f64,
            1.0,
            "0",
            Some(EDIT_COUNT_TIP),
        ),
        html::div()
            .class("controls-row")
            .child(buttons.map(|(signal, label, tip)| toggle_button(signal, label, Some(tip)))),
    ))
}

fn slider_with_controls(
    name: &'static str,
    signal: RwSignal<f64>,
    min: f64,
    max: f64,
    step: f64,
    zero: &'static str,
    tip: Option<&'static str>,
) -> impl IntoView {
    let slider = view! { <thaw::Slider value={signal} min={min} max={max} step={step} show_stops=false style="width: 100%;" />};
    let plus = html::button()
        .class("btn")
        .attr("disabled", signal.map(move |v| v == &max))
        .on(ev::click, move |_| {
            signal.update(|c| {
                if *c < max - step {
                    *c += step
                }
            })
        })
        .child("+");
    let minus = html::button()
        .class("btn")
        .attr("disabled", signal.map(move |v| v == &min))
        .on(ev::click, move |_| {
            signal.update(|c| {
                if *c > min + step {
                    *c -= step
                }
            })
        })
        .child("-");

    html::div()
        .class("control-edits")
        .class(("tooltip", tip.is_some()))
        .attr("data-text", tip)
        .child((
            html::span().class("control-label").child(name),
            slider,
            html::div().class("increment-buttons").child((
                minus,
                html::span()
                    .class("slider-value")
                    .child(signal.map(move |v| {
                        if *v == 0.0 {
                            zero.to_owned()
                        } else {
                            v.to_string()
                        }
                    })),
                plus,
            )),
        ))
}

fn render_full_phenotype_info(
    pgs_id: PgsId,
    init: &OutputInit,
    data: Load<Option<ArcRwSignal<OutputGetScores>>>,
    loading: Signal<bool>,
    normalise: bool,
) -> impl IntoView + use<> {
    fn label_value(label: &str, value: impl IntoView) -> impl IntoView {
        html::div().class("info-item").child((
            html::span().class("info-label").child(label.to_string()),
            html::span().class("info-value").child(value),
        ))
    }
    fn render_stats_section(t: &'static str, stats: Stats, loading: Signal<bool>) -> impl IntoView {
        html::div()
            .class("card")
            .class(("loading", loading))
            .child((
                html::h2().child(t),
                html::div().class("info-grid").child((
                    html::div().class("stats-row").child((
                        label_value("Mean: ", format!("{:.4}", stats.mean)),
                        label_value("Std Dev: ", format!("{:.4}", stats.std_dev)),
                    )),
                    html::div().class("stats-row").child((
                        label_value("Min: ", format!("{:.4}", stats.min)),
                        label_value("Max: ", format!("{:.4}", stats.max)),
                    )),
                )),
            ))
    }
    fn render_stats_section_placeholder(
        t: &'static str,
        msg: &'static str,
        loading: Signal<bool>,
    ) -> impl IntoView {
        html::div()
            .class("card")
            .class(("loading", loading))
            .child((
                html::h2().child(t),
                html::div().class("info-grid").child((
                    html::div()
                        .class("stats-row")
                        .child((label_value("Mean: ", msg), label_value("Std Dev: ", msg))),
                    html::div()
                        .class("stats-row")
                        .child((label_value("Min: ", msg), label_value("Max: ", msg))),
                )),
            ))
    }
    fn view_development_sample(sample: &ScoreDevelopmentSample) -> impl IntoView + use<> {
        html::div().class("population-item").child((
            html::h3().child(sample.stage_of_pgs_development.to_string()),
            html::div().class("population-details").child((
                sample.number_of_individuals.map(|n| {
                    html::div()
                        .class("population-stat")
                        .child(format!("Individuals: {n}",))
                }),
                sample.number_of_cases.map(|n| {
                    html::div()
                        .class("population-stat")
                        .child(format!("Cases: {n}"))
                }),
                sample.number_of_controls.map(|n| {
                    html::div()
                        .class("population-stat")
                        .child(format!("Controls: {n}"))
                }),
                (!sample.broad_ancestry_category.is_empty()).then(|| {
                    html::div().class("population-stat").child(format!(
                        "Ancestry: {}",
                        sample.broad_ancestry_category.clone()
                    ))
                }),
            )),
        ))
    }
    fn view_metric(metric: &PerformanceMetric) -> impl IntoView + use<> {
        html::div().class("performance-metric").child((
            html::h3().child(format!("Metric ID: {}", metric.id)),
            html::div().class("metric-details").child((
                (!metric.reported_trait.is_empty())
                    .then(|| label_value("Reported Trait: ", metric.reported_trait.clone())),
                (!metric.covariates_included_in_the_model.is_empty()).then(|| {
                    label_value(
                        "Covariates: ",
                        metric.covariates_included_in_the_model.clone(),
                    )
                }),
                (!metric.auroc.is_empty()).then(|| label_value("AUROC: ", metric.auroc.clone())),
                (!metric.concordance_statistic_c_index.is_empty()).then(|| {
                    label_value("C-index: ", metric.concordance_statistic_c_index.clone())
                }),
                (!metric.hazard_ratio.is_empty())
                    .then(|| label_value("Hazard Ratio (HR): ", metric.hazard_ratio.clone())),
                (!metric.odds_ratio.is_empty())
                    .then(|| label_value("Odds Ratio (OR): ", metric.odds_ratio.clone())),
                (!metric.beta.is_empty()).then(|| label_value("Beta: ", metric.beta.clone())),
                (!metric.other_metrics.is_empty())
                    .then(|| label_value("Other Metrics: ", metric.other_metrics.clone())),
                (!metric.other_relevant_information.is_empty()).then(|| {
                    label_value("Other Info: ", metric.other_relevant_information.clone())
                }),
                metric.publication_pmid.map(|pmid| {
                    html::div().class("info-item").child((
                        html::span().class("info-label").child("Publication: "),
                        html::span().class("info-value").child((
                            html::a()
                                .href(pmid.url())
                                .attr("target", "_blank")
                                .child(format!("PMID: {pmid}")),
                            (!metric.publication_doi.is_empty()).then(|| {
                                (
                                    " | ",
                                    html::a()
                                        .href(format!("https://doi.org/{}", metric.publication_doi))
                                        .attr("target", "_blank")
                                        .child(format!("DOI: {}", metric.publication_doi)),
                                )
                            }),
                        )),
                    ))
                }),
            )),
        ))
    }
    fn view_publication(pub_info: &Publication) -> impl IntoView + use<> {
        html::div().class("info-item").child((
            html::span()
                .class("info-label")
                .child(format!("{}: ", pub_info.first_author)),
            html::span().class("info-value").child((
                pub_info.title.clone(),
                " (",
                pub_info.pmid.map(|pmid| {
                    html::a()
                        .href(pmid.url())
                        .attr("target", "_blank")
                        .child(format!("PMID: {pmid}"))
                }),
                (pub_info.pmid.is_some() && !pub_info.doi.is_empty()).then_some(" | "),
                (!pub_info.doi.is_empty()).then(|| {
                    html::a()
                        .href(format!("https://doi.org/{}", pub_info.doi))
                        .attr("target", "_blank")
                        .child(format!("DOI: {}", pub_info.doi))
                }),
                ")",
            )),
        ))
    }

    let Metadata {
        cohorts: _,
        evaluation_sample_sets: _,
        performance_metrics,
        score_development_samples,
        scores,
        efo_traits: _,
        publications,
    } = init.metadata.0.clone();

    // Get the first score's information for display
    let score = scores.first().cloned().unwrap();

    let never_loading = Signal::derive(move || false);

    (
        // Statistics
        {
            let (title, stats) = match normalise {
                true => ("Normalised Statistics", init.norm_stats),
                false => ("Original Statistics", init.stats),
            };
            match stats {
                Some(stats) => render_stats_section(title, stats, never_loading).into_any(),
                None => render_stats_section_placeholder(title, "N/A", never_loading).into_any(),
            }
        },
        move || match &data {
            Load::Ready(Some(data)) => {
                render_stats_section("Edited Statistics", data.get().stats, loading).into_any()
            }
            Load::Ready(None) => {
                render_stats_section_placeholder("Edited Statistics", "N/A", loading).into_any()
            }
            Load::Loading => {
                render_stats_section_placeholder("Edited Statistics", "Loading…", loading)
                    .into_any()
            }
        },
        // Basic Information
        html::div().class("card").child((
            html::h2().child("Basic Information"),
            html::div().class("info-grid").child((
                html::div().class("info-item").child((
                    html::span().class("info-label").child("PGS ID: "),
                    html::a().href(pgs_id.url()).child(pgs_id.to_string()),
                )),
                label_value("PGS Name: ", score.name.clone()),
                label_value("Reported Trait: ", score.reported_trait.clone()),
                label_value("Mapped Traits: ", score.mapped_traits_efo_label.clone()),
                label_value("Number of Variants: ", score.number_of_variants),
                label_value(
                    "Variant Weight Type: ",
                    format!("{}", score.type_of_variant_weight),
                ),
                label_value("Development Method: ", score.pgs_development_method.clone()),
                label_value(
                    "Original Genome Build: ",
                    score.original_genome_build.clone(),
                ),
            )),
            if !score
                .pgs_development_details_and_relevant_parameters
                .is_empty()
            {
                Some(
                    html::div().class("info-grid").child(label_value(
                        "Development Details: ",
                        score
                            .pgs_development_details_and_relevant_parameters
                            .clone(),
                    )),
                )
            } else {
                None
            },
        )),
        // Population-specific Data
        html::div().class("card").child((
            html::h2().child("Score Development & Evaluation"),
            html::div().class("info-grid").child((
                label_value(
                    "Ancestry Distribution (GWAS): ",
                    score
                        .ancestry_distribution_source_of_variant_associations_gwas
                        .clone(),
                ),
                label_value(
                    "Ancestry Distribution (Development): ",
                    score
                        .ancestry_distribution_score_development_and_training
                        .clone(),
                ),
                label_value(
                    "Ancestry Distribution (Evaluation): ",
                    score.ancestry_distribution_pgs_evaluation.clone(),
                ),
            )),
            // Development samples information
            (!score_development_samples.is_empty()).then(|| {
                html::div().class("population-grid").child(
                    score_development_samples
                        .iter()
                        .filter(|s| s.score_id == score.id)
                        .map(view_development_sample)
                        .collect_view(),
                )
            }),
        )),
        // Performance Metrics
        (!performance_metrics.is_empty()).then(|| {
            html::div().class("card").child((
                html::h2().child("Performance Metrics"),
                html::div().class("info-grid").child(
                    performance_metrics
                        .iter()
                        .filter(|metric| metric.evaluated_score == score.id)
                        .map(view_metric)
                        .collect_view(),
                ),
            ))
        }),
        // Publications
        (!publications.is_empty()).then(|| {
            html::div().class("card").child((
                html::h2().child("Publications"),
                html::div()
                    .class("info-grid")
                    .child(publications.iter().map(view_publication).collect_view()),
            ))
        }),
    )
}

fn render_edit_analysis_panel(
    init: OutputInit,
    edits: Signal<Load<Option<ArcRwSignal<OutputGetEditAnalysis>>>>,
    edit_count: Signal<isize>,
) -> impl IntoView {
    let loading = edits.map(|o| !o.is_ready());
    let edits = edits.keep_if(|o| o.is_ready());

    let sorted_edits = Signal::derive(move || {
        let mut data = edits.get().ready()??.get();
        let OutputGetEditAnalysis { edits } = &mut data;
        edits.sort_unstable_by_key(|VariantInfo { association, .. }| {
            let simplified = association
                .0
                .clone()
                .association
                .simplified(parse_contig)
                .unwrap();
            let effect_weight = match simplified.effect {
                pgs_catalog::simplified::Effect::Additive { effect_weight, .. } => effect_weight,
                pgs_catalog::simplified::Effect::DosageSpecific {
                    dosage_2_weight, ..
                } => dosage_2_weight,
                pgs_catalog::simplified::Effect::Dominant { effect_weight } => effect_weight,
                pgs_catalog::simplified::Effect::Recessive { effect_weight } => effect_weight,
            };
            NotNan::new(-effect_weight.abs()).unwrap()
        });
        Some(data)
    });

    (
        html::div()
            .class("card card-chart")
            .class(("loading", loading))
            .child(move || match edits.get() {
                Load::Loading => html::div()
                    .class("loading-text")
                    .child("Loading…")
                    .into_any(),
                Load::Ready(None) => html::div()
                    .class("loading-text")
                    .child("Too many edits.")
                    .into_any(),
                Load::Ready(Some(_)) if edit_count.is(0).get() => html::div()
                    .class("loading-text")
                    .child("Apply at least one edit.")
                    .into_any(),
                Load::Ready(Some(analysis)) => render_edit_frequency_histogram(
                    edit_count,
                    init.sample_count,
                    analysis.map(|a| a.edits.iter().map(|v| v.count).collect()),
                )
                .into_any(),
            }),
        move || {
            // Rendering long lists can be kind of slow, this is a horrible hack that just renders it in chunks,
            // breaking up the work and keeping the interface responsive rather than implementing a proper virtual list.
            // But also, why is it slow at all? 500 elements is not that much.

            let len = sorted_edits.with(|s| s.as_ref().map(|s| s.edits.len()).unwrap_or(0));

            let loading = {
                let l = ArcRwSignal::new(loading.get_untracked());
                loading.for_each_after_first({
                    let l = l.clone();
                    move |v| l.set(*v)
                });
                l
            };

            let pieces = (len / 20) + 1;

            (0..pieces)
                .map(|i| {
                    let start = (i * 20).min(len);
                    let end = (start + 20).min(len);

                    let loading = loading.clone();

                    lazy(
                        move || {
                            sorted_edits.with(|sorted_edits| {
                                let sorted_edits = sorted_edits.as_ref()?;
                                Some(
                                    sorted_edits.edits[start..end]
                                        .iter()
                                        .map(|variant_info| {
                                            edit_table_line(variant_info, loading.clone().into())
                                        })
                                        .collect_view(),
                                )
                            })
                        },
                        i,
                    )
                })
                .collect_view()
        },
    )
}

/// Poor man's lazy loading.
fn lazy<T>(view: impl FnOnce() -> T + Send + Sync + 'static, i: usize) -> impl IntoView
where
    T: IntoView + 'static,
{
    let trigger = ArcTrigger::new();
    let task = Task::new_local({
        let trigger = trigger.clone();
        async move {
            utile::time::sleep(Duration::from_millis((i as u64) * 10)).await;
            trigger.notify();
        }
    });
    let mut view = Some(view);
    let mut first = true;
    move || {
        trigger.track();
        let _task = &task;

        if first {
            first = false;
            None
        } else {
            Some(view.take().unwrap()())
        }
    }
}

fn edit_table_line(variant_info: &VariantInfo, loading: ArcSignal<bool>) -> impl IntoView + use<> {
    let VariantInfo {
        count: _,
        association:
            AsJson(Association {
                reference_allele,
                association,
                ..
            }),
        ploidy_dosage,
    } = variant_info;

    let SimplifiedHarmonizedStudyAssociation {
        rs_id,
        chr,
        pos,
        effect_allele,
        other_allele: _,
        locus_name: _,
        kind: _,
        imputation_method: _,
        variant_description: _,
        inclusion_criteria: _,
        effect,
        allelefrequency_effect: _,
        allelefrequency_effect_european: _,
        allelefrequency_effect_asian: _,
        allelefrequency_effect_african: _,
        allelefrequency_effect_hispanic: _,
        variant_type: _,
        source: _,
    } = association.clone().simplified(parse_contig).unwrap();

    // Format dosage/ploidy information
    let dosage_info = if ploidy_dosage.is_empty() {
        "No dosage data".to_owned()
    } else {
        ploidy_dosage
            .iter()
            .take(5) // Show top 5 most common dosage/ploidy combinations
            .map(|(VariantSampleInfo { dosage, ploidy }, count)| {
                let s = if *count == 1 { "" } else { "s" };
                format!("{dosage} of {ploidy} ({count} sample{s})")
            })
            .collect::<Vec<_>>()
            .join(", ")
    };

    let effect = match effect {
        pgs_catalog::simplified::Effect::Additive {
            effect_weight,
            or: _,
            hr: _,
        } => format!("{effect_weight:.4}"),
        pgs_catalog::simplified::Effect::DosageSpecific {
            dosage_0_weight,
            dosage_1_weight,
            dosage_2_weight,
            effect_weight: _,
            or: _,
        } => format!("{dosage_0_weight:.4}/{dosage_1_weight:.4}/{dosage_2_weight:.4}"),
        pgs_catalog::simplified::Effect::Dominant { effect_weight } => {
            format!("{effect_weight:.4} (dominant)")
        }
        pgs_catalog::simplified::Effect::Recessive { effect_weight } => {
            format!("{effect_weight:.4} (recessive)")
        }
    };

    let gnomad_link = format!(
        "https://gnomad.broadinstitute.org/variant/{chr}-{pos}-{reference_allele}-{effect_allele}?dataset=gnomad_r4"
    );

    let ref_effect_match = reference_allele.len() == effect_allele.len()
        && reference_allele
            .iter()
            .zip(effect_allele.iter())
            .all(|(a, b)| a == b);

    html::div()
        .class("card edit-item")
        .class(("loading", loading))
        .child((
            html::div().class("edit-details").child((
                html::div().class("edit-primary").child((
                    html::span()
                        .class("edit-position")
                        .child(format!("{chr}:{pos}")),
                    if ref_effect_match {
                        html::span()
                            .class("edit-alleles")
                            .child(format!("{effect_allele}"))
                            .into_any()
                    } else {
                        (
                            html::span()
                                .class("edit-alleles")
                                .child(format!("{reference_allele} → {effect_allele}")),
                            html::a()
                                .class("gnomad-link")
                                .attr("href", gnomad_link)
                                .attr("target", "_blank")
                                .attr("rel", "noopener noreferrer")
                                .child("gnomAD"),
                        )
                            .into_any()
                    },
                    rs_id.map(|rs_id| {
                        html::span().child((
                            "rsID: ",
                            html::a()
                                .attr("href", format!("https://www.ncbi.nlm.nih.gov/snp/{rs_id}"))
                                .attr("target", "_blank")
                                .attr("rel", "noopener noreferrer")
                                .child(format!("{rs_id}")),
                        ))
                    }),
                )),
                html::div().class("edit-dosage").child((
                    html::span().class("dosage-label").child("Dosage×Ploidy: "),
                    html::span().class("dosage-info").child(dosage_info),
                )),
            )),
            html::div().child(html::span().class("count-badge").child(effect)),
        ))
}

fn toggle_button(
    signal: RwSignal<bool>,
    label: &'static str,
    tip: Option<&'static str>,
) -> impl IntoView {
    html::button()
        .class("toggle-btn")
        .class(("active", signal))
        .class(("tooltip", tip.is_some()))
        .attr("data-text", tip)
        .on(ev::click, move |_| signal.flip())
        .child(label)
}

impl PageState {
    fn new(pgs_id: PgsId) -> Self {
        let initial_state = ArcRwSignal::new(None);

        let bridge = PgsWorkerStruct::spawner()
            .callback({
                let initial_state = initial_state.clone();
                move |output| match output {
                    Output::Init(init) => initial_state.set(Some(init)),
                    Output::GetScores(_) | Output::GetEditAnalysis(_) => unreachable!(),
                }
            })
            .spawn("/pgs.js");

        bridge.send(Input::Set {
            pgs_id,
            origin: leptos::prelude::window().location().origin().unwrap(),
        });

        let bridge = Arc::new(SendOption::new_local(Some(bridge)));

        Self {
            pgs_id,

            bridge,

            init: initial_state,

            edit_count: RwSignal::new(0),
            normalise: RwSignal::new(false),

            show_stats: RwSignal::new(false),
            full_range: RwSignal::new(false),
        }
    }

    fn scores_signal(&self) -> Signal<Load<Option<ArcRwSignal<OutputGetScores>>>> {
        let bridge = self.bridge.clone();

        let edit_count = self.edit_count;
        let normalise = self.normalise;

        let data = RwSignal::new(Load::Loading);

        Signal::derive(move || Input::GetScores {
            edit_count: edit_count.get(),
            normalise: normalise.get(),
        })
        .rate_limit_leading(Duration::from_millis(10))
        .map_async(move |input| {
            let input = input.clone();
            let bridge = bridge.clone();
            async move {
                match query_bridge(&bridge, input).await {
                    Output::Init(_) | Output::GetEditAnalysis(_) => unreachable!(),
                    Output::GetScores(scores) => scores,
                }
            }
        })
        .for_each_immediate({
            let mut score_signal: Option<ArcRwSignal<_>> = None;

            move |scores| match scores {
                Load::Loading => data.set(Load::Loading),
                Load::Ready(None) => data.set(Load::Ready(None)),
                Load::Ready(Some(scores)) => {
                    let scores = scores.clone();
                    let s = match score_signal.clone() {
                        Some(s) => {
                            s.set(scores);
                            s
                        }
                        None => ArcRwSignal::new(scores),
                    };
                    score_signal = Some(s.clone());
                    if !data.read_untracked().is_ready() {
                        data.set(Load::Ready(Some(s.clone())));
                    }
                }
            }
        });

        data.into()
    }

    fn edits_signal(&self) -> Signal<Load<Option<ArcRwSignal<OutputGetEditAnalysis>>>> {
        let bridge = self.bridge.clone();

        let edit_count = self.edit_count;
        let normalise = self.normalise;

        let data = RwSignal::new(Load::Loading);

        Signal::derive(move || Input::GetEditAnalysis {
            edit_count: edit_count.get(),
            normalise: normalise.get(),
        })
        .rate_limit_leading(Duration::from_millis(10))
        .map_async(move |input| {
            let input = input.clone();
            let bridge = bridge.clone();
            async move {
                match query_bridge(&bridge, input).await {
                    Output::Init(_) | Output::GetScores(_) => unreachable!(),
                    Output::GetEditAnalysis(edits) => edits,
                }
            }
        })
        .for_each_immediate({
            let mut edits_signal: Option<ArcRwSignal<_>> = None;

            move |edits| match edits {
                Load::Loading => data.set(Load::Loading),
                Load::Ready(None) => data.set(Load::Ready(None)),
                Load::Ready(Some(edits)) => {
                    let edits = edits.clone();
                    let s = match edits_signal.clone() {
                        Some(s) => {
                            s.set(edits);
                            s
                        }
                        None => ArcRwSignal::new(edits),
                    };
                    edits_signal = Some(s.clone());
                    if !data.read_untracked().is_ready() {
                        data.set(Load::Ready(Some(s.clone())));
                    }
                }
            }
        });

        data.into()
    }
}

async fn query_bridge(bridge: &SendOption<WorkerBridge<PgsWorkerStruct>>, input: Input) -> Output {
    let (tx, rx) = futures::channel::oneshot::channel::<Output>();
    let bridge = {
        let callback = Mutex::new(Some(move |output| {
            let _ = tx.send(output);
        }));
        SendOption::new_local(Some((*bridge).as_ref().unwrap().fork(Some(
            move |output| {
                log::info!("[App] Received output");
                callback.lock().unwrap().take().unwrap()(output);
            },
        ))))
    };
    let start = web_time::Instant::now();
    log::info!("[App] Requesting {input:?}");
    bridge.as_ref().unwrap().send(input.clone());
    let output = rx.await.unwrap();
    log::info!(
        "[App] Received output for {input:?} in {:?}ms",
        start.elapsed().as_millis()
    );
    output
}
