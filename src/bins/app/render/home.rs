use leptos::{
    IntoView, ev,
    html::{self},
    prelude::*,
};

use leptos_ext::signal::ReadSignalExt;
use pan_ukbb::PhenotypeManifestEntry;

pub fn home() -> impl IntoView {
    let search = RwSignal::new(String::new());
    let min_variants = RwSignal::new(0usize);
    let min_date = RwSignal::new(String::new());

    html::div().class("page-content manifest").child((
        html::h1().child("Index"),
        html::section().inner_html(include_str!("intro.html")),
        html::div().class("search-bar").child((
            html::input()
                .attr("type", "text")
                .attr("placeholder", "Search scores…")
                .class("search-input")
                .on(ev::input, move |ev| {
                    search.set(event_target_value(&ev));
                }),
            html::div().class("filter-field").child((
                html::span().class("field-label").child("Min variants: "),
                html::input()
                    .attr("type", "number")
                    .attr("min", "0")
                    .attr("placeholder", "0")
                    .class("search-input variants-input")
                    .on(ev::input, move |ev| {
                        let v = event_target_value(&ev);
                        min_variants.set(v.parse().unwrap_or(0));
                    }),
            )),
            html::div().class("filter-field").child((
                html::span().class("field-label").child("Released after: "),
                html::input()
                    .attr("type", "date")
                    .class("search-input date-input")
                    .on(ev::input, move |ev| {
                        min_date.set(event_target_value(&ev));
                    }),
            )),
        )),
        move || {
            let metadata: ArcRwSignal<_> = crate::METADATA.with(|m| (**m).clone());
            metadata.with(|metadata| match metadata {
                Some(Ok(metadata)) => metadata
                    .scores
                    .iter()
                    .map(|score| metadata_entry(score, search, min_variants, min_date))
                    .collect_view()
                    .into_any(),
                Some(Err(e)) => {
                    log::error!("[App] Error loading metadata: {}", e);
                    html::p().child("Error loading the metadata.").into_any()
                }
                None => html::p().child("Loading…").into_any(),
            })
        },
        || {
            let manifest: ArcRwSignal<_> = crate::MANIFEST.with(|m| (**m).clone());
            manifest.with(|manifest| match manifest {
                Some(Ok(manifest)) => manifest
                    .iter()
                    .filter(|p| analysis::util::passes_qc(p))
                    .map(manifest_entry)
                    .collect_view()
                    .into_any(),
                Some(Err(e)) => {
                    log::error!("[App] Error loading manifest: {}", e);
                    html::p().child("Error loading the manifest.").into_any()
                }
                None => html::p().child("Loading…").into_any(),
            });
            {} // Hide
        },
    ))
}
fn manifest_entry(entry: &PhenotypeManifestEntry) -> impl IntoView + use<> {
    let PhenotypeManifestEntry {
        trait_type,
        phenocode,
        pheno_sex,
        description,
        description_more,
        category,
        n_cases_full_cohort_both_sexes,
        n_cases_hq_cohort_both_sexes,
        filename,
        ..
    } = entry;

    let identifier = filename.strip_suffix(".tsv.bgz").unwrap_or(filename);
    let url = format!("/pan_ukbb/{identifier}");

    let phenocode_url = format!("https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={phenocode}");

    html::div()
        .class("card")
        .child(html::h2().child(html::a().href(url).child(description.clone())))
        .child(
            html::div().class("entry-details").child((
                html::div().class("detail-row").child((
                    html::a()
                        .href(phenocode_url)
                        .child(format!("[{phenocode}]")),
                    format!("\u{00A0}({trait_type}) ({pheno_sex})"),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Category: "),
                    html::span()
                        .class("detail-value")
                        .child(category.clone().unwrap_or_default()),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Cases: "),
                    html::span().class("detail-value").child(format!(
                        "{n_cases_full_cohort_both_sexes} (HQ: {})",
                        n_cases_hq_cohort_both_sexes.unwrap_or(0)
                    )),
                )),
                description_more
                    .clone()
                    .map(|d| html::span().class("detail-value").child(d)),
            )),
        )
}

fn metadata_entry(
    score: &pgs_catalog::metadata::Score,
    search: RwSignal<String>,
    min_variants: RwSignal<usize>,
    min_date: RwSignal<String>,
) -> impl IntoView + use<> {
    let search_text = score_search_text(score);
    let filter_variants = score.number_of_variants;
    let filter_date = score.release_date.clone();

    let matches_search = Signal::derive(move || {
        let query = search.get().to_lowercase();
        query.is_empty() || search_text.contains(&query)
    });
    let meets_variants = Signal::derive(move || filter_variants >= min_variants.get());
    let meets_date = Signal::derive(move || {
        let d = min_date.get();
        d.is_empty() || filter_date.as_str() >= d.as_str()
    });
    let visible =
        Signal::derive(move || matches_search.get() && meets_variants.get() && meets_date.get());

    let pgs_catalog::metadata::Score {
        id,
        name,
        reported_trait,
        mapped_traits_efo_label,
        mapped_traits_efo_id: _,
        pgs_development_method,
        pgs_development_details_and_relevant_parameters: _,
        original_genome_build: _,
        number_of_variants,
        number_of_interaction_terms: _,
        type_of_variant_weight: _,
        pgs_publication_id: _,
        publication_pmid,
        publication_doi,
        score_and_results_match_the_original_publication: _,
        ancestry_distribution_source_of_variant_associations_gwas: _,
        ancestry_distribution_score_development_and_training: _,
        ancestry_distribution_pgs_evaluation: _,
        ftp_link: _,
        release_date,
        license_and_terms_of_use: _,
    } = score;

    let url = format!("/pgs_catalog/{id}");

    html::div()
        .class("card")
        .attr("hidden", visible.not())
        .child(html::h2().child(html::a().href(url).child(name.clone())))
        .child(
            html::div().class("entry-details").child((
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("PGS ID: "),
                    html::a().href(id.url()).child(id.to_string()),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Trait: "),
                    html::span()
                        .class("detail-value")
                        .child(reported_trait.clone()),
                )),
                (!mapped_traits_efo_label.is_empty()).then(|| {
                    html::div().class("detail-row").child((
                        html::span().class("detail-label").child("EFO Traits: "),
                        html::span()
                            .class("detail-value")
                            .child(mapped_traits_efo_label.clone()),
                    ))
                }),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Method: "),
                    html::span()
                        .class("detail-value")
                        .child(pgs_development_method.clone()),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Variants: "),
                    html::span()
                        .class("detail-value")
                        .child(*number_of_variants),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Released: "),
                    html::span()
                        .class("detail-value")
                        .child(release_date.clone()),
                )),
                html::div().class("detail-row").child((
                    html::span().class("detail-label").child("Publication: "),
                    html::span().class("detail-value").child({
                        let mut links = Vec::new();

                        if let Some(pmid) = publication_pmid {
                            links.push(
                                html::a()
                                    .href(pmid.url())
                                    .target("_blank")
                                    .child(format!("[PMID: {pmid}]"))
                                    .into_any(),
                            );
                        }

                        if !publication_doi.is_empty() {
                            links.push(
                                html::a()
                                    .href(format!("https://doi.org/{}", publication_doi))
                                    .target("_blank")
                                    .child("[DOI]")
                                    .into_any(),
                            );
                        }

                        if links.len() == 2 {
                            // Insert a space between the two links
                            links.insert(1, html::span().child(" ").into_any());
                        }

                        links
                    }),
                )),
            )),
        )
}

fn score_search_text(score: &pgs_catalog::metadata::Score) -> String {
    let pgs_catalog::metadata::Score {
        id,
        name,
        reported_trait,
        mapped_traits_efo_label,
        mapped_traits_efo_id,
        pgs_development_method,
        pgs_development_details_and_relevant_parameters,
        original_genome_build,
        number_of_variants,
        number_of_interaction_terms,
        type_of_variant_weight,
        pgs_publication_id,
        publication_pmid,
        publication_doi,
        score_and_results_match_the_original_publication,
        ancestry_distribution_source_of_variant_associations_gwas,
        ancestry_distribution_score_development_and_training,
        ancestry_distribution_pgs_evaluation,
        ftp_link,
        release_date,
        license_and_terms_of_use,
    } = score;

    let fields: [String; _] = [
        id.to_string(),
        name.to_string(),
        reported_trait.to_string(),
        mapped_traits_efo_label.join(" "),
        mapped_traits_efo_id.join(" "),
        pgs_development_method.to_string(),
        pgs_development_details_and_relevant_parameters.to_string(),
        original_genome_build.to_string(),
        number_of_variants.to_string(),
        number_of_interaction_terms.to_string(),
        type_of_variant_weight.to_string(),
        pgs_publication_id.to_string(),
        publication_pmid.map(|p| p.to_string()).unwrap_or_default(),
        publication_doi.to_string(),
        score_and_results_match_the_original_publication.to_string(),
        ancestry_distribution_source_of_variant_associations_gwas.join(" "),
        ancestry_distribution_score_development_and_training.join(" "),
        ancestry_distribution_pgs_evaluation.join(" "),
        ftp_link.to_string(),
        release_date.to_string(),
        license_and_terms_of_use.to_string(),
    ];
    fields.join("\n").to_lowercase()
}
