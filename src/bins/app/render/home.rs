use leptos::{
    html::{self},
    prelude::*,
    IntoView,
};

use pan_ukbb::PhenotypeManifestEntry;

pub fn home() -> impl IntoView {
    html::div().class("page-content manifest").child((
        html::h1().child("Index"),
        html::p().inner_html(include_str!("intro.html")),
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
            })
        },
        || {
            let metadata: ArcRwSignal<_> = crate::METADATA.with(|m| (**m).clone());
            metadata.with(|metadata| match metadata {
                Some(Ok(metadata)) => metadata
                    .scores
                    .iter()
                    .map(metadata_entry)
                    .collect_view()
                    .into_any(),
                Some(Err(e)) => {
                    log::error!("[App] Error loading metadata: {}", e);
                    html::p().child("Error loading the metadata.").into_any()
                }
                None => html::p().child("Loading…").into_any(),
            })
        },
    ))
}
fn manifest_entry(entry: &PhenotypeManifestEntry) -> impl IntoView {
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

fn metadata_entry(score: &pgs_catalog::metadata::Score) -> impl IntoView {
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
