use biocore::{
    dna::{DnaSequence, IupacDnaBase, IupacDnaSequence},
    location::ContigRange,
};
use futures::{StreamExt, stream};
use std::{path::PathBuf, sync::LazyLock};

use genomes1000::{
    DiploidGenotype, Genomes1000Fs, Genotype, GenotypePhasing, HaploidGenotype,
    contig::GRCh38Contig, simplified::SimplifiedRecord, source::Genomes1000Resource,
};
use pgs_catalog::{HarmonizedStudyAssociation, PgsId};
use resource::{
    RawResourceExt,
    fs::{FsCache, FsCacheEntry},
};

use analysis::{
    pgs_scores::{Association, PgsCatalogScores, parse_contig},
    util::InspectEvery,
};

type FastaReader = biocore::fasta::IndexedFastaReader<std::io::BufReader<std::fs::File>>;

static OUTPUT_FOLDER: LazyLock<FsCache> = LazyLock::new(|| {
    let path = PathBuf::from("./data/pgs_catalog");
    let path = path.canonicalize().unwrap();
    log::info!("[PGS Catalog] Output folder: {path:?}");
    FsCache::new(path)
});

fn all_metadata_output_path() -> FsCacheEntry {
    OUTPUT_FOLDER.entry("metadata/all.json.br")
}
fn metadata_output_path(id: PgsId) -> FsCacheEntry {
    OUTPUT_FOLDER.entry(format!("metadata/{id}.json.br"))
}
fn scores_output_path(id: PgsId) -> FsCacheEntry {
    OUTPUT_FOLDER.entry(format!("scores/{id}.json.br"))
}

#[tokio::main]
async fn main() {
    const CONCURRENCY: usize = 5;

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .filter_module("reqwest", log::LevelFilter::Info)
        .filter_module("hyper_util", log::LevelFilter::Info)
        .filter_module("resource", log::LevelFilter::Warn)
        .init();

    if !all_metadata_output_path().try_exists().unwrap() {
        log::info!("[PGS Catalog][all] Writing metadata");
        let metadata = pgs_catalog::metadata::Metadata::load_all().await.unwrap();
        all_metadata_output_path()
            .write_file_with(|file| {
                let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
                serde_json::to_writer(file, &metadata).unwrap();
                Ok(())
            })
            .unwrap();
        log::info!("[PGS Catalog][all] Finished writing metadata");
    }

    let pgs_ids = pgs_catalog::metadata::Metadata::load_all()
        .await
        .unwrap()
        .scores
        .into_iter()
        .map(|s| s.id)
        .filter(|id| {
            let exists = scores_output_path(*id).try_exists().unwrap();
            if exists {
                log::info!("[PGS Catalog][{id}] Found cached value");
            }
            !exists
        })
        .inspect(|id| log::info!("[PGS Catalog][{id}] Processing"))
        .map(run_pgs_scores)
        .map(|fut| async move { tokio::spawn(fut).await.unwrap() });

    stream::iter(pgs_ids)
        .buffer_unordered(CONCURRENCY)
        .for_each(|()| async {})
        .await;
}

async fn run_pgs_scores(id: PgsId) {
    std::fs::create_dir_all(scores_output_path(id).as_ref().parent().unwrap()).unwrap();

    if !metadata_output_path(id).try_exists().unwrap() {
        log::info!("[PGS Catalog][{id}] Writing metadata");
        let metadata = pgs_catalog::metadata::Metadata::load(id).await.unwrap();
        metadata_output_path(id)
            .write_file_with(|file| {
                let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
                serde_json::to_writer(file, &metadata).unwrap();
                Ok(())
            })
            .unwrap();
        log::info!("[PGS Catalog][{id}] Finished writing metadata");
    }

    if scores_output_path(id).try_exists().unwrap() {
        log::info!("[PGS Catalog][{id}] Found cached value");
        return;
    }

    log::info!("[PGS Catalog][{id}] Scoring");

    let associations = pgs_catalog::HarmonizedStudy::load_associations(
        pgs_catalog::PgsCatalogResource::HarmonizedStudy {
            id,
            build: pgs_catalog::GenomeBuild::GRCh38,
        }
        .log_progress()
        .with_global_fs_cache()
        .ensure_cached_async()
        .await
        .unwrap()
        .decompressed()
        .buffered(),
    )
    .unwrap()
    .map(|a| a.unwrap());

    let mut scores = score(id, associations).await;

    scores.finalize();

    log::info!("[PGS Catalog][{id}] Writing");

    scores_output_path(id)
        .write_file_with(|file| {
            let file = brotli::CompressorWriter::new(file, 4096, 9, 20);
            serde_json::to_writer(file, &scores).unwrap();
            Ok(())
        })
        .unwrap();

    log::info!("[PGS Catalog][{id}] Finished scoring");
}

async fn score(
    id: PgsId,
    associations: impl Iterator<Item = HarmonizedStudyAssociation>,
) -> PgsCatalogScores {
    let mut scores = PgsCatalogScores::new().await;

    let mut genomes1000 = Genomes1000Fs::new().await.unwrap();
    let mut grch38: FastaReader = load_grch38_reference_genome().await.unwrap();

    let sample_names = genomes1000.sample_names().to_vec();

    let associations = associations.inspect_every(100_000, |i, _| {
        if i != 0 {
            log::info!("[pgs_scores][{id}] Processing association {i}");
        }
    });

    'association: for association in associations {
        let simplified = association.clone().simplified(parse_contig);
        let simplified = match simplified {
            Ok(simplified) => simplified,
            Err(e) => {
                scores.push_error(e);
                continue;
            }
        };

        let at = simplified.at();
        let at_genomes1000 = at.map_contig(to_genomes1000_contig);

        let association = Association::new(
            grch38
                .query::<IupacDnaBase, _>(&ContigRange {
                    contig: simplified.chr,
                    at: simplified.pos - 1
                        ..(simplified.pos - 1
                            + u64::try_from(simplified.effect_allele.len()).unwrap()),
                })
                .unwrap(),
            association,
        );

        let mut not_ref_genotypes = vec![Genotype::Missing; sample_names.len()];
        for record in genomes1000
            .query_simplified(&at_genomes1000.into())
            .unwrap()
        {
            if record.at() != at_genomes1000 {
                // In theory this could still overlap, but in practice all the
                // tooling assumes this.
                continue;
            }

            if record.alternate_allele == simplified.effect_allele {
                // If we have found the effect allele, we are done.

                // log::info!("[pgs_scores][{id}] Found match for {at:?}");
                scores.push_variant(association, &simplified, &sample_names, &record);
                continue 'association;
            }

            // Otherwise aggregate all the non-reference genotypes so that we can
            // infer the true reference genotype.
            for (i, sample) in record.samples.iter().enumerate() {
                match sample {
                    Genotype::Missing => {}
                    Genotype::Haploid(HaploidGenotype { value: 0 }) => {}
                    Genotype::Diploid(DiploidGenotype {
                        left: 0,
                        right: 0,
                        phasing: _,
                    }) => {}

                    Genotype::Haploid(HaploidGenotype { value }) => {
                        match &mut not_ref_genotypes[i] {
                            Genotype::Missing => not_ref_genotypes[i] = *sample,
                            Genotype::Haploid(HaploidGenotype {
                                value: non_ref_value,
                            }) => {
                                *non_ref_value += *value;
                            }
                            Genotype::Diploid(_) => {
                                unreachable!("{sample:?}")
                            }
                        }
                    }
                    Genotype::Diploid(DiploidGenotype {
                        left,
                        right,
                        phasing: _,
                    }) => match &mut not_ref_genotypes[i] {
                        Genotype::Missing => not_ref_genotypes[i] = *sample,
                        Genotype::Diploid(DiploidGenotype {
                            left: non_ref_left,
                            right: non_ref_right,
                            phasing: _,
                        }) => {
                            *non_ref_left += *left;
                            *non_ref_right += *right;
                        }
                        Genotype::Haploid(_) => {
                            unreachable!("{sample:?}")
                        }
                    },
                }
            }
        }

        if eq_seq(&association.reference_allele, &simplified.effect_allele) {
            // log::info!("[pgs_scores][{id}] Found reference fallback match for {at:?}");

            // The effect allele is the reference allele, so we need to retrieve the
            // reference genotypes by inverting the aggregation of all the non-reference genotypes.

            let samples = not_ref_genotypes
                .into_iter()
                .enumerate()
                .map(|(i, g)| match g {
                    Genotype::Missing => {
                        // No non-reference alleles were seen for this sample
                        match simplified
                            .chr
                            .ploidy(genomes1000.pedigree(&sample_names[i]).unwrap().sex)
                        {
                            0 => Genotype::Missing,
                            1 => Genotype::Haploid(HaploidGenotype { value: 1 }),
                            2 => Genotype::Diploid(DiploidGenotype {
                                left: 1,
                                phasing: GenotypePhasing::Phased,
                                right: 1,
                            }),
                            _ => unreachable!(),
                        }
                    }

                    // We use saturating sub, because sometimes both A -> T and A -> CCT might
                    // be in the data (note that CTT ends in T).
                    Genotype::Haploid(haploid_genotype) => Genotype::Haploid(HaploidGenotype {
                        value: 1u8.saturating_sub(haploid_genotype.value),
                    }),
                    Genotype::Diploid(diploid_genotype) => Genotype::Diploid(DiploidGenotype {
                        left: 1u8.saturating_sub(diploid_genotype.left),
                        right: 1u8.saturating_sub(diploid_genotype.right),
                        phasing: diploid_genotype.phasing,
                    }),
                })
                .collect();

            let record = SimplifiedRecord {
                contig: to_genomes1000_contig(simplified.chr),
                position: simplified.pos,
                reference_allele: simplified.effect_allele.clone(),
                alternate_allele: simplified.effect_allele.clone(),
                quality: None,
                filter: "".to_string(),
                samples,
            };
            scores.push_fallback(association, &simplified, &sample_names, &record);
        } else {
            // log::info!("[pgs_scores][{id}] No match found for {at:?}");
            scores.push_missing(association, &simplified);
        }
    }

    log::info!(
        "[pgs_scores][{id}] {} found; {} found_reference; {} not_found; {}/{} errors",
        scores.found,
        scores.found_reference,
        scores.not_found,
        scores.pre_processing_errors.len(),
        scores.scoring_errors.len(),
    );

    scores
}

fn eq_seq(a: &IupacDnaSequence, b: &DnaSequence) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(a, b)| a == b)
}

async fn load_grch38_reference_genome()
-> std::io::Result<biocore::fasta::IndexedFastaReader<std::io::BufReader<std::fs::File>>> {
    let resource = Genomes1000Resource::grch38_reference_genome()
        .log_progress()
        .with_global_fs_cache()
        .ensure_cached_async()
        .await?;
    let index_resource = Genomes1000Resource::grch38_reference_genome_index()
        .log_progress()
        .with_global_fs_cache()
        .ensure_cached_async()
        .await?;

    genomes1000::load_grch38_reference_genome(resource.buffered(), index_resource).await
}

fn to_genomes1000_contig(c: GRCh38Contig) -> genomes1000::contig::GRCh38Contig {
    let c = c.as_ref();
    if let Some(c) = genomes1000::contig::GRCh38Contig::new(c) {
        return c;
    }
    if let Some(c) = genomes1000::contig::GRCh38Contig::new(&format!("chr{c}")) {
        return c;
    }
    match c {
        "MT" => genomes1000::contig::GRCh38Contig::MT,

        _ => {
            panic!("[pgs_scores] No contig found for {c:?}")
        }
    }
}
