use biocore::{dna::IupacDnaSequence, location::ContigPosition};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use std::{
    cell::LazyCell,
    cmp::Ordering,
    collections::BTreeMap,
    ops::Range,
    sync::{Arc, atomic::AtomicU64},
};

use genomes1000::{
    contig::GRCh38Contig,
    pedigree::{Pedigree, Sex},
    simplified::SimplifiedRecord,
};
use pgs_catalog::{
    HarmonizedStudyAssociation,
    simplified::{SimplificationError, SimplifiedHarmonizedStudyAssociation},
};
use utile::collections::counting_set::CountingBTreeSet;

use crate::util::{MArc, Stats, load_pedigrees, mean, std_dev};

const EDIT_COUNT: usize = 100;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct PgsCatalogScores {
    pub scores: BTreeMap<String, PgsCatalogScore>,

    #[serde(with = "utile::serde_ext::as_vec")]
    pub pre_processing_errors: CountingBTreeSet<SimplificationError>,

    pub found: usize,
    pub found_reference: usize,
    pub not_found: usize,

    #[serde(with = "utile::serde_ext::as_vec")]
    pub scoring_errors: CountingBTreeSet<ScoringError>,

    /// Filled-in at the end before serialization.
    pub associations: BTreeMap<PgsCatalogAssociationKey, Association>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[derive(Serialize, Deserialize)]
pub struct PgsCatalogScore {
    pub score: NotNan<f64>,

    pub best: Vec<PgsCatalogVariant>,
    pub worst: Vec<PgsCatalogVariant>,

    pub sex: Sex,
    pub population: String,
    pub superpopulation: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct PgsCatalogVariant {
    pub dosage: u8,
    pub ploidy: u8,

    pub max_edit: NotNan<f64>,
    pub min_edit: NotNan<f64>,

    pub key: PgsCatalogAssociationKey,

    #[serde(skip)]
    pub _handle: Option<Arc<Association>>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub struct Association {
    /// An internal ID for de-duplication. Unique only within a single [PgsCatalogScores].
    pub key: PgsCatalogAssociationKey,
    pub reference_allele: IupacDnaSequence,
    pub association: HarmonizedStudyAssociation,
}
impl Association {
    pub fn new(
        reference_allele: IupacDnaSequence,
        association: HarmonizedStudyAssociation,
    ) -> Self {
        Self {
            key: PgsCatalogAssociationKey::next(),
            reference_allele,
            association,
        }
    }
}

impl PgsCatalogScores {
    pub async fn new() -> Self {
        let mut scores: BTreeMap<String, PgsCatalogScore> = BTreeMap::new();
        load_pedigrees()
            .await
            .into_iter()
            .map(|p| (p.id.clone(), PgsCatalogScore::new(p)))
            .collect_into(&mut scores);

        Self {
            scores,

            pre_processing_errors: CountingBTreeSet::new(),
            scoring_errors: CountingBTreeSet::new(),

            not_found: 0,
            found_reference: 0,
            found: 0,

            associations: BTreeMap::new(),
        }
    }

    pub fn push_missing(
        &mut self,
        association: Association,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
    ) {
        self.not_found += 1;
        let association = &mut MArc::new(association);
        for score in self.scores.values_mut() {
            score.push_missing(association, simplified);
        }
    }
    pub fn push_variant(
        &mut self,
        association: Association,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        self.found += 1;
        self._push_variant(association, simplified, sample_names, record);
    }
    pub fn push_fallback(
        &mut self,
        association: Association,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        self.found_reference += 1;
        self._push_variant(association, simplified, sample_names, record);
    }
    fn _push_variant(
        &mut self,
        association: Association,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
        sample_names: &[String],
        record: &SimplifiedRecord,
    ) {
        let association = &mut MArc::new(association);
        for (name, sample) in sample_names.iter().zip(&record.samples) {
            let score = self.scores.get_mut(name).unwrap();
            if let Some(ploidy) = sample.ploidy()
                && ploidy != record.contig.ploidy(score.sex)
            {
                self.scoring_errors
                    .increment(ScoringError::InvalidGenotype {
                        at: record.at(),
                        sample: name.clone(),
                    });
                continue;
            }
            score.push_variant(
                association,
                simplified,
                sample.dosage(1),
                record.contig.ploidy(score.sex),
            );
        }
    }
    pub fn push_error(&mut self, error: SimplificationError) {
        self.pre_processing_errors.increment(error);
    }

    pub fn finalize(&mut self) {
        for score in self.scores.values_mut() {
            for variant in &mut score.best {
                self.associations.insert(variant.key, variant.association());
            }
            for variant in &mut score.worst {
                self.associations.insert(variant.key, variant.association());
            }
        }
    }

    pub fn len(&self) -> usize {
        self.scores.len()
    }
    pub fn is_empty(&self) -> bool {
        self.scores.is_empty()
    }

    pub fn scores(&self) -> impl Iterator<Item = &PgsCatalogScore> {
        self.scores.values()
    }

    pub fn full_score_range(&self) -> Option<Range<f64>> {
        let (min, max) = self.scores.values().map(|s| s.full_score_range()).fold(
            (None, None),
            |(min, max): (Option<f64>, Option<f64>), range| {
                let min = min.unwrap_or(range.start).min(range.start);
                let max = max.unwrap_or(range.end).max(range.end);
                (Some(min), Some(max))
            },
        );
        Some(min?..max?)
    }

    pub fn edited_scores(&self, n: isize) -> impl Iterator<Item = Option<f64>> + use<'_> {
        self.scores.values().map(move |s| s.edited_score(n))
    }

    pub fn max_edit_count(&self) -> isize {
        self.scores()
            .map(|s| s.best.len() as isize)
            .max()
            .unwrap_or(0)
    }
    pub fn min_edit_count(&self) -> isize {
        self.scores()
            .map(|s| s.worst.len() as isize)
            .min()
            .unwrap_or(0)
    }

    pub fn stats(&self) -> Option<Stats> {
        Some(Stats {
            mean: *self.mean()?,
            std_dev: *self.std_dev().unwrap(),
            min: *self.min().unwrap(),
            max: *self.max().unwrap(),
        })
    }
    pub fn mean(&self) -> Option<NotNan<f64>> {
        Some(NotNan::new(mean(self.scores().map(|s| *s.score))?).unwrap())
    }
    pub fn std_dev(&self) -> Option<NotNan<f64>> {
        Some(NotNan::new(std_dev(*self.mean()?, self.scores().map(|s| *s.score))?).unwrap())
    }
    pub fn min(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score).min()
    }
    pub fn max(&self) -> Option<NotNan<f64>> {
        self.scores().map(|s| s.score).max()
    }

    pub fn normalised(&self) -> Self {
        if self.is_empty() {
            return self.clone();
        }

        let mean = self.mean().unwrap();
        let std_dev = self.std_dev().unwrap();

        Self {
            scores: self
                .scores
                .clone()
                .into_iter()
                .map(|(k, v)| (k, v.normalised(mean, std_dev)))
                .collect(),

            pre_processing_errors: self.pre_processing_errors.clone(),
            found: self.found,
            found_reference: self.found_reference,
            not_found: self.not_found,
            scoring_errors: self.scoring_errors.clone(),

            associations: self
                .associations
                .clone()
                .into_iter()
                .map(|(k, s)| (k, s.normalised(std_dev)))
                .collect(),
        }
    }
}

impl PgsCatalogScore {
    fn new(pedigree: Pedigree) -> Self {
        Self {
            score: NotNan::new(0.0).unwrap(),

            best: Vec::with_capacity(EDIT_COUNT + 1),
            worst: Vec::with_capacity(EDIT_COUNT + 1),

            sex: pedigree.sex,
            population: pedigree.population,
            superpopulation: pedigree.superpopulation,
        }
    }
    fn push_missing(
        &mut self,
        a: &mut MArc<Association>,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
    ) {
        let ploidy = match (self.sex, simplified.chr) {
            (Sex::Male, GRCh38Contig::X) => 1,
            (Sex::Female, GRCh38Contig::X) => 2,
            (Sex::Male, GRCh38Contig::Y) => 1,
            (Sex::Female, GRCh38Contig::Y) => 0,
            (_, GRCh38Contig::MT) => 1,
            (_, _) => 2,
        };
        self.push_variant(a, simplified, 0, ploidy);
    }
    fn push_variant(
        &mut self,
        a: &mut MArc<Association>,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
        dosage: u8,
        ploidy: u8,
    ) {
        {
            let a = &**a;
            assert!(
                dosage <= ploidy,
                "dosage: {dosage}, ploidy: {ploidy}, s: {a:?}"
            );
            assert!(ploidy <= 2, "dosage: {dosage}, ploidy: {ploidy}, s: {a:?}");
        }
        self.score += simplified.effect.score(dosage, ploidy);

        let max = simplified.max_edit(dosage, ploidy);
        let min = simplified.min_edit(dosage, ploidy);

        let v = LazyCell::new(move || PgsCatalogVariant::new(a.arc(), simplified, dosage, ploidy));

        if self.best_bound() < max {
            let i = self.best.partition_point(|v| v.max_edit >= max);
            if i < EDIT_COUNT {
                self.best.insert(i, v.clone());
                self.best.truncate(EDIT_COUNT);
            }
        }
        if min < self.worst_bound() {
            let i = self.worst.partition_point(|v| v.min_edit <= min);
            if i < EDIT_COUNT {
                self.worst.insert(i, v.clone());
                self.worst.truncate(EDIT_COUNT);
            }
        }
    }
    fn best_bound(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.best.len() < EDIT_COUNT {
            zero
        } else {
            self.best.last().map(|v| v.max_edit).unwrap_or(zero)
        }
    }
    fn worst_bound(&self) -> NotNan<f64> {
        let zero = NotNan::new(0.).unwrap();
        if self.worst.len() < EDIT_COUNT {
            zero
        } else {
            self.worst.last().map(|v| v.min_edit).unwrap_or(zero)
        }
    }
}
impl PgsCatalogScore {
    pub fn full_score_range(&self) -> Range<f64> {
        let min = self.edited_score(-(self.worst.len() as isize)).unwrap();
        let max = self.edited_score(self.best.len() as isize).unwrap();
        min..max
    }
    pub fn edited_score(&self, n: isize) -> Option<f64> {
        Some(
            self.score
                + match Ord::cmp(&n, &0) {
                    Ordering::Greater => self
                        .best
                        .get(..n as usize)?
                        .iter()
                        .map(|v| *v.max_edit)
                        .sum::<f64>(),
                    Ordering::Less => self
                        .worst
                        .get(..-n as usize)?
                        .iter()
                        .map(|v| *v.min_edit)
                        .sum::<f64>(),
                    Ordering::Equal => 0.,
                },
        )
    }

    pub fn normalised(&self, mean: NotNan<f64>, std_dev: NotNan<f64>) -> Self {
        Self {
            score: (self.score - mean) / std_dev,

            best: self.best.iter().map(|v| v.normalised(std_dev)).collect(),
            worst: self.worst.iter().map(|v| v.normalised(std_dev)).collect(),

            ..self.clone()
        }
    }
}

impl PgsCatalogVariant {
    pub fn new(
        a: Arc<Association>,
        simplified: &SimplifiedHarmonizedStudyAssociation<GRCh38Contig>,
        dosage: u8,
        ploidy: u8,
    ) -> Self {
        Self {
            dosage,
            ploidy,

            max_edit: simplified.max_edit(dosage, ploidy),
            min_edit: simplified.min_edit(dosage, ploidy),

            key: a.key,

            _handle: Some(a),
        }
    }
    pub fn association(&self) -> Association {
        let s: &Association = self._handle.as_ref().unwrap();
        s.clone()
    }
    pub fn normalised(&self, std_dev: NotNan<f64>) -> Self {
        assert_eq!(None, self._handle);
        Self {
            max_edit: self.max_edit / std_dev,
            min_edit: self.min_edit / std_dev,

            ..self.clone()
        }
    }
}

impl Association {
    pub fn normalised(self, std_dev: NotNan<f64>) -> Self {
        Self {
            key: self.key,
            reference_allele: self.reference_allele,
            association: self.association.normalised(std_dev),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[derive(Serialize, Deserialize)]
pub enum ScoringError {
    InvalidGenotype {
        at: ContigPosition<GRCh38Contig>,
        sample: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PgsCatalogAssociationKey(u64);
impl PgsCatalogAssociationKey {
    pub fn next() -> Self {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        Self(COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed))
    }
}
impl std::fmt::Display for PgsCatalogAssociationKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
// Round-trip as a string so it can be used as a JSON key.
impl Serialize for PgsCatalogAssociationKey {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.0.to_string())
    }
}
impl<'de> Deserialize<'de> for PgsCatalogAssociationKey {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        Ok(Self(s.parse().unwrap()))
    }
}

pub fn parse_contig(c: String) -> GRCh38Contig {
    if let Some(c) = genomes1000::contig::GRCh38Contig::new(&c) {
        return c;
    }
    if let Some(c) = genomes1000::contig::GRCh38Contig::new(&format!("chr{c}")) {
        return c;
    }
    match &*c {
        "MT" => genomes1000::contig::GRCh38Contig::MT,

        _ => {
            panic!("[pgs_scores] No contig found for {c:?}")
        }
    }
}
