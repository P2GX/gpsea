import pytest

import hpotk

from gpsea.model import *
from gpsea.model.genome import *
from gpsea.preprocessing import ProteinMetadataService

from gpsea.analysis.predicate.genotype import ProteinPredicates


@pytest.fixture(scope="package")
def protein_metadata_service() -> ProteinMetadataService:
    response = ProteinMetadata(
        protein_id="pt:xyz",
        label="xyz_label",
        protein_features=(
            ProteinFeature.create(
                FeatureInfo(name="MOCK_REPEAT", region=Region(55, 80)),
                FeatureType.REPEAT,
            ),
            ProteinFeature.create(
                FeatureInfo(name="MOCK_DOMAIN", region=Region(30, 50)),
                FeatureType.DOMAIN,
            ),
        ),
    )
    return MockProteinMetadataService(response)


@pytest.fixture(scope="package")
def protein_predicates(
    protein_metadata_service: ProteinMetadataService,
) -> ProteinPredicates:
    return ProteinPredicates(protein_metadata_service)


class MockProteinMetadataService(ProteinMetadataService):

    def __init__(self, response: ProteinMetadata):
        self._response = response

    def annotate(self, protein_id: str) -> ProteinMetadata:
        return self._response


@pytest.fixture(scope="package")
def sample_labels() -> SampleLabels:
    return SampleLabels("jim")


@pytest.fixture(scope="package")
def missense_variant(
    genome_build: GenomeBuild,
    sample_labels: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=100,
                    end=101,
                    strand=Strand.POSITIVE,
                ),
                ref="C",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(40, 41),
            ),
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:abc",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(VariantEffect.INTRON_VARIANT,),
                affected_exons=None,
                protein_id=None,
                hgvsp=None,
                protein_effect_coordinates=None,
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HETEROZYGOUS),
    )


@pytest.fixture(scope="package")
def frameshift_variant(
    genome_build: GenomeBuild,
    sample_labels: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=110,
                    end=113,
                    strand=Strand.POSITIVE,
                ),
                ref="CCC",
                alt="C",
                change_length=-2,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(VariantEffect.FRAMESHIFT_VARIANT,),
                affected_exons=(5,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(43, 44),
            ),
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:abc",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(VariantEffect.INTRON_VARIANT,),
                affected_exons=None,
                protein_id=None,
                hgvsp=None,
                protein_effect_coordinates=None,
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HETEROZYGOUS),
    )


@pytest.fixture(scope="package")
def structural_variant(
    sample_labels: SampleLabels,
) -> Variant:
    return Variant(
        variant_info=VariantInfo(
            sv_info=ImpreciseSvInfo(
                structural_type=hpotk.TermId.from_curie(
                    "SO:1000029"
                ),  # chromosomal_deletion
                variant_class=VariantClass.DEL,
                gene_id="HGNC:21316",
                gene_symbol="ANKRD11",
            ),
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="ANKRD11",
                tx_id="NM_013275.6",
                hgvs_cdna=None,
                is_preferred=True,
                variant_effects=(VariantEffect.TRANSCRIPT_ABLATION,),
                affected_exons=range(13),  # I counted 13 exons
                protein_id=None,
                hgvsp=None,
                protein_effect_coordinates=None,
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HETEROZYGOUS),
    )


@pytest.fixture(scope="package")
def patient_w_missense(
    sample_labels: SampleLabels,
    missense_variant: Variant,
) -> Patient:
    return Patient(
        labels=sample_labels,
        phenotypes=(),
        diseases=(),
        variants=(
            missense_variant,
        ),
    )


@pytest.fixture(scope="package")
def patient_w_frameshift(
    sample_labels: SampleLabels,
    frameshift_variant: Variant,
) -> Patient:
    return Patient(
        labels=sample_labels,
        phenotypes=(),
        diseases=(),
        variants=(
            frameshift_variant,
        ),
    )
