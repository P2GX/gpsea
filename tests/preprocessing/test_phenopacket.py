import os

import hpotk
import pytest

from hpotk.validate import ValidationRunner
from stairval import Level
from stairval.notepad import create_notepad

from google.protobuf.json_format import Parse
from phenopackets.schema.v2.core.interpretation_pb2 import GenomicInterpretation
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket

from gpsea.model import VariantClass
from gpsea.model.genome import GenomeBuild, Strand

from gpsea.preprocessing import VVHgvsVariantCoordinateFinder
from gpsea.preprocessing import (
    VariantCoordinateFinder,
    PhenopacketVariantCoordinateFinder,
)
from gpsea.preprocessing import (
    FunctionalAnnotator,
    ImpreciseSvFunctionalAnnotator,
    DefaultImpreciseSvFunctionalAnnotator,
    configure_default_functional_annotator,
)
from gpsea.preprocessing import PhenopacketPatientCreator, PhenopacketOntologyTermOnsetParser
from gpsea.preprocessing import VVMultiCoordinateService


class TestPhenopacketVariantCoordinateFinder:
    @pytest.fixture(scope="class")
    def fpath_test_genomic_interpretations(
        self,
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, "pp_genomic_interpretations")

    @pytest.fixture(scope="class")
    def hgvs_vc_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(genome_build)

    @pytest.fixture(scope="class")
    def pp_vc_finder(
        self,
        genome_build: GenomeBuild,
        hgvs_vc_finder: VariantCoordinateFinder,
    ) -> PhenopacketVariantCoordinateFinder:
        return PhenopacketVariantCoordinateFinder(genome_build, hgvs_vc_finder)

    @pytest.mark.online
    @pytest.mark.parametrize(
        "pp_name, contig, start, end, ref, alt, change_length",
        [
            (
                "deletion_test.json",
                "16",
                89284128,
                89284134,
                "CTTTTT",
                "C",
                -5,
            ),
            (
                "insertion_test.json",
                "16",
                89280828,
                89280829,
                "C",
                "CA",
                1,
            ),
            (
                "missense_test.json",
                "16",
                89279134,
                89279135,
                "G",
                "C",
                0,
            ),
            (
                "missense_hgvs_test.json",
                "16",
                89279134,
                89279135,
                "G",
                "C",
                0,
            ),
            (
                "duplication_test.json",
                "16",
                89279849,
                89279850,
                "G",
                "GC",
                1,
            ),
            (
                "delinsert_test.json",
                "16",
                89284600,
                89284602,
                "GG",
                "A",
                -1,
            ),
            (
                "CVDup_test.json",
                "16",
                89_284_522,
                89_373_231,
                "N",
                "<DUP>",
                88_709,
            ),
            (
                "CVDel_test.json",
                "16",
                89_217_280,
                89_506_042,
                "N",
                "<DEL>",
                -288_762,
            ),
        ],
    )
    def test_find_coordinates(
        self,
        pp_name: str,
        contig: str,
        start: int,
        end: int,
        ref: str,
        alt: str,
        change_length: int,
        fpath_test_genomic_interpretations: str,
        pp_vc_finder: PhenopacketVariantCoordinateFinder,
    ):
        fpath_pp = os.path.join(fpath_test_genomic_interpretations, pp_name)
        gi = read_genomic_interpretation_json(fpath_pp)

        vc = pp_vc_finder.find_coordinates(gi)

        assert vc is not None

        assert vc.chrom == contig
        assert vc.start == start
        assert vc.end == end
        assert vc.region.strand == Strand.POSITIVE
        assert vc.ref == ref
        assert vc.alt == alt
        assert vc.change_length == change_length

    def test_find_large_structural(
        self,
        fpath_test_genomic_interpretations: str,
        pp_vc_finder: PhenopacketVariantCoordinateFinder,
    ):
        fpath_pp = os.path.join(fpath_test_genomic_interpretations, "chromosomal_deletion.ANKRD11.json")
        gi = read_genomic_interpretation_json(fpath_pp)

        vc = pp_vc_finder.find_coordinates(gi)
        assert vc is None


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())


class TestPhenopacketPatientCreator:
    @pytest.fixture(scope="class")
    def functional_annotator(
        self,
        fpath_project_dir: str,
    ) -> FunctionalAnnotator:
        fpath_cache_dir = os.path.join(fpath_project_dir, ".gpsea_cache")
        fpath_variant_cache_dir = os.path.join(fpath_cache_dir, "variant_cache")
        os.makedirs(fpath_variant_cache_dir, exist_ok=True)

        return configure_default_functional_annotator(
            ann_source="VEP",
            cache_dir=fpath_variant_cache_dir,
        )

    @pytest.fixture(scope="class")
    def imprecise_sv_functional_annotator(
        self,
        genome_build: GenomeBuild,
    ) -> ImpreciseSvFunctionalAnnotator:
        return DefaultImpreciseSvFunctionalAnnotator(
            gene_coordinate_service=VVMultiCoordinateService(
                genome_build=genome_build,
            ),
        )

    @pytest.fixture(scope="class")
    def variant_coordinate_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(
            genome_build=genome_build,
        )

    @pytest.fixture(scope="class")
    def onset_term_parser(self) -> PhenopacketOntologyTermOnsetParser:
        return PhenopacketOntologyTermOnsetParser.default_parser()

    @pytest.fixture(scope="class")
    def patient_creator(
        self,
        hpo: hpotk.MinimalOntology,
        validation_runner: ValidationRunner,
        genome_build: GenomeBuild,
        functional_annotator: FunctionalAnnotator,
        imprecise_sv_functional_annotator: ImpreciseSvFunctionalAnnotator,
        variant_coordinate_finder: VariantCoordinateFinder,
        onset_term_parser: PhenopacketOntologyTermOnsetParser,
    ) -> PhenopacketPatientCreator:
        return PhenopacketPatientCreator(
            hpo=hpo,
            validator=validation_runner,
            build=genome_build,
            functional_annotator=functional_annotator,
            imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
            hgvs_coordinate_finder=variant_coordinate_finder,
            term_onset_parser=onset_term_parser,
        )

    @pytest.fixture(scope="class")
    def phenopacket(
        self,
        fpath_phenopacket_dir: str,
    ) -> Phenopacket:
        fpath_pp = os.path.join(fpath_phenopacket_dir, "PMID_30968594_individual_1.json")

        with open(fpath_pp) as fh:
            return Parse(fh.read(), Phenopacket())

    def test_phenopacket_patient_creator(
        self,
        phenopacket: Phenopacket,
        patient_creator: PhenopacketPatientCreator,
    ):
        notepad = create_notepad("A phenopacket")
        patient = patient_creator.process(phenopacket, notepad)

        assert patient is not None

        # No issues
        assert not notepad.has_errors_or_warnings(include_subsections=True)

        # Individual credentials are OK
        assert patient.labels.label_summary() == "individual 1[PMID_30968594_individual_1]"
        assert patient.sex.is_male()
        assert patient.age is not None
        assert patient.age.days == pytest.approx(334.8125)
        assert patient.age.is_postnatal

        assert patient.vital_status is not None
        assert patient.vital_status.is_deceased
        assert patient.vital_status.age_of_death is not None
        assert patient.vital_status.age_of_death.days == pytest.approx(365.25)
        assert patient.vital_status.age_of_death.is_postnatal

        # 6 Phenotype features
        assert len(patient.phenotypes) == 6
        assert tuple(p.identifier.value for p in patient.phenotypes) == (
            "HP:0000953",
            "HP:0030088",
            "HP:0003154",
            "HP:0008163",
            "HP:0000870",
            "HP:0025133",
        )
        assert tuple(p.is_present for p in patient.phenotypes) == (
            True,
            True,
            True,
            True,
            True,
            False,
        )
        # Check onset of Hyperpigmentation of the skin `HP:0000953`
        hyperpigmentation_of_the_skin = patient.phenotype_by_id("HP:0000953")
        assert hyperpigmentation_of_the_skin is not None
        # Expecting congenital onset
        onset = hyperpigmentation_of_the_skin.onset
        assert onset is not None
        assert onset.days == pytest.approx(0.0)
        assert onset.is_postnatal

        # 6 measurements
        assert len(patient.measurements) == 6
        assert tuple(m.identifier.value for m in patient.measurements) == (
            "LOINC:1668-3",
            "LOINC:2986-8",
            "LOINC:2141-0",
            "LOINC:2143-6",
            "LOINC:2842-3",
            "LOINC:2243-4",
        )

        units = tuple(m.unit.value for m in patient.measurements)
        assert units == (
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
        )

        values = tuple(m.test_result for m in patient.measurements)
        assert values == (800.0, 127.0, 180.2, 116.6, 52.93, 23.71)

        # a disease
        assert len(patient.diseases) == 1

        disease = patient.diseases[0]
        assert disease.identifier.value == "OMIM:201910"
        assert disease.is_present is True

        assert disease.onset is not None
        assert disease.onset.is_postnatal is True
        assert disease.onset.days == pytest.approx(20.0)

        # variants
        assert len(patient.variants) == 3
        snp = patient.variants[0]
        assert snp.variant_info.has_variant_coordinates()
        snp_vc = snp.variant_info.variant_coordinates
        assert snp_vc is not None
        assert snp_vc.chrom == "6"
        assert snp_vc.start == 32_040_420
        assert snp_vc.end == 32_040_421
        assert snp_vc.ref == "C"
        assert snp_vc.alt == "T"
        assert snp_vc.change_length == 0
        assert snp_vc.variant_class == VariantClass.SNV

        imprecise_sv = patient.variants[1]
        assert imprecise_sv.variant_info.has_sv_info()
        sv_vi = imprecise_sv.variant_info.sv_info
        assert sv_vi is not None
        assert sv_vi.gene_id == "HGNC:2600"
        assert sv_vi.gene_symbol == "CYP21A2"
        assert sv_vi.structural_type.value == "SO:1000029"  # `chromosomal_deletion`
        assert sv_vi.variant_class == VariantClass.DEL

        imprecise_tra = patient.variants[2]
        assert imprecise_tra.variant_info.has_sv_info()
        tra_vi = imprecise_tra.variant_info.sv_info
        assert tra_vi is not None
        assert tra_vi.gene_id == "HGNC:24650"
        assert tra_vi.gene_symbol == "EHMT1"
        assert tra_vi.structural_type.value == "SO:1000044"  # `chromosomal_translocation`
        assert tra_vi.variant_class == VariantClass.TRANSLOCATION

    def test_individual_with_no_genotype(
        self,
        phenopacket: Phenopacket,
        patient_creator: PhenopacketPatientCreator,
    ):
        """
        Check if we get an error if no genotype info is available.
        """
        pp = Phenopacket()
        pp.CopyFrom(phenopacket)
        del pp.interpretations[:]  # clear variants

        notepad = create_notepad("no-gt")

        _ = patient_creator.process(pp=pp, notepad=notepad)

        errors = tuple(notepad.errors())

        assert len(errors) == 1
        error = errors[0]
        assert error.level == Level.ERROR
        assert error.message == "Individual PMID_30968594_individual_1 has no genotype data (variants) to work with"
        assert error.solution == "Add variants or remove the individual from the analysis"

    def test_individual_with_no_phenotype(
        self,
        phenopacket: Phenopacket,
        patient_creator: PhenopacketPatientCreator,
    ):
        """
        Check if we get an error if no phenotype info is available.
        """
        pp = Phenopacket()
        pp.CopyFrom(phenopacket)
        # clear phenotype
        del pp.phenotypic_features[:]
        del pp.diseases[:]
        del pp.measurements[:]

        notepad = create_notepad("no-gt")

        _ = patient_creator.process(pp=pp, notepad=notepad)

        errors = tuple(notepad.errors())

        assert len(errors) == 1
        error = errors[0]
        assert error.level == Level.ERROR
        assert (
            error.message
            == "Individual PMID_30968594_individual_1 has no phenotype data (HPO, a diagnosis, measurement) to work with"
        )
        assert (
            error.solution == "Add HPO terms, a diagnosis, or measurements, or remove the individual from the analysis"
        )
