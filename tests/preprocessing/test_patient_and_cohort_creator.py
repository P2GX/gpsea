import os
import io

import hpotk
from hpotk.validate import ValidationRunner
import pytest

from gpsea.model.genome import GenomeBuild
from gpsea.preprocessing import FunctionalAnnotator, ImpreciseSvFunctionalAnnotator, VariantCoordinateFinder
from gpsea.preprocessing import VVHgvsVariantCoordinateFinder, DefaultImpreciseSvFunctionalAnnotator
from gpsea.preprocessing import PhenopacketPatientCreator
from gpsea.preprocessing import VVMultiCoordinateService
from gpsea.preprocessing import CohortCreator, load_phenopacket_folder
from gpsea.preprocessing import configure_default_functional_annotator


class TestPhenopacketCohortCreator:
    @pytest.fixture
    def functional_annotator(
        self,
    ) -> FunctionalAnnotator:
        return configure_default_functional_annotator(
            ann_source="VEP",
        )

    @pytest.fixture
    def imprecise_sv_functional_annotator(
        self,
        genome_build: GenomeBuild,
    ) -> ImpreciseSvFunctionalAnnotator:
        return DefaultImpreciseSvFunctionalAnnotator(
            gene_coordinate_service=VVMultiCoordinateService(
                genome_build=genome_build,
            ),
        )

    @pytest.fixture
    def variant_coordinate_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(
            genome_build=genome_build,
        )

    @pytest.fixture
    def patient_creator(
        self,
        hpo: hpotk.MinimalOntology,
        validation_runner: ValidationRunner,
        genome_build: GenomeBuild,
        functional_annotator: FunctionalAnnotator,
        imprecise_sv_functional_annotator: ImpreciseSvFunctionalAnnotator,
        variant_coordinate_finder: VariantCoordinateFinder,
    ) -> PhenopacketPatientCreator:
        return PhenopacketPatientCreator(
            hpo=hpo,
            validator=validation_runner,
            build=genome_build,
            functional_annotator=functional_annotator,
            imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
            hgvs_coordinate_finder=variant_coordinate_finder,
        )

    @pytest.fixture
    def phenopacket_cohort_creator(
        self,
        patient_creator: PhenopacketPatientCreator,
    ) -> CohortCreator:
        return CohortCreator(
            patient_creator=patient_creator,
        )

    def test_cohort_creator(
        self,
        fpath_test_dir: str,
        phenopacket_cohort_creator: CohortCreator,
    ):
        folder = os.path.join(fpath_test_dir, "preprocessing", "data", "dup_id_test_data")
        _, results = load_phenopacket_folder(folder, phenopacket_cohort_creator)

        outfile = io.StringIO()
        results.summarize(outfile)

        actual_lines = outfile.getvalue().splitlines(keepends=False)

        expected = (
            " ·Patient ID/s Pat_1[PMID_12345], Pat_2[PMID_67890] have a duplicate. "
            "Please verify every patient has an unique ID."
        )
        assert expected in actual_lines
