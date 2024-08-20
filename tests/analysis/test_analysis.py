import pathlib
import typing

import hpotk
import numpy as np
import pandas as pd
import pytest

from genophenocorr.analysis import (
    apply_predicates_on_patients,
    CohortAnalysis, CohortAnalysisConfiguration,
    configure_cohort_analysis,
)

from genophenocorr.model import *
from genophenocorr.model.genome import *
from genophenocorr.analysis.predicate import GenotypePolyPredicate, PatientCategories
from genophenocorr.analysis.predicate.phenotype import PhenotypePolyPredicate
from genophenocorr.analysis.predicate.genotype import VariantPredicate


def test_apply_predicates_on_patients(
        suox_cohort: Cohort,
        suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
        suox_gt_predicate: GenotypePolyPredicate,
):
    categories, n_usable, counts = apply_predicates_on_patients(
        patients=suox_cohort.all_patients,
        pheno_predicates=suox_pheno_predicates,
        gt_predicate=suox_gt_predicate,
    )

    assert isinstance(categories, typing.Collection)
    assert PatientCategories.YES in categories
    assert PatientCategories.NO in categories
    assert len(categories) == 2

    assert isinstance(n_usable, pd.Series)
    assert len(n_usable) == 5

    assert isinstance(counts, typing.Mapping)
    assert len(counts) == 5

    seizure_counts = counts[hpotk.TermId.from_curie('HP:0001250')]
    assert np.array_equal(seizure_counts.to_numpy(), np.array([[17, 11], [7, 0]]))


class TestCohortAnalysis:

    @pytest.fixture
    def cohort_analysis(
        self,
        suox_cohort: Cohort,
        hpo: hpotk.MinimalOntology,
        tmp_path: pathlib.Path,
    ) -> CohortAnalysis:
        return configure_cohort_analysis(
            cohort=suox_cohort,
            hpo=hpo,
            cache_dir=str(tmp_path),
        )

    def test_analysis_passes_if_variant_predicate_always_returns_false(
        self,
        cohort_analysis: CohortAnalysis,
        always_false_variant_predicate: VariantPredicate,
    ):
        results = cohort_analysis.compare_hpo_vs_genotype(
            predicate=always_false_variant_predicate
        )
        assert results is not None

    def test_analysis_explodes_if_no_phenotypes_are_left_for_analysis(
        self,
        hpo: hpotk.MinimalOntology,
        degenerated_cohort: Cohort,
        tmp_path: pathlib.Path,
        always_false_variant_predicate: VariantPredicate,
    ):
        config = CohortAnalysisConfiguration()
        config.hpo_mtc_strategy()
        cohort_analysis = configure_cohort_analysis(
            hpo=hpo,
            cohort=degenerated_cohort,
            cache_dir=str(tmp_path),
            config=config,
        )

        with pytest.raises(ValueError) as e:
            _ = cohort_analysis.compare_hpo_vs_genotype(
                predicate=always_false_variant_predicate,
            )

        assert e.value.args[0] == 'No phenotypes are left for the analysis after MTC filtering step'
