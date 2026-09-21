import hpotk
import pytest

from gpsea.analysis.clf import (
    HpoClassifier,
)
from gpsea.model import Patient, Phenotype


class TestHpoPredicate:
    @staticmethod
    def make_patient() -> Patient:
        return Patient.from_raw_parts(
            labels="TEST",
            phenotypes=(
                Phenotype.from_raw_parts("HP:0001250", is_observed=False),  # No Seizure
                Phenotype.from_raw_parts("HP:0001166", is_observed=True),  # Yes Arachnodactyly
            ),
        )

    @staticmethod
    def make_empty() -> Patient:
        return Patient.from_raw_parts(
            labels="EMPTY",
            phenotypes=(),
        )

    @pytest.mark.parametrize(
        "query,expected",
        [
            (
                "HP:0001167",  # Abnormal finger morphology
                "Yes",
            ),
            (
                "HP:0001166",  # Arachnodactyly
                "Yes",
            ),
            (
                "HP:0012638",  # Abnormal nervous system physiology
                None,
            ),
            (
                "HP:0001250",  # Seizure
                "No",
            ),
            (
                "HP:0033259",  # Non-motor seizure
                "No",
            ),
            (
                "HP:0001640",  # Cardiomegaly
                None,
            ),
        ],
    )
    def test_hpo_predicate__one_term_individual(
        self,
        hpo: hpotk.MinimalOntology[hpotk.TermId, hpotk.MinimalTerm],
        query: str,
        expected: str | None,
    ):
        predicate = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie(query),
            missing_implies_phenotype_excluded=False,
        )

        patient = TestHpoPredicate.make_patient()

        actual = predicate.test(patient)

        if expected is None:
            assert actual is None
        else:
            assert actual is not None
            assert actual.category.name == expected

    @pytest.mark.parametrize(
        "query,expected",
        [
            (
                "HP:0001167",  # Abnormal finger morphology
                "Yes",
            ),
            (
                "HP:0001166",  # Arachnodactyly
                "Yes",
            ),
            (
                "HP:0012638",  # Abnormal nervous system physiology
                "No",
            ),
            (
                "HP:0001250",  # Seizure
                "No",
            ),
            (
                "HP:0033259",  # Non-motor seizure
                "No",
            ),
            (
                "HP:0001640",  # Cardiomegaly
                "No",
            ),
        ],
    )
    def test_hpo_predicate__one_term_individual__missing_implies_excluded(
        self,
        hpo: hpotk.MinimalOntology[hpotk.TermId, hpotk.MinimalTerm],
        query: str,
        expected: str | None,
    ):
        predicate = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie(query),
            missing_implies_phenotype_excluded=True,
        )

        patient = TestHpoPredicate.make_patient()

        actual = predicate.test(patient)

        if expected is None:
            assert actual is None
        else:
            assert actual is not None
            assert actual.category.name == expected

    def test_hpo_predicate__empty(
        self,
        hpo: hpotk.MinimalOntology[hpotk.TermId, hpotk.MinimalTerm],
    ):
        """
        An individual with no terms is either assigned into no category (`None`)
        or into the "No" category if missing implies excluded.
        """
        patient = TestHpoPredicate.make_empty()

        predicate = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001250"),
            missing_implies_phenotype_excluded=False,
        )
        assert predicate.test(patient) is None

        predicate = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001250"),
            missing_implies_phenotype_excluded=True,
        )
        actual = predicate.test(patient)
        assert actual is not None
        assert actual.category.name == "No"
