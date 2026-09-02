import hpotk
import pytest

from gpsea.analysis.clf import DiseasePresenceClassifier
from gpsea.model import Cohort, Patient


class TestDiseasePresencePredicate:
    @pytest.mark.parametrize(
        "patient_id, patient_category",
        [
            ("HetSingleVar", "Yes"),
            ("HomoVar", "No"),
        ],
    )
    def test_disease_predicate(
        self,
        patient_id: str,
        patient_category: str,
        toy_cohort: Cohort,
    ):
        patient = find_patient(patient_id, toy_cohort)
        disease_id = hpotk.TermId.from_curie("OMIM:148050")
        predicate = DiseasePresenceClassifier(disease_id)
        actual = predicate.test(patient)

        assert actual is not None
        assert actual.phenotype == disease_id
        assert actual.category.name == patient_category


def find_patient(pat_id: str, cohort: Cohort) -> Patient:
    for pat in cohort.all_patients:
        if pat.patient_id == pat_id:
            return pat
    raise ValueError(f"Could not find patient {pat_id}")
