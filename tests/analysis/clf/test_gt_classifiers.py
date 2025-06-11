import typing
import pytest

from gpsea.analysis.clf import frozen_classifier, GenotypeClassifier
from gpsea.model import Patient, Cohort


class TestFrozenClassifier:
    @pytest.fixture(scope="class")
    def gt_classifier(
        self,
        genesis_family: typing.Sequence[Patient],
    ) -> GenotypeClassifier:
        return frozen_classifier(
            samples=genesis_family,
            codes=(0, 0, 1),
            labels={
                0: "Saw paradise",
                1: "Did not see paradise",
            },
        )

    def test_classifier_properties(
        self,
        gt_classifier: GenotypeClassifier,
    ):
        assert all(
            cl in gt_classifier.class_labels
            for cl in ("Saw paradise", "Did not see paradise")
        )
        assert len(gt_classifier.class_labels) == 2

        assert gt_classifier.name == "Frozen genotype classifier"
        assert (
            gt_classifier.description
            == "Classify the individuals based on provided genotype codes"
        )
        assert gt_classifier.variable_name == "Genotype code"

    def test_test(
        self,
        gt_classifier: GenotypeClassifier,
        adam: Patient,
        eve: Patient,
        cain: Patient,
    ):
        cat = gt_classifier.test(adam)
        assert cat is not None
        assert cat.category.name == "Saw paradise"

        cat = gt_classifier.test(eve)
        assert cat is not None
        assert cat.category.name == "Saw paradise"

        cat = gt_classifier.test(cain)
        assert cat is not None
        assert cat.category.name == "Did not see paradise"

    def test_test_fails_on_unknown_individual(
        self,
        gt_classifier: GenotypeClassifier,
        walt: Patient,
    ):
        with pytest.raises(ValueError) as e:
            gt_classifier.test(walt)

        assert e.value.args == ("Unexpected patient Walt",)

    def test_can_be_created_from_a_cohort(
        self,
        genesis_family: typing.Sequence[Patient],
    ):
        cohort = Cohort.from_patients(
            members=genesis_family,
        )
        
        gt_clf = frozen_classifier(
            samples=cohort,
            codes=(0, 0, 1),
        )
        
        assert isinstance(gt_clf, GenotypeClassifier)

    def test_reports_code_count_mismatch(self):
        with pytest.raises(AssertionError) as e:
            frozen_classifier(
                samples=(),
                codes=(0,),
            )
        
        assert e.value.args == ("Sample count 0 must match the code count 1",)
