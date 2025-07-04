import abc
import dataclasses
import typing

from stairval.notepad import Notepad

from gpsea.model import Patient, Cohort

T = typing.TypeVar("T")
"""
The input for `PatientCreator`.

It can be any object that contains the patient data (e.g. a phenopacket).
"""


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    """
    `PatientCreator` can create a `Patient` from some input `T`.
    """

    @abc.abstractmethod
    def process(
        self,
        item: T,
        notepad: Notepad,
    ) -> typing.Optional[Patient]:
        pass


@dataclasses.dataclass()
class CohortCreatorOptions:
    """
    Options for :class:`~gpsea.preprocessing.CohortCreator`.
    """

    keep_individuals_with_no_hpo: bool = False
    keep_individuals_with_no_variants: bool = False


class CohortCreator(typing.Generic[T]):
    """
    `CohortCreator` creates a cohort from the provided `inputs`,
    subjecting the cohort members to Q/C and filtering.

    Cohort creator uses :class:`~gpsea.preprocessing.PatientCreator`
    to map each cohort member into a :class:`~gpsea.model.Patient`.
    The cohort creator is generic over the cohort member type `T`
    and all that matters is if the inner :class:`~gpsea.preprocessing.PatientCreator`
    can map `T` into a :class:`~gpsea.model.Patient`.


    **Q/C checks**

    The members are checked for duplicates and the duplicates are reported into the `notepad`.
    The `notepad` also retains the issues found by :class:`~gpsea.preprocessing.PatientCreator`.

    Note, a cohort *is* created even in presence of Q/C errors.


    **Filtering**

    The following filters are applied after mapping `T` to cohort members:

    * filter out the individuals who have 0 phenotypes, controlled by :class:`gpsea.preprocessing.CohortCreatorOptions.keep_individuals_with_no_hpo`
    * filter out the individuals who have 0 variants, controlled by :class:`gpsea.preprocessing.CohortCreatorOptions.keep_individuals_with_no_variants`


    **Cohort member order**

    Cohort creator guarantees stable order of the cohort members, i.e. iterating over :meth:`~gpsea.model.Cohort`
    yields the cohort members in the same order as they were seen in the `inputs` iterable.

    :param patient_creator: an instance of :class:`~gpsea.preprocessing.PatientCreator` to map `T` into :class:`~gpsea.model.Patient`.
    :param options: cohort creator options or `None` if default options should be used.
    """

    def __init__(
        self,
        patient_creator: PatientCreator[T],
        options: typing.Optional[CohortCreatorOptions] = None,
    ):
        # Check that we're getting a `PatientCreator`.
        # Unfortunately, we cannot check that `T`s of `PatientCreator` and `CohortCreator` actually match
        # due to Python's loosey-goosey nature.
        assert isinstance(patient_creator, PatientCreator)
        self._pc = patient_creator
        if options is None:
            self._options = CohortCreatorOptions()
        else:
            assert isinstance(options, CohortCreatorOptions)
            self._options = options

    def process(
        self,
        inputs: typing.Iterable[T],
        notepad: Notepad,
    ) -> Cohort:
        """
        Process the `inputs` into a :class:`~gpsea.model.Cohort`
        and write any Q/C issues into the `notepad`.
        """
        patients = []
        patient_labels = set()
        duplicate_pat_labels = set()

        for i, pp in enumerate(inputs):
            sub = notepad.add_subsection(f"patient #{i}")
            patient = self._pc.process(pp, sub)
            if patient is not None:
                if patient.labels in patient_labels:
                    duplicate_pat_labels.add(patient.labels)
                patient_labels.add(patient.labels)
                patients.append(patient)

        if len(duplicate_pat_labels) > 0:
            label_summaries = [d.label_summary() for d in duplicate_pat_labels]
            label_summaries.sort()
            notepad.add_error(
                f"Patient ID/s {', '.join(label_summaries)} have a duplicate",
                "Please verify every patient has an unique ID.",
            )

        return Cohort.from_patients(
            members=filter(self._passes_filtering, patients),
        )

    def _passes_filtering(
        self,
        patient: Patient,
    ) -> bool:
        if len(patient.phenotypes) == 0 and not self._options.keep_individuals_with_no_hpo:
            return False
        elif len(patient.variants) == 0 and not self._options.keep_individuals_with_no_variants:
            return False
        else:
            return True
