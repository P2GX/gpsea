import abc
import typing

import hpotk

from genophenocorr.model import Patient

from .._api import PolyPredicate, PatientCategory, PatientCategories, Categorization, IntegerCountCategory

P = typing.TypeVar('P')
"""
Phenotype entity of interest, such as :class:`hpotk.model.TermId`, representing an HPO term or an OMIM/MONDO term.

However, phenotype entity can be anything as long as it is :class:`typing.Hashable` and comparable
(have `__eq__` and `__lt__` magic methods).
"""


class PhenotypeCategorization(typing.Generic[P], Categorization):
    """"
    On top of the attributes of the `Categorization`, `PhenotypeCategorization` keeps track of the target phenotype `P`.
    """

    def __init__(
            self,
            category: PatientCategory,
            phenotype: P,
    ):
        super().__init__(category)
        self._phenotype = phenotype

    @property
    def phenotype(self) -> P:
        return self._phenotype

    def __repr__(self) -> str:
        return f"PhenotypeCategorization(" \
               f"category={self._category}, " \
               f"phenotype={self._phenotype})"

    def __str__(self) -> str:
        return repr(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, PhenotypeCategorization) \
            and self._category == other._category \
            and self._phenotype == other._phenotype

    def __hash__(self) -> int:
        return hash((self._category, self._phenotype))


class PhenotypePolyPredicate(typing.Generic[P], PolyPredicate[PhenotypeCategorization[P]], metaclass=abc.ABCMeta):
    """
    `PhenotypePolyPredicate` investigates a patient in context of one or more phenotype categories `P`.

    Each phenotype category `P` can be a :class:`hpotk.model.TermId` representing an HPO term or an OMIM/MONDO term.

    Most of the time, only one category is investigated, and :attr:`phenotype` will return a sequence with
    one element only (e.g. *Arachnodactyly* `HP:0001166`). However, multiple categories can be sought as well,
    such as when the predicate bins the patient into one of discrete diseases
    (e.g. Geleophysic dysplasia, Marfan syndrome, ...). Then the predicate should return a sequence of disease
    identifiers.
    """

    @property
    @abc.abstractmethod
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the phenotype entities of interest.
        """
        pass


class PropagatingPhenotypePredicate(PhenotypePolyPredicate[hpotk.TermId]):
    """
    `PropagatingPhenotypePredicate` tests if a patient is annotated with an HPO term.

    Note, `query` must be a term of the provided `hpo`!

    :param hpo: HPO object
    :param query: the HPO term to test
    :param missing_implies_phenotype_excluded: `True` if lack of an explicit annotation implies term's absence`.
    """

    def __init__(self, hpo: hpotk.MinimalOntology,
                 query: hpotk.TermId,
                 missing_implies_phenotype_excluded: bool = False):
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._query = hpotk.util.validate_instance(query, hpotk.TermId, 'phenotypic_feature')
        self._query_label = self._hpo.get_term(query)
        assert self._query_label is not None, f'Query {query} is in HPO'
        self._missing_implies_phenotype_excluded = hpotk.util.validate_instance(missing_implies_phenotype_excluded, bool,
                                                                      'missing_implies_phenotype_excluded')
        self._phenotype_observed = PhenotypeCategorization(
            category=PatientCategories.YES,
            phenotype=self._query,
        )
        self._phenotype_excluded = PhenotypeCategorization(
            category=PatientCategories.NO,
            phenotype=self._query,
        )

    def get_question(self) -> str:
        return f'Is {self._query_label} present in the patient?'

    @property
    def phenotypes(self) -> typing.Sequence[hpotk.TermId]:
        # We usually test just a single HPO term, so we return a tuple with a single member.
        return self._query,

    def get_categorizations(self) -> typing.Sequence[PhenotypeCategorization[hpotk.TermId]]:
        return self._phenotype_observed, self._phenotype_excluded

    def test(self, patient: Patient) -> typing.Optional[PhenotypeCategorization[hpotk.TermId]]:
        """An HPO TermID is given when initializing the class.
        Given a Patient class, this function tests whether the patient has the
        given phenotype.

        Args:
            patient (Patient): A Patient class representing a patient.

        Returns:
            typing.Optional[PhenotypeCategorization[P]]: PhenotypeCategorization,
                                                        either "YES" if the phenotype
                                                        is listed and is not excluded, or
                                                        "NO" if the phenotype is listed and excluded,
                                                        otherwise will return None.
                                                        Unless _missing_implies_phenotype_excluded is True, then
                                                        will return "NO" if the phenotype is listed and excluded
                                                        or not listed.
        """
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_observed:
                if any(self._query == anc for anc in self._hpo.graph.get_ancestors(phenotype, include_source=True)):
                    return self._phenotype_observed
            else:
                if self._missing_implies_phenotype_excluded:
                    return self._phenotype_excluded
                else:
                    if any(self._query == desc for desc in
                           self._hpo.graph.get_descendants(phenotype, include_source=True)):
                        return self._phenotype_excluded

        return None

    def __repr__(self):
        return f'PropagatingPhenotypeBooleanPredicate(query={self._query})'


class DiseasePresencePredicate(PhenotypePolyPredicate[hpotk.TermId]):
    """
    `DiseasePresencePredicate` tests if the patient was diagnosed with a disease.

    The predicate tests if the patient's diseases include a disease ID formatted as a :class:`hpotk.model.TermId`.

    :param disease_id_query: the Disease ID to test
    """

    def __init__(self, disease_id_query: hpotk.TermId):
        self._query = hpotk.util.validate_instance(disease_id_query, hpotk.TermId, 'disease_id_query')

        self._diagnosis_present = PhenotypeCategorization(
            category=PatientCategories.YES,
            phenotype=disease_id_query,
        )
        self._diagnosis_excluded = PhenotypeCategorization(
            category=PatientCategories.NO,
            phenotype=disease_id_query,
        )

    def get_question(self) -> str:
        return f'Was {self._query} diagnosed in the patient?'

    @property
    def phenotypes(self) -> typing.Sequence[hpotk.TermId]:
        # We usually test just a single Disease, so we return a tuple with a single member.
        return self._query,

    def get_categorizations(self) -> typing.Sequence[PhenotypeCategorization[hpotk.TermId]]:
        return self._diagnosis_present, self._diagnosis_excluded

    def test(self, patient: Patient) -> typing.Optional[PhenotypeCategorization[hpotk.TermId]]:
        """
        Test if the `patient` was diagnosed with a disease.

        Args:
            patient (Patient): a `patient` instance to be tested.

        Returns:
            typing.Optional[PhenotypeCategorization[P]]: PhenotypeCategorization,
                                                        either "YES" if the phenotype
                                                        is listed and is not excluded, or
                                                        "NO" if the disease is not listed
                                                        or is excluded.
        """
        self._check_patient(patient)

        for dis in patient.diseases:
            if dis.is_present and dis.identifier == self._query:
                return self._diagnosis_present

        return self._diagnosis_excluded

    def __repr__(self):
        return f'DiseasePresencePredicate(query={self._query})'


class CountingPhenotypePredicate(PhenotypePolyPredicate[int]):
    """
    `CountingPhenotypePredicate` bins the patient
    according to the count of present phenotypes that are either
    an exact match to the `query` terms or their descendants. For instance, we may want to count whether an individual has
    brain, liver, kidney, and skin abormalities. In the case, the query would include the corresponding terms (e.g.,
    Abnormal brain morphology HP:0012443). An individual can then have between 0 and 4. This predicate is intended
    to be used with the Mann Whitney U test.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        query: typing.Iterable[typing.Union[str, hpotk.TermId]],
    ):
        self._hpo = hpo
        # `query` is an iterable of either CURIEs (e.g. `HP:0001250`) or `TermId`s.
        # We take curies as a convenience to the user. Ensure we map a curie to `TermId` using `TermId.from_curie`.
        self._query = set()
        for q in query:
            if isinstance(q,str):
                q = hpotk.TermId.from_curie(q)
            if not isinstance(q, hpotk.TermId):
                raise ValueError(f"query argument must be iterable of hpotk TermId's or strings but we found {type(q)}")
            if not q in self._hpo:
                raise ValueError(f"The query {q} was not found in the HPO")
            self._query.add(q)
        # the query terms must not include a term and its ancestor
        for q in self._query:
            ancs = self._hpo.graph.get_ancestors(q, include_source=False)
            for anc in ancs:
                if anc in self._query:
                    raise ValueError(f"Both {q} and its ancestor term {anc} were found in the query, but query terms must not include a term and its ancestor")
        self._max_terms = len(self._query)

    def test(self, patient: Patient) -> typing.Optional[PhenotypeCategorization[int]]:
        """
        return the count (number) of terms in the target set (self._query) that have matching terms (exact matches or descendants) in the patient.
        Do not double count if the patient has two terms (e.g., two different descendants) of one of the query terms.
        """
        self._check_patient(patient)
        matching_termid_set = set() ## use a set to avoid double-counting. We add terms from self._query that are matched
        for pf in patient.present_phenotypes():
            hpo_id = pf.identifier
            ancs = self._hpo.graph.get_ancestors(hpo_id, include_source=True)
            for anc in ancs:
                if anc in self._query:
                    matching_termid_set.add(anc)
        ## Note the categories are only needed for 2x2 or 2x3 tests. They
        ## are not needed for Mann Whitney U test. The following is just to satisfy the API.
        countCategorization = PatientCategory(cat_id=424242, name="counting")
        return PhenotypeCategorization(phenotype=len(matching_termid_set), category=countCategorization)

    @property
    def phenotypes(self) -> typing.Sequence[P]:
        """
        The phenotype ranges from 0 (none of the target HPO ids in self._query match) to self._max_terms (all of the target HPO ids match)
        We return a list of corresponding integers, e.g.,  [0, .. , 5]
        """
        return list(range(self._max_terms + 1))

    def get_question(self) -> str:
        # What does this predicate answer?
        return "how many of the target HPO terms (or their descendants) are does the individual display?"

    def get_categorizations(self) -> typing.Sequence[PhenotypeCategorization[int]]:
        # e.g. [0, .. , 5]
        # This should return descriptions of all categories that the predicate can produce.
        return [PhenotypeCategorization(IntegerCountCategory(p), p) for p in self.phenotypes]
