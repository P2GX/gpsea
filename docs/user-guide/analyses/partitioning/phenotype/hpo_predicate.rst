
.. _hpo-classifier:


HPO classifier
==============

When testing for presence or absence of an HPO term, the :class:`~gpsea.analysis.clf.HpoClassifier`
leverages the :ref:`true-path-rule` to take advantage of the HPO hierarchy.
In result, an individual annotated with a term is implicitly annotated with all its ancestors.
For instance, an individual annotated with `Ectopia lentis <https://hpo.jax.org/browse/term/HP:0001083>`_
is also annotated with `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_,
`Abnormal anterior eye segment morphology <https://hpo.jax.org/browse/term/HP:0004328>`_,
`Abnormal eye morphology <https://hpo.jax.org/browse/term/HP:0012372>`_, ...

Similarly, all descendants of a term, whose presence was specifically excluded in an individual,
are implicitly excluded.


Example
-------

Here we show how to set up :class:`~gpsea.analysis.clf.HpoClassifier`
to test for a presence of `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_.

We need to load :class:`~hpotk.ontology.MinimalOntology` with HPO data to access the HPO hierarchy:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

and now we can set up the classifier to test for presence of *Abnormal lens morphology*:

>>> from gpsea.analysis.clf import HpoClassifier
>>> query = hpotk.TermId.from_curie('HP:0000517')
>>> pheno_clf = HpoClassifier(
...     hpo=hpo,
...     query=query,
... )
>>> pheno_clf.name
'HPO Classifier'
>>> pheno_clf.description
'Test for presence of Abnormal lens morphology [HP:0000517]'
>>> pheno_clf.class_labels
('Yes', 'No')



missing_implies_phenotype_excluded
----------------------------------

In many cases, published reports of clinical data about individuals with rare diseases
describe phenotypic features that were observed, but do not provide
a comprehensive list of features that were explicitly excluded.
By default, GPSEA will only include features that are recorded as observed or excluded in a phenopacket.

However, setting ``missing_implies_excluded=True`` will cause "n/a" entries to be set to "excluded".
We provide this option for exploration but do not recommend its use
for the final analysis unless the assumption behind it is known to be true.



Classifiers for all cohort phenotypes
=====================================

Constructing phenotype classifiers for all HPO terms of a cohort sounds a bit tedious.
The :func:`~gpsea.analysis.clf.prepare_classifiers_for_terms_of_interest`
function cuts down the tedium.


Example
-------

For a phenopacket collection (e.g. 156 patients with mutations in *TBX5* gene included in Phenopacket Store version `0.1.18`)

>>> from ppktstore.registry import configure_phenopacket_registry
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store(release='0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets('TBX5'))
>>> len(phenopackets)
156

processed into a cohort

>>> from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets
>>> cohort_creator = configure_caching_cohort_creator(hpo)
>>> cohort, _ = load_phenopackets(phenopackets, cohort_creator)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Individuals Processed: ...


we can create HPO classifiers for testing all 369 HPO terms used in the cohort:

>>> from gpsea.analysis.clf import prepare_classifiers_for_terms_of_interest
>>> pheno_clfs = prepare_classifiers_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
... )
>>> len(pheno_clfs)
369

and subject the predicates into further analysis, such as :class:`~gpsea.analysis.pcats.HpoTermAnalysis`.
