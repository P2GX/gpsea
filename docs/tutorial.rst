.. _tutorial:


########
Tutorial
########

.. doctest::
  :hide:

  >>> from gpsea import _overwrite

Here we present an example genotype-phenotype (G/P) analysis with GPSEA.
We assume GPSEA was installed as described in the :ref:`setup` section
and we encourage the users to execute the tutorial code in a Jupyter notebook
or a similar interactive environment.


**********
Background
**********

In this tutorial, we analyze a cohort of individuals with pathogenic variants in *TBX5* leading to
`Holt-Oram syndrome MIM:142900 <https://omim.org/entry/142900>`_.

Holt-Oram syndrome is an autosomal dominant disorder characterized by
upper limb defects, congenital heart defects, and arrhythmias (`PMID:38336121 <https://pubmed.ncbi.nlm.nih.gov/38336121/>`_).
It has been observed in the literature that congenital defects of the ventricular and atrial septum
are more common in the truncating than in the missense variants (`PMID:30552424 <https://pubmed.ncbi.nlm.nih.gov/30552424/>`_).
Additionally, upper limb defects are more frequent in patients with protein-truncating variants (`PMID:38336121 <https://pubmed.ncbi.nlm.nih.gov/38336121/>`_).

We curated the literature and created a `GA4GH phenopacket <https://pubmed.ncbi.nlm.nih.gov/35705716/>`_
for each affected individual.
The phenopackets are published in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.


********
Analysis
********

A typical GPSEA analysis will consist of several steps. Starting with a collection of phenopackets,
we perform input Q/C and functional variant annotation to prepare a cohort.
With the cohort on hand, we generate reports with summary statistics, variant distributions,
and most common Human Phenotype Ontology (HPO) terms or measurements.
We then configure the methods for partitioning the cohort into genotype and phenotype classes,
to test for possible associations between the classes.
We finalize the analysis by statistical testing and evaluation of the results.


.. contents:: GPSEA analysis workflow
  :depth: 2
  :local:


Inputs and Q/C
==============


Gene and transcript
-------------------

For the analysis, the `MANE <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ transcript
(i.e., the "main" biomedically relevant transcript of a gene) should be chosen unless
there is a specific reason not to (which should occur rarely if at all).

In the case of *TBX5* the MANE transcript is `NM_181486.4`. Note that the trascript identifier (`NM_181486`) and the version (`4`) are both required.
A good way to find the MANE transcript is to search on the gene symbol (e.g., *TBX5*) in `ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`_ and to
choose a variant that is specifically located in the gene. The MANE transcript will be displayed here (e.g., `NM_181486.4(TBX5):c.1221C>G (p.Tyr407Ter)
<https://www.ncbi.nlm.nih.gov/clinvar/variation/495227/>`_).

We additionally need the corresponding protein identifier.
A good way to find this is to search on the transcript id in `NCBI Nucleotide <https://www.ncbi.nlm.nih.gov/nuccore/>`_.
In our case, search on `NM_181486.4` will bring us to `this page <https://www.ncbi.nlm.nih.gov/nuccore/NM_181486.4>`_.
If we search within this page for `"NP_"`, this will bring us to the
corresponding protein accession `NP_852259.1`.

>>> cohort_name = 'TBX5'
>>> tx_id = 'NM_181486.4'
>>> px_id = 'NP_852259.1'


Human Phenotype Ontology
------------------------

The analysis in this tutorial needs to access Human Phenotype Ontology (HPO).
We use `HPO toolkit <https://ielis.github.io/hpo-toolkit/stable/>`_
to load the version `v2024-07-01` of HPO:


>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

.. tip::

  Use the latest HPO release by omitting the `release` option from the loader method.


Load phenopackets
-----------------

Now we will load the samples to analyze. We will use the cohort of 156 individuals with mutations in *TBX5*
whose clinical signs and symptoms were encoded into HPO terms
and stored in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.

>>> from ppktstore.registry import configure_phenopacket_registry
>>> phenopacket_registry = configure_phenopacket_registry()
>>> with phenopacket_registry.open_phenopacket_store("0.1.20") as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets(cohort_name))
>>> len(phenopackets)
156

We loaded 156 phenopackets which need further preprocessing to prepare for the analysis.
We will compute functional annotations for the mutations and then include the individuals into
a :class:`~gpsea.model.Cohort`:

>>> from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets
>>> cohort_creator = configure_caching_cohort_creator(hpo)
>>> cohort, validation = load_phenopackets(  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
...     phenopackets=phenopackets,
...     cohort_creator=cohort_creator,
... )
Individuals Processed: ...

and we will check that there are no Q/C issues:

>>> validation.summarize()  # doctest: +SKIP
Validated under none policy
No errors or warnings were found

We loaded the patient data into a `cohort` which is ready for the next steps.

.. seealso::

  Here we show how to create a :class:`~gpsea.model.Cohort` from phenopackets.
  See :ref:`input-data` section to learn how to create a cohort from another inputs.


Explore cohort
==============

Once the genotype and phenotype has been standardized, we can generate reports
to gain insight for the cohort data.


Show cohort summary
-------------------

The cohort summary report provides an overview about
the most common HPO terms, variants, diseases, and variant effects:

>>> from gpsea.view import CohortViewer
>>> viewer = CohortViewer(hpo)
>>> report = viewer.process(cohort=cohort, transcript_id=tx_id)
>>> report  # doctest: +SKIP

.. raw:: html
  :file: report/tbx5_cohort_info.html

.. doctest::
  :hide:

  >>> if _overwrite: report.write('docs/report/tbx5_cohort_info.html')


Plot distribution of variants with respect to the protein sequence
------------------------------------------------------------------

We can use :class:`~gpsea.view.CohortArtist` to plot the distribution of variants
with respect to the encoded protein on
a Matplotlib `Axes <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.html>`_:

>>> import matplotlib.pyplot as plt
>>> from gpsea.view import configure_default_cohort_artist
>>> cohort_artist = configure_default_cohort_artist()
>>> fig, ax = plt.subplots(figsize=(15, 8))
>>> cohort_artist.draw_protein(
...     cohort=cohort,
...     protein_id=px_id,
...     ax=ax,
... )

.. image:: /img/tutorial/tbx5_protein_diagram.png
   :alt: TBX5 protein diagram
   :align: center
   :width: 600px

.. doctest::
  :hide:

  >>> if _overwrite: fig.tight_layout(); fig.savefig('docs/img/tutorial/tbx5_protein_diagram.png')

The diagram plots the location of the variants with respect to the protein sequence.
The variant location is represented by a "lollipop".
The lollipop color represents the predicted variant effect and the lollipop size corresponds to the allele count within the cohort.
The diagram also highlights the protein features (domains, repeats, etc.).


.. _show-cohort-variants:

Summarize all variant alleles
-----------------------------

We can prepare a table of all variant alleles that occur in the cohort.

Each table row corresponds to a single allele and lists the variant key,
the predicted effect on the transcript (*cDNA*) and protein of interest,
the variant effects, and the number of patients who present
with one or more variant alleles (*Count*):

>>> from gpsea.view import CohortVariantViewer
>>> viewer = CohortVariantViewer(tx_id=tx_id)
>>> report = viewer.process(cohort=cohort)
>>> report  # doctest: +SKIP

.. raw:: html
  :file: report/tbx5_all_variants.html

.. doctest:: tutorial
  :hide:

  >>> if _overwrite: report.write('docs/report/tbx5_all_variants.html')


Partition the cohort by genotype and phenotype
==============================================

Testing for a genotype-phenotype association uses genotype and phenotype as variables.
In GPSEA, the variable value for an individual is computed
either by a :class:`~gpsea.analysis.clf.Classifier`
or by a :class:`~gpsea.analysis.pscore.PhenotypeScorer`.
A `Classifier` assigns the individual into a class,
whereas a `PhenotypeScorer` computes a continuous score.
The classifiers and scorers are applied on all individuals of the cohort
and the resulting variable distributions are then assessed by a statistical test.

In GPSEA, genotype is always treated as a class
and a genotype `Classifier` is a prerequisite for each analysis.
However, there is much more flexibility on the phenotype part,
where either a `Classifier` or a `PhenotypeScorer` can be used to compute the values,
depending on the analysis goals.

In this tutorial section, we first configure a `Classifier` for assigning
the individuals into a genotype class,
and we follow with generating classifiers for testing the presence or exclusion
of HPO terms in the individuals.


Partition by genotype
---------------------

In context of the tutorial, we assign each cohort member into a class
depending on presence of a single allele of a missense or truncating variant
(e.g. frameshift, stop gain, or splice site region):

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate import variant_effect, anyof
>>> from gpsea.analysis.clf import monoallelic_classifier
>>> is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id)
>>> truncating_effects = (
...    VariantEffect.TRANSCRIPT_ABLATION,
...    VariantEffect.TRANSCRIPT_TRANSLOCATION,
...    VariantEffect.FRAMESHIFT_VARIANT,
...    VariantEffect.START_LOST,
...    VariantEffect.STOP_GAINED,
...    VariantEffect.SPLICE_DONOR_VARIANT,
...    VariantEffect.SPLICE_ACCEPTOR_VARIANT,
...    # more effects could be listed here ...
... )
>>> is_truncating = anyof(variant_effect(e, tx_id) for e in truncating_effects)
>>> gt_clf = monoallelic_classifier(
...     a_predicate=is_missense,
...     b_predicate=is_truncating,
...     a_label="Missense", b_label="Truncating",
... )
>>> gt_clf.class_labels
('Missense', 'Truncating')

This is a lot of code, and detailed explanations and examples are available in the :ref:`partitioning` section.
For now, it is enough to know that the `gt_clf` will assign the individuals
into `Missense` or `Truncating` class. The individuals with the number of missense (or truncating) variants
different than one will be omitted from the analysis.


Partition by phenotype
----------------------

We use HPO terms to assign the individuals into phenotype classes,
according to the term's presence or exclusion.
The testing leverages the :ref:`true-path-rule` of ontologies.

We now prepare the classifiers for assigning into phenotype classes:

>>> from gpsea.analysis.clf import prepare_classifiers_for_terms_of_interest
>>> pheno_clfs = prepare_classifiers_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
... )


Multiple testing correction
---------------------------

By default, GPSEA performs a test for each HPO term used to annotate at least one individual in the cohort,
and there are 369 such terms in *TBX5* cohort:

>>> len(pheno_clfs)
369

However, testing multiple hypothesis on the same dataset increases the chance of receiving false positive result.
Luckily, GPSEA simplifies the application of an appropriate multiple testing correction.

For general use, we recommend using a combination
of a *phenotype MT filter* (:class:`~gpsea.analysis.mtc_filter.PhenotypeMtcFilter`) with a *multiple testing correction*.
Phenotype MT filter chooses the HPO terms to test according to several heuristics, which
reduce the multiple testing burden and focus the analysis
on the most interesting terms (see :ref:`Independent filtering for HPO <hpo-if-filter>` for more info).
Then the multiple testing correction, such as Bonferroni or Benjamini-Hochberg,
is used to control the family-wise error rate or the false discovery rate.
See :ref:`mtc` for more information.

>>> from gpsea.analysis.pcats import configure_hpo_term_analysis
>>> analysis = configure_hpo_term_analysis(hpo)

:func:`~gpsea.analysis.pcats.configure_hpo_term_analysis` configures the analysis
that uses Independent filtering for HPO (:class:`~gpsea.analysis.mtc_filter.IfHpoFilter`) for selecting HPO terms of interest,
Fisher Exact test for computing nominal p values, and Benjamini-Hochberg for multiple testing correction.


Statistical testing
===================

Now we can perform the testing and evaluate the results.

>>> result = analysis.compare_genotype_vs_phenotypes(
...     cohort=cohort,
...     gt_clf=gt_clf,
...     pheno_clfs=pheno_clfs,
... )
>>> result.total_tests
30

We only tested 30 HPO terms. This is despite the individuals being collectively annotated with
369 direct and indirect HPO terms

>>> len(result.phenotypes)
369

We can show the reasoning behind *not* testing 339 (`369 - 30`) HPO terms
by exploring the phenotype MTC filtering report:

>>> from gpsea.view import MtcStatsViewer
>>> mtc_viewer = MtcStatsViewer()
>>> mtc_report = mtc_viewer.process(result)
>>> mtc_report  # doctest: +SKIP

.. raw:: html
  :file: report/tbx5_truncating_vs_missense.mtc_report.html

.. doctest:: tutorial
  :hide:

  >>> if _overwrite: mtc_report.write('docs/report/tbx5_truncating_vs_missense.mtc_report.html')


and these are the tested HPO terms ordered by the p value corrected with the Benjamini-Hochberg procedure:

>>> from gpsea.view import summarize_hpo_analysis
>>> summary_df = summarize_hpo_analysis(hpo, result)
>>> summary_df  # doctest: +SKIP

.. csv-table:: *TBX5* truncating vs. missense
   :file: report/tbx5_truncating_vs_missense.csv
   :header-rows: 1

.. doctest:: tutorial
  :hide:

  >>> if _overwrite: summary_df.to_csv('docs/report/tbx5_truncating_vs_missense.csv')

We see that several HPO terms are significantly associated
with presence of a truncating variant in *TBX5*.
For example, `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_
was observed in 31/60 (52%) patients with a missense variant
but it was observed in 29/29 (100%) patients with a truncating variant.
Fisher exact test computed a p value of 5.61e\ :sup:`-7`
and the p value corrected by Benjamini-Hochberg procedure
is 9.55e\ :sup:`-6`.


**********
Conclusion
**********

We showed the high-level structure of genotype-phenotype association analysis using GPSEA
and we found an association between truncating *TBX5* variants
and `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_.

This is just one of many analysis types that are possible with GPSEA.
Please refer to the user guide (next section) to learn more.
