.. _custom-components:

#################
Custom components
#################

GPSEA aims to stay useful in the long run. Therefore, we took a great care
to adhere to "good software development practices" and we framed GPSEA functionalities
as a set of loosely coupled components. As a rule of thumb, each component corresponds
to a Python abstract base class which is then extended by the builtin components
and can also be extended by the future components, to serve both common or exotic use cases.

The abstract base classes define the component API.
Per guidelines in Python's :mod:`abc` module, the abstract classes use :class:`abc.ABCMeta` as a metaclass
and the class API consists of methods annotated with the :func:`abc.abstractmethod` decorator.
These decorated methods must be overridden in the subclasses.

The following sections provide guidance for customizing the most commonly used GPSEA components.


.. _custom-phenotype-scorer:

****************
Phenotype scorer
****************

:class:`~gpsea.analysis.pscore.PhenotypeScorer` computes a phenotype score for an individual.
The phenotype score is a `float` with range :math:`(-\infty, \infty)` where `NaN` indicates
that a score cannot be computed (e.g. the lab measurement value was not ascertained for the individual).

Here we show an example of a toy phenotype scorer
for using body mass index (BMI) as a phenotype score.

>>> import typing
>>> from gpsea.model import Patient
>>> from gpsea.analysis.pscore import PhenotypeScorer
>>> class BmiScorer(PhenotypeScorer):  # ❶
...    
...     def __init__(
...         self,
...         id2bmi: typing.Mapping[str, float], # ❷
...     ):
...         self._id2bmi = id2bmi
...     
...     @property
...     def name(self) -> str: # ❸
...         return "BMI phenotype scorer"
...    
...     @property
...     def description(self) -> str: # ❹
...         return "Body mass index used as a phenotype score"
...    
...     @property
...     def variable_name(self) -> str: # ❺
...         return "BMI"
...     
...     def score(self, patient: Patient) -> float: # ❻
...         try:
...             return self._id2bmi[patient.labels.label]
...         except KeyError:
...             return float('nan')

The ``BmiScorer`` must extend :class:`~gpsea.analysis.pscore.PhenotypeScorer`
to be used as a phenotype scorer (❶).
The scorer needs a mapping (e.g. a Python ``dict``) with `label` → `BMI` for the analyzed individuals (❷).
We assume the user will pre-compute the BMI values.

Then, the scorer must expose several properties, including ``name``, ``description``,
and the ``variable_name`` it operates on (❸❹❺).
GPSEA uses the properties to describe the scorer in reports or visualizations.
We should always aim for short and concise descriptions.

The most important part of the scorer is the `score` method (❻).
As stated above, the scorer is expected to compute a numerical value or `NaN`
if the individual should be excluded from the analysis.
In the case of BMI scorer, the BMI is retrieved from the ``id2bmi`` dictionary.
If the BMI is missing, `NaN` is returned and the individual is omitted from the analysis.


.. _custom-variant-predicate:

*****************
Variant predicate
*****************

A :class:`~gpsea.analysis.predicate.VariantPredicate` tests
if a variant meets a certain criterion (e.g. variant is a deletion, variant is annotated wrt. a transcript of interest)
in order to assign the individual harboring the variant into a genotype class.
GPSEA ships with an array of builtin predicates (see :mod:`gpsea.analysis.predicate` module)
that should cover the most commonly needed cases.

However, since it is unlikely that the builtin predicates cover *all* cases,
GPSEA allows to define custom variant predicates. Here we show how to create one.

As an example, we show how to create a predicate for checking if the variant affects a glycine residue
in a transcript of interest.

>>> from gpsea.model import Variant, VariantEffect
>>> from gpsea.analysis.predicate import VariantPredicate
>>> class AffectsGlycinePredicate(VariantPredicate): # ❶
...     def __init__( # ❷
...         self,
...         tx_id: str,
...     ):
...         self._tx_id = tx_id
...         self._aa_code = "Gly"
... 
...     @property
...     def name(self) -> str: # ❸
...         return "Affects Glycine"
...    
...     @property
...     def description(self) -> str: # ❹
...         return "affects a glycine residue"
...    
...     @property
...     def variable_name(self) -> str: # ❺
...         return "affected aminoacid residue"
...     
...     def test(self, variant: Variant) -> bool: # ❻
...         tx_ann = variant.get_tx_anno_by_tx_id(self._tx_id)
...         if tx_ann is not None:
...            hgvsp = tx_ann.hgvsp
...            if hgvsp is not None:
...                return hgvsp.startswith(f"p.{self._aa_code}")
...         return False
...     
...     def __eq__(self, value: object) -> bool: # ➐
...         return isinstance(value, AffectsGlycinePredicate) and self._tx_id == value._tx_id
...     
...     def __hash__(self) -> int: # ❽
...         return hash((self._tx_id,))
...     
...     def __repr__(self) -> str: # ❾
...         return str(self)
...     
...     def __str__(self) -> str: # ➓
...         return f"AffectsGlycinePredicate(tx_id={self._tx_id})"


The ``AffectsGlycinePredicate`` must extend :class:`~gpsea.analysis.predicate.VariantPredicate` to work with GPSEA (❶).
We ask the user to provide the transcript accession `str` and we set the target aminoacid code to glycine ``Gly`` (❷).

.. note::

    Clearly, to test for change of *any* aminoacid
    with only a slight rewrite of the predicate's constructor.
    We will leave this as an exercise for the interested readers.

Like in the :ref:`custom-phenotype-scorer` above, we provide metadata required for reports and visualizations (❸❹❺).

The ``test`` method includes the most important logic of the predicate (❻).
In this specific case, we retrieve the :class:`~gpsea.model.TranscriptAnnotation`
with the functional annotation data for the transcript of interest,
and we test if the HGVS protein indicates that the reference aminoacid is glycine.

.. note::

    We recommend using an Integrated Development Environment (IDE) such as PyCharm or VS Code to design the predicate.
    On top of autocompletion and syntax checking features, an IDE simplifies accessing the properties and methods of objects.
    In case of :class:`~gpsea.model.Variant`, an IDE will help us discover its ``get_tx_anno_by_tx_id`` method,
    realize that it returns either :class:`~gpsea.model.TranscriptAnnotation` or ``None``,
    and retrieve the functional annotation of the variant with respect to transcript's protein sequence
    from the ``hgvsp`` field.

Last, we override ``__eq__()`` and ``__hash__()`` (required, ➐❽) as well as ``__repr__()`` and ``__str__()`` (recommended, ❾➓).
