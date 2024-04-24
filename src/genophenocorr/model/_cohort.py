import typing
from collections import Counter

from ._base import SampleLabels
from ._phenotype import Phenotype, Disease
from ._variant import Variant


class Patient:
    """A class that represents an individual patient

    Attributes:
        patient_id (SampleLabels): The patient identifiers
        phenotypes (Sequence[Phenotype]): A list of Phenotype objects
        variants (Sequence[Variant]): A list of Variant objects
    """

    def __init__(self, labels: SampleLabels,
                 phenotypes: typing.Iterable[Phenotype],
                 diseases: typing.Iterable[Disease],
                 variants: typing.Iterable[Variant]):
        """Constructs all necessary attributes for a Patient object

        Args:
            labels (string): A string unique to this Patient object
            phenotypes (Iterable[Phenotype]): A list of Phenotype objects
            variants (Iterable[Variant]): A list of Variant objects
        """
        self._labels = labels
        self._phenotypes = tuple(phenotypes)
        self._diseases = tuple(diseases)
        self._variants = tuple(variants)

    @property
    def patient_id(self) -> str:
        """
        Returns:
            string: Patient ID unique to this Patient object
        """
        return self._labels.label_summary()

    @property
    def labels(self) -> SampleLabels:
        """
        Get the sample identifiers.
        """
        return self._labels

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        """
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects associated with this Patient object
        """
        return self._phenotypes

    @property
    def diseases(self) -> typing.Sequence[Disease]:
        return self._diseases

    @property
    def variants(self) -> typing.Sequence[Variant]:
        """
        Returns:
            Sequence[Variant]: A list of Variant objects associated with this Patient object
        """
        return self._variants

    def present_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over *present* phenotypes of the patient.
        """
        return filter(lambda p: p.is_observed, self._phenotypes)

    def excluded_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over *excluded* phenotypes of the patient.
        """
        return filter(lambda p: p.is_excluded, self._phenotypes)
    
    def present_diseases(self) -> typing.Iterator[Disease]:
        return filter(lambda d: d.is_present, self._diseases)
    
    def excluded_diseases(self) -> typing.Iterator[Disease]:
        return filter(lambda d: not d.is_present, self._diseases)

    def __str__(self) -> str:
        return (f"Patient("
                f"labels:{self._labels}, "
                f"variants:{self.variants}, "
                f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, "
                f"diseases:{[dis.identifier for dis in self.diseases]}")

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return (isinstance(other, Patient)
                and self.patient_id == other.patient_id
                and self.variants == other.variants
                and self.phenotypes == other.phenotypes
                and self.diseases == other.diseases)

    def __hash__(self) -> int:
        return hash((self.patient_id, self.variants, self.phenotypes, self.diseases)) 


class Cohort(typing.Sized):

    @staticmethod
    def from_patients(members: typing.Sequence[Patient], include_patients_with_no_HPO: bool = False, include_patients_with_no_variants: bool = False):
        """
        Create a cohort from a sequence of patients.

        :param members: a sequence of cohort members.
        :return: the cohort
        """
        cohort_variants, cohort_phenotypes, cohort_diseases = set(), set(), set()
        var_counts, pheno_count, diseases_count = Counter(), Counter(), Counter()
        members = set(members)
        excluded_members = []
        for patient in members:
            if len(patient.phenotypes) == 0 and not include_patients_with_no_HPO:
                excluded_members.append(patient)
                continue
            if len(patient.variants) == 0 and not include_patients_with_no_variants:
                excluded_members.append(patient)
                continue
            cohort_phenotypes.update(patient.phenotypes)
            cohort_variants.update(patient.variants)
            cohort_diseases.update(patient.diseases)
            var_counts.update([var.variant_coordinates.variant_key for var in patient.variants])
            pheno_count.update(pheno.identifier.value for pheno in patient.present_phenotypes())
            diseases_count.update(dis.identifier.value for dis in patient.present_diseases())
        all_counts = {'patients': len(members), 'variants': var_counts, 'phenotypes': pheno_count, 'diseases': diseases_count} 
        return Cohort(members, cohort_phenotypes, cohort_variants, cohort_diseases,
                      all_counts, excluded_members) 

    """This class creates a collection of patients and makes it easier to determine overlapping diseases,
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping.

    Attributes:
        all_patients (Sequence[Patient]): A set of all Patient objects in the Cohort
        all_phenotypes (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
        all_variants (Sequence[Variant]): A set of all Variant objects in the Cohort
        all_transcripts (Sequence[string]): A set of all transcript IDs referenced in all the Variant objects
        total_patient_count (integer): The total number of Patient objects
    Methods:
        list_all_patients(): A list of all patient IDs
        list_all_phenotypes(top:Optional[integer]): A list of all the top phenotype IDs (or all IDs if top is None) and the count of how many patients have it.
        list_all_variants(top:Optional[integer]): A list of all the top variants (or all variants if top is None) and the count of how many patients have it.
        list_all_proteins(top:Optional[integer]): A list of all the top protein IDs (or all IDs if top is None) and the count of how many patients have it.
        list_data_by_tx(transcript:Optional[string]): A list and count of all the variants effects found for all transcripts or a given transcript if transcript is not None.
    """

    def __init__(self, patient_set: typing.Set[Patient], phenotype_set, variant_set, disease_set, counts_dict, excluded_members,
                 recessive=False):
        """Constructs all necessary attributes for a Cohort object

        Args:
            patient_set (Sequence[Patient]): A set of all Patient objects in the Cohort
            phenotype_set (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
            variant_set (Sequence[Variant]): A set of all Variant objects in the Cohort
            protein_set (Sequence[ProteinMetadata]): A set of all ProteinMetadata objects in the Cohort
            counts_dict (Dictionary{String, Counter()}): A Dictionary with counts for Phenotypes, Variant, and Proteins objects represented in all Patients
            recessive (boolean, Optional): True if the Cohort is focused on a recessive allele. Defaults to False.
        """
        if not isinstance(patient_set, set):
            raise ValueError(f'`patient_set` must be a set but got {type(patient_set)}')
        else:
            self._patient_set = frozenset(patient_set)

        self._phenotype_set = phenotype_set
        self._variant_set = variant_set
        self._disease_set = disease_set
        self._all_counts_dict = counts_dict
        self._excluded_members = excluded_members
        self._recessive = recessive

    @property
    def all_patients(self) -> typing.FrozenSet[Patient]:
        """
        Returns:
            set: A frozen set of all the Patient objects in the Cohort
        """
        return self._patient_set

    @property
    def all_phenotypes(self):
        """
        Returns:
            set: A set of all the Phenotype objects in the Cohort
        """
        return self._phenotype_set
    
    @property
    def all_diseases(self):
        return self._disease_set

    @property
    def all_variants(self) -> typing.Collection[Variant]:
        """
        Returns:
            set: A set of all the Variant objects in the Cohort
        """
        return self._variant_set

    @property
    def all_transcripts(self):
        """
        Returns:
            set: A set of all the transcript IDs in the Cohort
        """
        all_trans = set()
        for var in self.all_variants:
            all_trans.update([trans.transcript_id for trans in var.tx_annotations])
        return all_trans

    @property
    def all_excluded_patients(self):
        return self._excluded_members

    @property
    def total_patient_count(self):
        """
        Returns:
            integer: The count of all the Patient objects
        """
        return self._all_counts_dict.get('patients')

    def list_all_patients(self):
        """
        Returns:
            list: A list of all the patient IDs in the Cohort
        """
        return [pat.patient_id for pat in self.all_patients]

    def list_all_phenotypes(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all phenotypes. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (phenotype ID, number of patients with that phenotype)
        """
        return self._all_counts_dict.get('phenotypes').most_common(top)
    
    def list_all_diseases(self, top=None):
        return self._all_counts_dict.get('diseases').most_common(top)

    def list_all_variants(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all variants. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (variant string, number of patients with that variant)
        """
        return self._all_counts_dict.get('variants').most_common(top)
    
    def list_all_proteins(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all proteins. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (protein ID string, number of variants with that protein)
        """
        prots = Counter()
        for pat in [pats for pats in self.all_patients if pats not in self.all_excluded_patients]:
            for var in pat.variants:
                for trans in var.tx_annotations:
                    prots.update([trans.protein_id])
        return prots.most_common(top)

    def list_data_by_tx(self, transcript=None):
        """
        Args:
            transcript (string, Optional): If not given, lists all transcripts. Otherwise, will only list the given transcript
        Returns:
            dictionary: Each transcript ID references a Counter(), with the variant effect as the key and total variants with that effect as the count value
        """
        if transcript is not None:
            var_type_dict = {transcript: Counter()} 
        else:
            var_type_dict = {tx_id: Counter() for tx_id in self.all_transcripts}
        for var in self.all_variants:
            for trans in var.tx_annotations:
                if trans.transcript_id in var_type_dict:
                    var_type_dict.get(trans.transcript_id).update([var_eff.name for var_eff in trans.variant_effects])
        too_small = []
        for tx_id, var_effect_counter in var_type_dict.items():
            if len(var_effect_counter) <= 1:
                too_small.append(tx_id)
        for tx_id in too_small:
            del var_type_dict[tx_id]
        return var_type_dict

    def get_excluded_ids(self):
        return [ex.patient_id for ex in self.all_excluded_patients]

    def get_excluded_count(self):
        return len(self.all_excluded_patients)

    def __len__(self) -> int:
        return len(self._patient_set)
