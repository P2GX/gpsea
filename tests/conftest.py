import os

import hpotk
import pytest

from ._protein_test_service import ProteinTestMetadataService
from genophenocorr.model import *
from genophenocorr.model.genome import GRCh38, GenomicRegion, Region, Strand


def pytest_addoption(parser):
    parser.addoption(
        "--runonline", action="store_true", default=False, help="run online tests"
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "online: mark test that require internet access to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runonline"):
        # --runonline given in cli: do not skip online tests
        return
    skip_online = pytest.mark.skip(reason="need --runonline option to run")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


@pytest.fixture(scope='session')
def toy_hpo() -> hpotk.MinimalOntology:
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data', 'hp.toy.json')
    return hpotk.load_minimal_ontology(path)


@pytest.fixture(scope='session')
def toy_validation_runner(toy_hpo: hpotk.MinimalOntology) -> hpotk.validate.ValidationRunner:
    validators = (
        hpotk.validate.ObsoleteTermIdsValidator(toy_hpo),
        hpotk.validate.AnnotationPropagationValidator(toy_hpo),
        hpotk.validate.PhenotypicAbnormalityValidator(toy_hpo)
    )
    return hpotk.validate.ValidationRunner(validators)

def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    return GenomicRegion(GRCh38.contig_by_name(contig), start, end, Strand.POSITIVE)

@pytest.fixture(scope='session')
def protein_test_service() -> ProteinTestMetadataService:
    return ProteinTestMetadataService()

@pytest.fixture
def toy_cohort() -> Cohort:

    prot_id = 'NP_037407.4'

    phenos = get_test_phenotypes()

    dup = Variant(VariantCoordinates(make_region("16", 89279849, 89279850), ref='G', alt='GC', change_length=1),
                  [
                      TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                                           prot_id, Region(2230, 2231))
                  ],
                  Genotypes.from_mapping({SampleLabels('HetSingleVar'): Genotype.HETEROZYGOUS}))
    indel = Variant(VariantCoordinates(make_region("16", 89284600, 89284602), ref='GG', alt='A', change_length=-1),
                    [
                        TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', False, [VariantEffect.FRAMESHIFT_VARIANT],
                                             [9], prot_id, Region(646, 647))
                    ],
                    Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_stop_gain = Variant(VariantCoordinates(make_region("16", 89280751, 89280752), ref='G', alt='T', change_length=0),
                            [
                                TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', False, [VariantEffect.STOP_GAINED], [9], prot_id,
                             Region(1929, 1930))],
                            Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_missense = Variant(VariantCoordinates(make_region("16", 89275127, 89275128), ref='G', alt='A', change_length=0),
                           [
                               TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', False, [VariantEffect.MISSENSE_VARIANT], [10],
                             prot_id, Region(2511, 2512))
                           ],
                           Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_frameshift = Variant(VariantCoordinates(make_region("16", 89279707, 89279725), ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17),
                             [
                                 TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', False, [VariantEffect.FRAMESHIFT_VARIANT],
                              [9], prot_id, Region(2272, 2278))
                             ],
                             Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_small = Variant(VariantCoordinates(make_region("16", 89279457, 89279459), ref='TG', alt='T', change_length=-1),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                             prot_id, Region(2360, 2362))
                        ],
                        Genotypes.from_mapping({SampleLabels('HomoVar'): Genotype.HOMOZYGOUS_ALTERNATE}))
    del_large = Variant(VariantCoordinates(make_region("16", 89_190_070, 89_439_815), ref='N', alt='<DEL>', change_length=-249_745),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', None, False,
                                 [VariantEffect.STOP_LOST, VariantEffect.FEATURE_TRUNCATION, VariantEffect.CODING_SEQUENCE_VARIANT, VariantEffect.FIVE_PRIME_UTR_VARIANT,
                                  VariantEffect.THREE_PRIME_UTR_VARIANT, VariantEffect.INTRON_VARIANT], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                 prot_id, None)
                        ],
                        Genotypes.from_mapping({SampleLabels('LargeCNV'): Genotype.HETEROZYGOUS}))

    patients = (
        Patient(SampleLabels('HetSingleVar'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_F'], phenos['focal_clonic_seizure_T']),
                variants=(dup,)
                ),
        Patient(SampleLabels('HetDoubleVar1'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['seizure_T'], phenos['spasticity_T']),
                variants=(indel, snv_stop_gain)
                ),
        Patient(SampleLabels('HetDoubleVar2'),
                phenotypes=(phenos['arachnodactyly_F'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=(snv_missense, del_frameshift)
                ),
        Patient(SampleLabels('HomoVar'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=(del_small,)
                ),
        Patient(SampleLabels('LargeCNV'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_F']),
                variants=(del_large,)
                ),
    )

    return Cohort.from_patients(patients)


def get_test_phenotypes():
    phenotypes = {}

    phenotypes['arachnodactyly_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", True)
    phenotypes['seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", True)
    phenotypes['focal_clonic_seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", True)
    phenotypes['spasticity_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", True)
    phenotypes['arachnodactyly_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", False)
    phenotypes['seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", False)
    phenotypes['spasticity_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", False)
    phenotypes['focal_clonic_seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", False)

    return phenotypes

