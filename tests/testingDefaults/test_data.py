import os
import hpotk

from genophenocorr.variant import Variant, VariantCoordinates, TranscriptAnnotation
from genophenocorr.phenotype import Phenotype
from genophenocorr.protein import ProteinMetadata, SimpleProteinFeature, FeatureType, FeatureInfo
from genophenocorr.patient import Patient
from genophenocorr.cohort import Cohort


def get_test_cohort() -> Cohort:

    prot = ProteinMetadata(protein_id='NP_037407.4', label='Ankyrin repeat domain-containing protein 11', 
            protein_features=(SimpleProteinFeature(feature_type=FeatureType.REPEAT, info=FeatureInfo(name='ANK 1', start=167, end=196)), 
            SimpleProteinFeature(feature_type=FeatureType.REPEAT, info=FeatureInfo(name='ANK 2', start=200, end=229)), 
            SimpleProteinFeature(feature_type=FeatureType.REPEAT, info=FeatureInfo(name='ANK 3', start=233, end=262)), 
            SimpleProteinFeature(feature_type=FeatureType.REPEAT, info=FeatureInfo(name='ANK 4', start=266, end=292)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=1, end=90)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=128, end=169)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=289, end=380)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=398, end=647)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=723, end=783)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=881, end=1043)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=1059, end=1393)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=1424, end=1710)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=1814, end=1836)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=1988, end=2019)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Disordered', start=2131, end=2406)), 
            SimpleProteinFeature(feature_type=FeatureType.REGION, info=FeatureInfo(name='Important for protein degradation', start=2369, end=2663))))

    HetSingleVar = [Variant('16_89279851_-/C', 'insertion', 
                VariantCoordinates(chrom=16, start=89279849, end=89279851, ref='G', alt='GC', change_length=1, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', ['frameshift_variant'], 
                [9], [prot], 2231, 2231)],genotype='heterozygous')]
    HetDoubleVar1 = [Variant('16_89284601_GG/A', 'indel', 
                VariantCoordinates(chrom=16, start=89284600, end=89284602, ref='GG', alt='A', change_length=-1, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', ['frameshift_variant'], 
                [9], [prot], 647, 647)],genotype='heterozygous'),
                Variant('16_89280752_G/T', 'SNV', 
                VariantCoordinates(chrom=16, start=89280751, end=89280752, ref='G', alt='T', change_length=0, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', ['stop_gained'], 
                [9], [prot], 1930, 1930)],genotype='heterozygous')]
    HetDoubleVar2 = [Variant('16_89275128_G/A', 'SNV', 
                VariantCoordinates(chrom=16, start=89275127, end=89275128, ref='G', alt='A', change_length=0, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', ['missense_variant'], 
                [10], [prot], 2512, 2512)],genotype='heterozygous'),
                Variant('16_89279708_AGTGTTCGGGGCGGGGCC/A', 'indel', 
                VariantCoordinates(chrom=16, start=89279707, end=89279725, ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', ['frameshift_variant'], 
                [9], [prot], 2273, 2278)],genotype='heterozygous')]
    HomoVar = [Variant('16_89279458_TG/T', 'indel', 
                VariantCoordinates(chrom=16, start=89279457, end=89279459, ref='TG', alt='T', change_length=-1, genotype='homozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', ['frameshift_variant'], 
                [9], [prot], 2361, 2362)],genotype='homozygous')]
    LargeCNV = [Variant('16_89190071_deletion', 'deletion', 
                VariantCoordinates(chrom=16, start=89190070, end=89439815, ref='N', alt='<DEL>', change_length=4, genotype='heterozygous'), 
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', None, ['stop_lost', 'feature_truncation', 'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'intron_variant'], 
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], [prot], None, None)],genotype='heterozygous')]

    phenos = get_test_phenotypes()

    patients = (
        Patient('HetSingleVar',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_F'], phenos['focal_clonic_seizure_T']),
                variants=HetSingleVar,
                proteins=[prot]
                ),
        Patient('HetDoubleVar1',
                phenotypes=(phenos['arachnodactyly_T'], phenos['seizure_T'], phenos['spasticity_T']),
                variants=HetDoubleVar1,
                proteins=[prot]
                ),
        Patient('HetDoubleVar2',
                phenotypes=(phenos['arachnodactyly_F'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=HetDoubleVar2,
                proteins=[prot]
                ),
        Patient('HomoVar',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=HomoVar,
                proteins=[prot]
                ),
        Patient('LargeCNV',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_F']),
                variants=LargeCNV,
                proteins=[prot]
                ),
    )

    return Cohort.from_patients(patients)


def get_test_phenotypes():
    hpo = get_toy_hpo()

    phenotypes = {}

    phenotypes['arachnodactyly_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), hpo.get_term('HP:0001166').name, True)
    phenotypes['seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), hpo.get_term('HP:0001250').name, True)
    phenotypes['focal_clonic_seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), hpo.get_term('HP:0002266').name, True)
    phenotypes['spasticity_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), hpo.get_term('HP:0001257').name, True)
    phenotypes['arachnodactyly_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), hpo.get_term('HP:0001166').name, False)
    phenotypes['seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), hpo.get_term('HP:0001250').name, False)
    phenotypes['spasticity_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), hpo.get_term('HP:0001257').name, False)
    phenotypes['focal_clonic_seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), hpo.get_term('HP:0002266').name, False)

    return phenotypes


def get_toy_hpo():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hp.toy.json')
    hpo = hpotk.ontology.load.obographs.load_ontology(str(path))
    return hpo