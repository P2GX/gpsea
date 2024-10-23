import pytest

from gpsea.analysis.predicate.genotype import VariantPredicates
from gpsea.model import *
from gpsea.model.genome import *


class TestVariantPredicates:

    def test_always_true_predicate(
        self,
        suox_cohort: Cohort,
    ):
        predicate = VariantPredicates.true()
        assert all(predicate.test(v) for v in suox_cohort.all_variants())

    @pytest.mark.parametrize(
        'effect, expected',
        [
            (VariantEffect.MISSENSE_VARIANT, True),
            (VariantEffect.SPLICE_DONOR_VARIANT, True),
            (VariantEffect.STOP_GAINED, False),
            (VariantEffect.STOP_LOST, False),
        ]
    )
    def test_variant_effect_predicate(
            self,
            missense_variant: Variant,
            effect: VariantEffect,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_effect(effect, tx_id='tx:xyz')

        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'variant_key, expected',
        [
            ('22_101_101_C_G', True),
            ('22_101_101_C_T', False),
        ]
    )
    def test_variant_key_predicate(
            self,
            missense_variant: Variant,
            variant_key: str,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_key(variant_key)

        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'exon, expected',
        [
            (1, False),
            (4, True),
            (5, False),
        ]
    )
    def test_exon_predicate(
            self,
            missense_variant: Variant,
            exon: int,
            expected: bool,
    ):
        predicate = VariantPredicates.exon(exon, tx_id='tx:xyz')

        assert predicate.test(missense_variant) == expected

    def test_exon_predicate_fails_on_invalid_exon(self):
        with pytest.raises(AssertionError) as e:
            VariantPredicates.exon(0, tx_id='tx:xyz')
        assert e.value.args[0] == '`exon` must be a positive `int`'

    @pytest.mark.parametrize(
        'tx_id, expected',
        [
            ('tx:xyz', True),
            ('other', False),
        ]
    )
    def test_transcript_predicate(
            self,
            missense_variant: Variant,
            tx_id: str,
            expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id)

        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'symbol, expected',
        [
            ('a_gene', True),
            ('b_gene', False),
        ]
    )
    def test_gene_predicate(
            self,
            missense_variant: Variant,
            symbol: str,
            expected: bool,
    ):
        predicate = VariantPredicates.gene(symbol)

        assert predicate.test(missense_variant) == expected

    def test_is_large_imprecise_sv(
        self,
        missense_variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_large_imprecise_sv()

        assert predicate.test(missense_variant) is False
        assert predicate.test(structural_variant) is True

    def test_is_structural_predicate(
        self,
        missense_variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_structural_variant()

        assert predicate.test(missense_variant) is False
        assert predicate.test(structural_variant) is True

    def test_structural_type(
        self,
        missense_variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.structural_type('SO:1000029')

        assert predicate.test(missense_variant) is False
        assert predicate.test(structural_variant) is True

    def test_variant_class(
        self,
        missense_variant: Variant,
    ):
        predicate = VariantPredicates.variant_class(VariantClass.SNV)

        assert predicate.test(missense_variant) is True

    def test_change_length(
        self,
        missense_variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.change_length('==', 0)
        
        # variant is an SNP
        assert predicate.test(missense_variant) is True

        # structural_variant is an imprecise DEL, hence False
        assert predicate.test(structural_variant) is False

    def test_change_length_is_false_for_imprecise_SVs_no_matter_what(
        self,
        structural_variant: Variant,
    ):
        for op in ("<", "<=", "==", "!=", ">=", ">"):
            predicate = VariantPredicates.change_length(op, 0)
            assert predicate.test(structural_variant) is False

    def test_structural_deletion(
        self,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_structural_deletion()

        assert predicate.test(structural_variant) is True


class TestProteinPredicates:

    @pytest.fixture(scope="class")
    def protein_metadata(self) -> ProteinMetadata:
        return ProteinMetadata(
            protein_id="pt:xyz",
            label="xyz_label",
            protein_features=(
                ProteinFeature.create(
                    FeatureInfo(name="MOCK_REPEAT", region=Region(55, 80)),
                    FeatureType.REPEAT,
                ),
                ProteinFeature.create(
                    FeatureInfo(name="MOCK_DOMAIN", region=Region(30, 50)),
                    FeatureType.DOMAIN,
                ),
            ),
            protein_length=100,
        )

    @pytest.mark.parametrize(
        'feature_type, expected',
        [
            (FeatureType.DOMAIN, True),
            (FeatureType.REPEAT, False),
        ]
    )
    def test_protein_feature_type(
            self,
            missense_variant: Variant,
            feature_type: FeatureType,
            protein_metadata: ProteinMetadata,
            expected: bool,
    ):
        predicate = VariantPredicates.protein_feature_type(
            feature_type=feature_type,
            protein_metadata=protein_metadata,
        )

        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'feature_id, expected',
        [
            ('MOCK_DOMAIN', True),
            ('REAL_DOMAIN', False),
            ('MOCK_REPEAT', False),
        ]
    )
    def test_protein_feature_id(
            self,
            missense_variant: Variant,
            feature_id: str,
            protein_metadata: ProteinMetadata,
            expected: bool,
    ):
        predicate = VariantPredicates.protein_feature(
            feature_id=feature_id,
            protein_metadata=protein_metadata,
        )

        assert predicate.test(missense_variant) == expected


class TestLogicalVariantPredicate:
    """
    Test that the AND and OR variant predicate combinators work as expected.
    """

    def test_equivalent_predicates_are_not_chained(self):
        a1 = VariantPredicates.gene(symbol='A')
        a2 = VariantPredicates.gene(symbol='A')

        assert a1 & a2 is a1
        assert a1 | a2 is a1

        assert a2 & a1 is a2
        assert a2 | a1 is a2

    @pytest.mark.parametrize(
        'left,right,expected',
        [
            ('tx:abc', 'tx:xyz', True),
            ('tx:abc', 'whatever', False),
            ('whatever', 'tx:xyz', False),
            ('whatever', 'whoever', False),
        ]
    )
    def test_und_predicate(
        self,
        missense_variant: Variant,
        left: str,
        right: str,
        expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id=left) & VariantPredicates.transcript(tx_id=right)

        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'left,right,expected',
        [
            ('tx:abc', 'tx:xyz', True),
            ('tx:abc', 'whatever', True),
            ('whatever', 'tx:xyz', True),
            ('whatever', 'whoever', False),
        ]
    )
    def test_or_predicate(
        self,
        missense_variant: Variant,
        left: str,
        right: str,
        expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id=left) | VariantPredicates.transcript(tx_id=right)
        
        assert predicate.test(missense_variant) == expected

    @pytest.mark.parametrize(
        'tx_id,expected',
        [
            ('tx:abc', False),
            ('whatever', True),
        ]
    )
    def test_inv_predicate(
        self,
        missense_variant: Variant,
        tx_id: str,
        expected: bool,
    ):
        predicate = ~VariantPredicates.transcript(tx_id)

        assert predicate.test(missense_variant) == expected
    
    def test_no_double_inv_happens(
        self,
    ):
        predicate = VariantPredicates.gene('FBN1')
        
        # Inverting a predicate must produce a new predicate.
        inv_predicate = ~predicate
        assert inv_predicate is not predicate

        # No double negation!
        inv_inv_predicate = ~inv_predicate
        assert inv_inv_predicate is predicate

    def test_empty_all_predicate_raises_error(
        self,
    ):
        with pytest.raises(ValueError) as e:
            empty = ()
            VariantPredicates.all(empty)
        assert e.value.args[0] == 'Predicates must not be empty!'
    
    def test_empty_any_predicate_raises_error(
        self,
    ):
        with pytest.raises(ValueError) as e:
            empty = ()
            VariantPredicates.any(empty)
        assert e.value.args[0] == 'Predicates must not be empty!'

    def test_logical_predicates_are_hashable(self):
        a = VariantPredicates.gene(symbol='A')
        b = VariantPredicates.gene(symbol='B')

        a_and_b = a & b
        assert isinstance(hash(a_and_b), int)

        a_or_b = a | b
        assert isinstance(hash(a_or_b), int)

        inv_a = ~a
        assert isinstance(hash(inv_a), int)
