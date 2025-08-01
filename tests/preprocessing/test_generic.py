import typing

import hpotk
import pytest

from gpsea.model import ImpreciseSvInfo, TranscriptCoordinates, VariantClass, VariantEffect
from gpsea.preprocessing import DefaultImpreciseSvFunctionalAnnotator, GeneCoordinateService


class MockGeneCoordinateService(GeneCoordinateService):
    def __init__(
        self,
        genes: typing.Mapping[str, typing.Sequence[TranscriptCoordinates]],
    ):
        self._genes = genes

    def fetch_for_gene(
        self,
        gene: str,
    ) -> typing.Sequence[TranscriptCoordinates]:
        return self._genes[gene]


class TestDefaultImpreciseSvFunctionalAnnotator:
    @pytest.fixture(scope="class")
    def gcs(
        self,
        suox_mane_tx_coordinates: TranscriptCoordinates,
    ) -> GeneCoordinateService:
        return MockGeneCoordinateService(
            {
                "HGNC:11460": (suox_mane_tx_coordinates,),
            }
        )

    @pytest.fixture(scope="class")
    def annotator(
        self,
        gcs: GeneCoordinateService,
    ) -> DefaultImpreciseSvFunctionalAnnotator:
        return DefaultImpreciseSvFunctionalAnnotator(gcs)

    def test_annotate(
        self,
        annotator: DefaultImpreciseSvFunctionalAnnotator,
    ):
        sv_info = ImpreciseSvInfo(
            structural_type=hpotk.TermId.from_curie("SO:1000029"),
            variant_class=VariantClass.DEL,
            gene_id="HGNC:11460",
            gene_symbol="SUOX",
        )

        tas = annotator.annotate(item=sv_info)

        assert len(tas) == 1

        ta = tas[0]

        assert ta.gene_symbol == "SUOX"
        assert ta.transcript_id == "NM_001032386.2"
        assert ta.hgvs_cdna is None
        assert not ta.is_preferred
        assert ta.variant_effects == (VariantEffect.TRANSCRIPT_ABLATION,)
        assert ta.overlapping_exons is not None and ta.overlapping_exons == (0, 1, 2, 3, 4)
        assert ta.protein_id is None
        assert ta.hgvsp is None
        assert ta.protein_effect_location is None
