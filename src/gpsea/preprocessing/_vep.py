# A module with classes that interact with Ensembl's REST API to fetch required data.
import logging
import typing

import requests

from gpsea.model import VariantCoordinates, TranscriptAnnotation, VariantEffect
from gpsea.model.genome import Region
from ._api import FunctionalAnnotator


class VepFunctionalAnnotator(FunctionalAnnotator):
    """
    `VepFunctionalAnnotator` uses the Variant Effect Predictor (VEP) REST API
    to perform functional variant annotation.
    """

    NONCODING_EFFECTS = {
        VariantEffect.UPSTREAM_GENE_VARIANT,
        VariantEffect.FIVE_PRIME_UTR_VARIANT,
        VariantEffect.NON_CODING_TRANSCRIPT_VARIANT,
        VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT,
        VariantEffect.SPLICE_ACCEPTOR_VARIANT,
        VariantEffect.SPLICE_DONOR_VARIANT,
        VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT,
        VariantEffect.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
        VariantEffect.INTRON_VARIANT,
        VariantEffect.THREE_PRIME_UTR_VARIANT,
        VariantEffect.DOWNSTREAM_GENE_VARIANT,
        VariantEffect.INTERGENIC_VARIANT,
    }
    """
    Non-coding variant effects where we do not complain if the functional annotation lacks the protein effects.
    """

    _REFSEQ_TX_ID_PREFIXES = (
        "NM_",
        "NR_",
        "NC_",
        "XM_",
        "XR_",
        "XC_",
    )

    def __init__(self, include_computational_txs: bool = False, timeout: float = 10.0):
        self._logger = logging.getLogger(__name__)
        self._url = (
            "https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1"
            "&domains=1&hgvs=1"
            "&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1"
            "&transcript_version=1&variant_class=1"
        )
        self._include_computational_txs = include_computational_txs
        self._timeout = timeout

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        response = self.fetch_response(variant_coordinates)
        return self.process_response(variant_coordinates.variant_key, response)

    def process_response(
        self,
        variant_key: str,
        response: typing.Mapping[str, typing.Any],
    ) -> typing.Sequence[TranscriptAnnotation]:
        annotations = []
        if "transcript_consequences" not in response:
            raise ValueError(
                "The VEP response for `%s` lacked the required `transcript_consequences` field. %s",
                variant_key,
                response,
            )
        for trans in response["transcript_consequences"]:
            annotation = self._process_item(trans)
            if annotation is not None:
                annotations.append(annotation)

        return annotations

    def _parse_variant_effect(self, effect: str) -> typing.Optional[VariantEffect]:
        effect = effect.upper()
        if effect == "5_PRIME_UTR_VARIANT":
            effect = "FIVE_PRIME_UTR_VARIANT"
        elif effect == "3_PRIME_UTR_VARIANT":
            effect = "THREE_PRIME_UTR_VARIANT"
        try:
            var_effect = VariantEffect[effect]
        except KeyError:
            # A missing variant effect, pls submit an issue to the GPSEA GitHub repository.
            raise ValueError("VariantEffect %s was not found in our record of possible effects.", effect)
        return var_effect

    def _process_item(self, item: typing.Dict) -> typing.Optional[TranscriptAnnotation]:
        """
        Parse one transcript annotation from the JSON response.
        """
        trans_id = item.get("transcript_id")
        assert isinstance(trans_id, str)
        if VepFunctionalAnnotator._seems_like_refseq_tx(trans_id):
            if not self._include_computational_txs and not trans_id.startswith("NM_"):
                # Skipping a computational transcript
                return None
        is_preferred = True if ("canonical" in item and item["canonical"] == 1) else False
        hgvs_cdna = item.get("hgvsc")
        var_effects = []
        consequences = item.get("consequence_terms")
        for con in consequences:
            var_effect = self._parse_variant_effect(con)
            if var_effect is not None:
                var_effects.append(var_effect)
        gene_name = item.get("gene_symbol")
        exons_affected = item.get("exon")
        if exons_affected is not None:
            exons_affected = exons_affected.split("/")[0].split("-")
            if len(exons_affected) == 2:
                exons_affected = range(int(exons_affected[0]), int(exons_affected[1]) + 1)
            exons_affected = (int(x) for x in exons_affected)

        protein_id = item.get("protein_id")
        hgvsp = item.get("hgvsp")
        protein_effect_start = item.get("protein_start")
        protein_effect_end = item.get("protein_end")
        if protein_effect_start is None or protein_effect_end is None:
            if not any(ve in VepFunctionalAnnotator.NONCODING_EFFECTS for ve in var_effects):
                self._logger.warning(
                    "Missing start/end coordinate for %s on protein %s. Protein effect will not be included.",
                    hgvs_cdna,
                    protein_id,
                )
            protein_effect = None
        else:
            # The coordinates are in 1-based system and we need 0-based.
            protein_effect_start = int(protein_effect_start) - 1
            protein_effect_end = int(protein_effect_end)
            protein_effect = Region(protein_effect_start, protein_effect_end)

        return TranscriptAnnotation(
            gene_name,
            trans_id,
            hgvs_cdna,
            is_preferred,
            var_effects,
            exons_affected,
            protein_id,
            hgvsp,
            protein_effect,
        )

    def fetch_response(
        self,
        variant_coordinates: VariantCoordinates,
    ) -> typing.Mapping[str, typing.Any]:
        """
        Get a `dict` with the response from the VEP REST API.
        Args:
            variant_coordinates: a query :class:`~gpsea.model.VariantCoordinates`.
        """
        api_url = self._url % (VepFunctionalAnnotator.format_coordinates_for_vep_query(variant_coordinates))
        r = requests.get(api_url, headers={"Accept": "application/json"}, timeout=self._timeout)
        # Throw an exception rather than errors so we can skip the variant in _phenopackets
        if not r.ok:
            self._logger.error("Expected a result but got an Error for variant: %s", variant_coordinates.variant_key)
            self._logger.error(r.text)
            raise ValueError("Expected a result but got an Error. See log for details.")
        results = r.json()
        if not isinstance(results, list):
            self._logger.error(results.get("error"))
            raise ValueError("Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logger.error(
                "Expected only one variant per request but received %s different variants.", len(results)
            )
            self._logger.error([result.id for result in results])
            raise ValueError(f"Expected only one variant per request but received {len(results)} different variants.")
        return results[0]

    @staticmethod
    def format_coordinates_for_vep_query(vc: VariantCoordinates) -> str:
        """
        Converts the 0-based VariantCoordinates to ones that will be interpreted
        correctly by VEP

        Example - an insertion/duplication of G after the given G at coordinate 3:
        1 2 3 4 5
        A C G T A

        0-based: 2 3 G GG       1-based: 3 G GG         VEP: 4 3 - G

        Args:
            vc (VariantCoordinates): A VariantCoordinates object
        Returns:
            str: The variant coordinates formatted to work with VEP
        """

        chrom = vc.chrom
        end = vc.end
        start = vc.start + 1
        alt = vc.alt
        if vc.is_structural():
            alt = vc.alt[1:-1]
            # TODO: Verify <INS> are working correctly
        else:
            if len(vc.ref) == 0 or len(vc.alt) == 0:
                raise ValueError("Trimmed alleles are not yet supported!")
            if len(vc.ref) == 1 and len(vc.alt) != 1:
                # INS/DUP
                start = start + 1  # we must "trim"
                alt = vc.alt[1:]
                # 100 AC AGT
                # MNV

        return f"{chrom}:{start}-{end}/{alt}"

    @staticmethod
    def _seems_like_refseq_tx(tx_id: str) -> bool:
        return len(tx_id) >= 3 and tx_id[:3] in VepFunctionalAnnotator._REFSEQ_TX_ID_PREFIXES
