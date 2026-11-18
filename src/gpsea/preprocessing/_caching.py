import logging
import typing

from gpsea.model import (
    ProteinMetadata,
    TranscriptAnnotation,
    TranscriptCoordinates,
    TranscriptInfoAware,
    VariantCoordinates,
)
from gpsea.model.exceptional import VariantCoordinateFinderExceptionData

from ._api import (
    FunctionalAnnotator,
    GeneCoordinateService,
    ProteinMetadataService,
    TranscriptCoordinateService,
    VariantCoordinateFinder,
    VariantCoordinateFinderException,
)
from ._cache import Cache

T = typing.TypeVar("T")

FallbackPolicy = typing.Literal["retry", "proceed"]
"""
How to handle an entry that could not be retrieved from the fallback source
on the first attempt due to issues unrelated to fallback timeouts (e.g. validation issues).

* `retry` stipulates reattempting to retrieve the value.
* `proceed` will *not* reobtain the value and report the results of the previous attempt.
"""


class CachingVariantCoordinateFinder(
    typing.Generic[T],
    VariantCoordinateFinder[T],
):
    def __init__(
        self,
        cache: Cache[VariantCoordinates | VariantCoordinateFinderExceptionData],
        fallback: VariantCoordinateFinder[T],
        key_func: typing.Callable[[T], str],
        fallback_policy: FallbackPolicy,
    ):
        self._logger = logging.getLogger(__name__ + "." + CachingVariantCoordinateFinder.__name__)
        self._cache = cache
        self._fallback = fallback
        self._key_func = key_func
        self._fallback_policy: FallbackPolicy = fallback_policy

    def find_coordinates(self, item: T) -> typing.Optional[VariantCoordinates]:
        key = self._key_func(item)
        val = self._cache.load_item(key)

        if val is None:  # cache miss
            return self._fetch_and_cache(key, item)
        else:
            if isinstance(val, VariantCoordinateFinderExceptionData):
                match self._fallback_policy:
                    case "proceed":
                        raise VariantCoordinateFinderException(
                            kind=val.kind,
                            context=val.context,
                            issues=val.issues,
                        )
                    case "retry":
                        return self._fetch_and_cache(key, item)
                    case _:
                        raise ValueError(f"Invalid failure policy {self._on_failure}")
            else:
                return val

    def _fetch_and_cache(self, key: str, item: T) -> typing.Optional[VariantCoordinates]:
        self._logger.debug("Fetching: %s", key)
        try:
            val = self._fallback.find_coordinates(item)
        except VariantCoordinateFinderException as e:
            if e.kind != "timeout":
                data = VariantCoordinateFinderExceptionData(
                    kind=e.kind,
                    context=e.context,
                    issues=e.issues,
                )
                self._cache.store_item(key, data)
            raise e
        if val is None:
            self._logger.debug("Unable to find the coordinates: %s", key)
        else:
            self._logger.debug("Received the coordinates: %s", key)
            self._cache.store_item(key, val)
        return val


class CachingProteinMetadataService(ProteinMetadataService):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        cache: Cache[ProteinMetadata],
        fallback: ProteinMetadataService,
    ):
        self._logger = logging.getLogger(__name__ + "." + CachingProteinMetadataService.__name__)
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, ProteinMetadataService)
        self._fallback = fallback

    def annotate(self, protein_id: str) -> ProteinMetadata:
        assert isinstance(protein_id, str)
        item = self._cache.load_item(protein_id)
        if item is None:  # cache miss
            self._logger.debug("Fetching: %s", protein_id)
            item = self._fallback.annotate(protein_id)
            self._cache.store_item(protein_id, item)
        else:
            self._logger.debug("Cached: %s", protein_id)

        return item


class CachingTranscriptCoordinateService(TranscriptCoordinateService):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        cache: Cache[TranscriptCoordinates],
        fallback: TranscriptCoordinateService,
    ):
        self._logger = logging.getLogger(__name__ + "." + CachingTranscriptCoordinateService.__name__)
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, TranscriptCoordinateService)
        self._fallback = fallback

    def fetch(
        self,
        tx: typing.Union[str, TranscriptInfoAware],
    ) -> TranscriptCoordinates:
        tx_id = self._parse_tx(tx)
        item = self._cache.load_item(tx_id)
        if item is None:  # cache miss
            self._logger.debug("Fetching: %s", tx_id)
            item = self._fallback.fetch(tx_id)
            self._cache.store_item(tx_id, item)
        else:
            self._logger.debug("Cached: %s", tx_id)

        return item


class CachingFunctionalAnnotator(FunctionalAnnotator):
    """A class that retrieves a Variant object if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(variant_coordinates:VariantCoordinates): Gets data and returns a Variant object for given variant coordinates
    """

    def __init__(
        self,
        cache: Cache[typing.Sequence[TranscriptAnnotation]],
        fallback: FunctionalAnnotator,
    ):
        self._logger = logging.getLogger(__name__ + "." + CachingFunctionalAnnotator.__name__)
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, FunctionalAnnotator)
        self._fallback = fallback

    @staticmethod
    def _create_cache_key(vc: VariantCoordinates) -> str:
        vk = vc.variant_key
        if len(vk) <= 50:
            return f"{vk}"
        else:
            # long INDELs in sequence notation
            return f"{vc.chrom}_{vc.start}_{vc.end}_{vc.variant_class}"

    def annotate(
        self,
        variant_coordinates: VariantCoordinates,
    ) -> typing.Sequence[TranscriptAnnotation]:
        assert isinstance(variant_coordinates, VariantCoordinates)

        cache_key = CachingFunctionalAnnotator._create_cache_key(variant_coordinates)
        annotations = self._cache.load_item(cache_key)
        if annotations is None:  # cache miss
            self._logger.debug("Fetching: %s", cache_key)
            annotations = self._fallback.annotate(variant_coordinates)
            self._cache.store_item(cache_key, annotations)
        else:
            self._logger.debug("Cached: %s", cache_key)

        return annotations


class CachingGeneCoordinateService(GeneCoordinateService):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        cache: Cache[typing.Sequence[TranscriptCoordinates]],
        fallback: GeneCoordinateService,
    ):
        self._logger = logging.getLogger(__name__ + "." + CachingGeneCoordinateService.__name__)
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, GeneCoordinateService)
        self._fallback = fallback

    def fetch_for_gene(
        self,
        gene: str,
    ) -> typing.Sequence[TranscriptCoordinates]:
        assert isinstance(gene, str)

        # `gene` may be treated as a file name by the cache.
        item = self._cache.load_item(gene)
        if item is None:  # cache miss
            self._logger.debug("Fetching: %s", gene)
            item = self._fallback.fetch_for_gene(gene)
            self._cache.store_item(gene, item)
        else:
            self._logger.debug("Cached: %s", gene)

        return item
