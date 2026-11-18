import abc
import io
import json
import os
import pickle
import typing
import uuid

from gpsea.io import GpseaJSONDecoder, GpseaJSONEncoder

from ._api import (
    Closeable,
)

T = typing.TypeVar("T")


class Cache(typing.Generic[T], metaclass=abc.ABCMeta):
    """
    Cache stores items locally to prevent having to consult a remote resource.

    The itemes are persisted across runs.
    """

    @abc.abstractmethod
    def load_item(self, identifier: str) -> typing.Optional[T]:
        pass

    @abc.abstractmethod
    def store_item(self, identifier: str, item: T):
        pass


class FilesystemCache(typing.Generic[T], Cache[T], metaclass=abc.ABCMeta):
    def __init__(
        self,
        data_dir: str,
    ):
        if not os.path.isdir(data_dir):
            raise ValueError(f"datadir {data_dir} must be an existing directory")
        self._datadir = data_dir

    @property
    def data_dir(self) -> str:
        return self._datadir

    def load_item(self, identifier: str) -> typing.Optional[T]:
        name = self._prepare_resource_path(identifier=identifier)
        if os.path.isfile(name):
            with open(name, "rb") as fh:
                return self._deserialize(fh)  # type: ignore
        else:
            return None

    def store_item(self, identifier: str, item: T):
        name = self._prepare_resource_path(identifier=identifier)
        with open(name, "wb") as fh:
            self._serialize(item, fh)  # type: ignore

    @abc.abstractmethod
    def _prepare_resource_path(self, identifier: str) -> str:
        pass

    @abc.abstractmethod
    def _serialize(self, item: T, fh: io.BytesIO):
        pass

    @abc.abstractmethod
    def _deserialize(self, fh: io.BytesIO) -> T:
        pass


class PicklingCache(typing.Generic[T], FilesystemCache[T]):
    def _prepare_resource_path(self, identifier: str) -> str:
        return os.path.join(self._datadir, f"{identifier}.pickle")

    def _serialize(self, item: T, fh: io.BytesIO):
        pickle.dump(item, fh)

    def _deserialize(self, fh: io.BytesIO) -> T:
        return pickle.load(fh)


class JsonCache(typing.Generic[T], FilesystemCache[T]):
    def __init__(
        self,
        data_dir: str,
        indent: typing.Optional[typing.Union[int, str]] = None,
    ):
        super().__init__(data_dir=data_dir)
        self._indent = indent

    def _prepare_resource_path(self, identifier: str) -> str:
        return os.path.join(self._datadir, f"{identifier}.json")

    def _serialize(
        self,
        item: T,
        fh: io.BytesIO,
    ):
        text_fh = JsonCache._wrap_with_text(fh)
        json.dump(
            item,
            text_fh,
            cls=GpseaJSONEncoder,
            indent=self._indent,
        )

    def _deserialize(
        self,
        fh: io.BytesIO,
    ) -> T:
        text_fh = JsonCache._wrap_with_text(fh)
        return json.load(fp=text_fh, cls=GpseaJSONDecoder)

    @staticmethod
    def _wrap_with_text(fh: io.BytesIO) -> typing.TextIO:
        return io.TextIOWrapper(fh)


class FileCache(typing.Generic[T], Cache[T], Closeable):
    """
    Sometimes we need to cache items that are not easy to store in files
    named by their description/key (e.g. HGVS variant annotation) due
    to use of characters that are forbidden in the file paths on the major
    platforms.

    Therefore, the ``FileCache`` maintains a mapping
    from the key to UUID and the UUID is used as a file name.

    The class loads the content of the mapping, if present,
    and stores it upon its destruction. Note, the class is *not* a context
    manager because the current preprocessing APIs do not allow such usage.

    The ``cache_dir`` must point to an already existing directory.
    An error is raised otherwise.
    """

    _MAPPING_FILE_NAME: str = "GPSEA_FILE_CACHE_MAPPING.json"

    def __init__(
        self,
        cache_dir: str,
    ) -> None:
        assert os.path.isdir(cache_dir)
        self._cache_dir = cache_dir
        self._cache_file = os.path.join(cache_dir, FileCache._MAPPING_FILE_NAME)
        if os.path.exists(self._cache_file):
            if os.path.isfile(self._cache_file):
                with open(self._cache_file, "r") as fh:
                    self._vals: typing.MutableMapping[str, str] = json.load(fh)
            else:
                raise ValueError(f"{self._cache_file} exists yet it is not a file")
        else:
            self._vals = {}

    def load_item(
        self,
        identifier: str,
    ) -> typing.Optional[T]:
        uid = self._vals.get(identifier)
        if uid is None:
            return None  # cache miss
        else:
            fname = FileCache._make_filename(uid)
            fpath = os.path.join(self._cache_dir, fname)
            with open(fpath, "r") as fh:
                return json.load(
                    fh,
                    cls=GpseaJSONDecoder,
                )

    def store_item(
        self,
        identifier: str,
        item: T,
    ):
        if identifier in self._vals:
            # already cached
            return
        else:
            uid = str(uuid.uuid4())
            fname = FileCache._make_filename(uid)
            fpath = os.path.join(self._cache_dir, fname)
            with open(fpath, "w") as fh:
                json.dump(
                    item,
                    fh,
                    cls=GpseaJSONEncoder,
                )
            self._vals[identifier] = uid

    @staticmethod
    def _make_filename(uid: str) -> str:
        return f"{uid}.json"

    def close(self):
        with open(self._cache_file, "w") as fh:
            json.dump(self._vals, fh)
