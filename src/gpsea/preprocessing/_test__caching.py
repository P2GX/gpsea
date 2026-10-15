import os
import re
import uuid

from gpsea.model import Age

from ._caching import FileCache, JsonCache, PicklingCache


class TestJsonCache:
    def test_can_create_cache(self, tmp_path):
        data_dir = str(tmp_path)
        cache = JsonCache(data_dir=data_dir)

        assert cache is not None

    def test_store_and_load_item_roundtrip(
        self,
        tmp_path,
    ):
        data_dir = str(tmp_path)
        cache = JsonCache(data_dir=data_dir)

        identifier = "an_age"
        item = Age.birth()

        # We start with clean cache ...
        assert len(os.listdir(cache.data_dir)) == 0

        cache.store_item(identifier, item)

        # ... and now there is one item in the cache.
        items = tuple(os.listdir(cache.data_dir))
        assert items == ("an_age.json",)

        # Let's check we can get the item back.
        loaded = cache.load_item(identifier)
        assert loaded is not None
        assert loaded == item


class TestPicklingCache:
    def test_can_create_cache(self, tmp_path):
        data_dir = str(tmp_path)
        cache = PicklingCache(data_dir=data_dir)

        assert cache is not None

    def test_store_and_load_item_roundtrip(
        self,
        tmp_path,
    ):
        data_dir = str(tmp_path)
        cache = PicklingCache(data_dir=data_dir)

        identifier = "an_age"
        item = Age.birth()

        # We start with clean cache ...
        assert len(os.listdir(cache.data_dir)) == 0

        cache.store_item(identifier, item)

        # ... and now there is one item in the cache.
        items = tuple(os.listdir(cache.data_dir))
        assert items == ("an_age.pickle",)

        # Let's check we can get the item back.
        loaded = cache.load_item(identifier)
        assert loaded is not None
        assert loaded == item


class TestFileCache:
    def test_can_create_cache(self, tmp_path):
        data_dir = str(tmp_path)
        cache = FileCache(data_dir)

        assert cache is not None

    def test_store_and_load_item_roundtrip(
        self,
        tmp_path,
    ):
        data_dir = str(tmp_path)
        cache = FileCache(data_dir)

        identifier = "whatever-id"
        item = Age.birth()

        # We start with clean cache ...
        assert len(os.listdir(data_dir)) == 0

        cache.store_item(identifier, item)

        # ... and now there is one item in the cache.
        vals = tuple(os.listdir(data_dir))
        assert len(vals) == 1
        val = vals[0]

        matcher = re.match(r"(?P<uid>[\w-]+).json", val)
        assert matcher is not None
        uid = matcher.group("uid")
        assert uuid.UUID(hex=uid)

        # Let's check we can get the item back.
        loaded = cache.load_item(identifier)
        assert loaded is not None
        assert loaded == item

        # Delete the cache ...
        del cache
        # ... and check that the file was created
        vals = tuple(os.listdir(data_dir))
        assert len(vals) == 2
        assert any(val == FileCache._MAPPING_FILE_NAME for val in vals)
