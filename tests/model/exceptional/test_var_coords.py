import json

from gpsea.io import GpseaJSONDecoder, GpseaJSONEncoder
from gpsea.model.exceptional import VariantCoordinateFinderExceptionData


class TestVariantCoordinateFinderException:
    def test_serialize_roundtrip(self):
        val = VariantCoordinateFinderExceptionData(
            kind="validation",
            context="validating ABC",
            issues=("Oh no", "Not again"),
        )

        encoded = json.dumps(val, cls=GpseaJSONEncoder, indent=2)
        assert (
            encoded
            == """{
  "kind": "validation",
  "context": "validating ABC",
  "issues": [
    "Oh no",
    "Not again"
  ]
}"""
        )
        decoded = json.loads(encoded, cls=GpseaJSONDecoder)

        assert decoded == val
