import typing
from dataclasses import dataclass

VariantCoordinateFinderExceptionKind = typing.Literal["validation", "timeout"]

@dataclass(frozen=True, slots=True)
class VariantCoordinateFinderExceptionData:
    kind: VariantCoordinateFinderExceptionKind | None
    context: str
    issues: typing.Sequence[str] | None
