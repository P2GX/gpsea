"""
Information regarding exceptional execution outcomes, such as error states encountered during preprocessing.
"""

from ._var_coords import VariantCoordinateFinderExceptionData, VariantCoordinateFinderExceptionKind

__all__ = [
    "VariantCoordinateFinderExceptionKind",
    "VariantCoordinateFinderExceptionData",
]
