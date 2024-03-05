import pytest

from pkg_resources import resource_filename

from google.protobuf.json_format import Parse

# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket, GenomicInterpretation

from genophenocorr.model.genome import GenomeBuild, GRCh38

from ._api import VariantCoordinateFinder
from ._phenopacket import PhenopacketVariantCoordinateFinder
from ._variant import VariantAnnotationCache, VarCachingFunctionalAnnotator
from ._vep import VepFunctionalAnnotator
from ._vv import VVHgvsVariantCoordinateFinder


@pytest.fixture
def build() -> GenomeBuild:
    return GRCh38


@pytest.fixture
def hgvs_vc_finder(build: GenomeBuild) -> VariantCoordinateFinder:
    return VVHgvsVariantCoordinateFinder(build)


@pytest.fixture
def pp_vc_finder(build: GenomeBuild,
                 hgvs_vc_finder: VariantCoordinateFinder) -> PhenopacketVariantCoordinateFinder:
    return PhenopacketVariantCoordinateFinder(build, hgvs_vc_finder)


@pytest.mark.online
@pytest.mark.parametrize("pp_path, expected",
                         [('test_data/deletion_test.json', '16_89284129_89284134_CTTTTT_C'),
                          ('test_data/insertion_test.json', '16_89280829_89280829_C_CA'),
                          ('test_data/missense_test.json', '16_89279135_89279135_G_C'),
                          ('test_data/missense_hgvs_test.json', '16_89279135_89279135_G_C'),
                          ('test_data/duplication_test.json', '16_89279850_89279850_G_GC'),
                          ('test_data/delinsert_test.json', '16_89284601_89284602_GG_A'),
                          ('test_data/CVDup_test.json', '16_89284523_89373231_DUP'),
                          ('test_data/CVDel_test.json', '16_89217281_89506042_DEL')
                          ])
def test_find_coordinates(pp_path, expected, pp_vc_finder):
    fname = resource_filename(__name__, pp_path)
    gi = read_genomic_interpretation_json(fname)

    vc, gt = pp_vc_finder.find_coordinates(gi)

    assert vc.variant_key == expected


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())

@pytest.fixture
def variant_annotator():
    return VepFunctionalAnnotator()

@pytest.fixture
def caching_annotator(variant_annotator, tmp_path):
    vac = VariantAnnotationCache(tmp_path)
    return VarCachingFunctionalAnnotator(vac, variant_annotator)

@pytest.mark.online
def test_caching_full_circle(caching_annotator, pp_vc_finder, variant_annotator):
    fname = resource_filename(__name__, 'test_data/missense_test.json')
    gi = read_genomic_interpretation_json(fname)
    var_coords, gt = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cache_anno_results = caching_annotator.annotate(var_coords)
    assert var_anno_results == cache_anno_results
    assert cache_anno_results == caching_annotator.annotate(var_coords)

@pytest.fixture
def oldfile_cache_annotator(variant_annotator):
    data = resource_filename(__name__, 'test_data/annotations')
    vac = VariantAnnotationCache(data)
    return VarCachingFunctionalAnnotator(vac, variant_annotator)


@pytest.mark.skip("Skipped until new Pickle files are generated")
def test_cache_from_older_file(oldfile_cache_annotator, pp_vc_finder, variant_annotator):
    fname = resource_filename(__name__, 'test_data/missense_test.json')
    gi = read_genomic_interpretation_json(fname)
    var_coords, gt = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cached_file_results = oldfile_cache_annotator.annotate(var_coords)
    assert var_anno_results == cached_file_results
