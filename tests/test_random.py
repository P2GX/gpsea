# Tests of random stuff asked for by the users.
import collections.abc

import hpotk
import pytest
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket

from gpsea.preprocessing import CohortCreator, configure_caching_cohort_creator_context, load_phenopacket_files


@pytest.fixture(scope="module")
def hpo() -> hpotk.MinimalOntology:
    store = hpotk.configure_ontology_store()
    return store.load_minimal_hpo(release="v2025-10-22")


@pytest.fixture(scope="module")
def cohort_creator(
    hpo: hpotk.MinimalOntology,
) -> collections.abc.Generator[CohortCreator[Phenopacket]]:
    with configure_caching_cohort_creator_context(hpo) as cohort_creator:
        yield cohort_creator


@pytest.mark.skip(reason="Just for interactive debugging")
def test_load_phenopacket(
    cohort_creator: CohortCreator,
):
    pps = ("dev/Mito/1-10011778Ff.json",)
    _cohort, qc = load_phenopacket_files(
        pp_files=iter(pps),
        cohort_creator=cohort_creator,
    )
    qc.summarize()
