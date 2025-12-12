# Tests of random stuff asked for by the users.

import hpotk
import pytest

from gpsea.preprocessing import CohortCreator, configure_caching_cohort_creator, load_phenopacket_files


@pytest.fixture(scope="module")
def hpo() -> hpotk.MinimalOntology:
    store = hpotk.configure_ontology_store()
    return store.load_minimal_hpo(release="v2025-10-22")


@pytest.fixture(scope="module")
def cohort_creator(
    hpo: hpotk.MinimalOntology,
) -> CohortCreator:
    return configure_caching_cohort_creator(
        hpo=hpo,
    )


@pytest.mark.skip(reason="Just for interactive debugging")
def test_load_phenopacket(
    cohort_creator: CohortCreator,
):
    pps = ("dev/Mito/1-10011778Ff.json",)
    _cohort, qc = load_phenopacket_files(
        pp_files=pps,
        cohort_creator=cohort_creator,
    )
    qc.summarize()
