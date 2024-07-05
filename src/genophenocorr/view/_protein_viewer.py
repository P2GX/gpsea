import typing

from jinja2 import Environment, PackageLoader
from collections import OrderedDict

from genophenocorr.model import Cohort
from genophenocorr.model.genome import Region
from ._protein_visualizable import ProteinVisualizable


class ProteinViewable:
    """
    TODO
    """
    def __init__(self) -> None:
        environment = Environment(loader=(PackageLoader('genophenocorr.view', 'templates')))
        self._cohort_template = environment.get_template("protein.html")
        
    def process(self, cohort: Cohort, pvis: ProteinVisualizable) -> str:
        context = self._prepare_context(cohort, pvis)
        return self._cohort_template.render(context)
    
    @staticmethod
    def _get_start(feat_list: typing.Iterable[typing.Any]) -> int:
        return feat_list[2].start
    
    def _prepare_context(self, cohort: Cohort, pvis: ProteinVisualizable) -> typing.Mapping[str, typing.Any]:
        protein_id = pvis.protein_id
        protein_features = []
        
        for i in range(len(pvis.protein_feature_names)):
            protein_features.append([pvis.protein_feature_names[i], pvis._protein_feature_types[i], Region(pvis.protein_feature_starts[i],pvis.protein_feature_ends[i]), 0])
        
        protein_features = sorted(protein_features, key=self._get_start)
            
        for var in cohort.all_variants():
            tx_anno = var.get_tx_anno_by_tx_id(pvis.transcript_id)
            if tx_anno is not None:
                for feat_list in protein_features:
                    if tx_anno.protein_effect_location is not None and tx_anno.protein_effect_location.overlaps_with(feat_list[2]):
                        feat_list[3] += 1
            
        return {
            'protein_id': protein_id,
            'protein_label': pvis.protein_metadata.label,
            'protein_features': protein_features
        }
    
