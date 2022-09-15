# recommender.py
# Recomender for running annotation predictions

import collections
import itertools
import libsbml
import numpy as np
import os
import pickle

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra

with open(os.path.join(cn.REF_DIR, 'kegg2rhea_bi.pickle'), 'rb') as handle:
  ref_kegg2rhea_bi = pickle.load(handle)
with open(os.path.join(cn.REF_DIR, 'rhea_all2bi.pkl'), 'rb') as f:
  ref_rhea2bi = pickle.load(f)


Recommendation = collections.namedtuple('Recommendation',
                                        ['id', 'credibility_score', 'candidates', 'urls'])


class Recommender(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier=cn.RHEA):
    # First of all, collect model information from libsbml model
    # and send the informaton to create species/reaction annotations
    if libsbml_fpath:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      # Create species_annotation instance
      exist_spec_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                                   for val in self.model.getListOfSpecies()}
      exist_spec_annotation_filt = {val:exist_spec_annotation_raw[val] for val in exist_spec_annotation_raw.keys() \
                                    if exist_spec_annotation_raw[val] is not None}
      species_names = {val.getId():val.name for val in self.model.getListOfSpecies()}
      species_tuple = (species_names, exist_spec_annotation_filt)
      self.species = sa.SpeciesAnnotation(inp_tuple=species_tuple)
      # Create reaction_annotation instance
      # Annotation of Rhea
      reac_dict_raw_rhea = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.RHEA) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_rhea = {k:reac_dict_raw_rhea[k] \
                                 for k in reac_dict_raw_rhea.keys() \
                                 if reac_dict_raw_rhea[k] is not None}
      reac_dict_format_rhea = {k:['RHEA:'+val for val in reac_dict_raw_filt_rhea[k]] \
                                 for k in reac_dict_raw_filt_rhea.keys()}
      reac_dict_rhea = dict()
      for one_id in reac_dict_format_rhea.keys():
        one_itm = list(set([ref_rhea2bi[val] for val in reac_dict_format_rhea[one_id] \
                   if val in ref_rhea2bi.keys()]))
        if len(one_itm) > 0:
          reac_dict_rhea[one_id] = one_itm
      # Annotation of KEGG (mapped to corresponding Rhea BI term) 
      reac_dict_raw_kegg = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.KEGG_REACTION) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_kegg = {k:reac_dict_raw_kegg[k] \
                                 for k in reac_dict_raw_kegg.keys() \
                                 if reac_dict_raw_kegg[k] is not None}
      reac_dict_kegg = {k:[ref_kegg2rhea_bi[val] \
                           for val in reac_dict_raw_filt_kegg[k] if val in ref_kegg2rhea_bi.keys()] \
                        for k in reac_dict_raw_filt_kegg.keys()}
      reac_exist_annotation = reac_dict_rhea
      for one_id in reac_dict_kegg.keys():
        if one_id in reac_exist_annotation.keys():
          reac_exist_annotation[one_id] = list(set(reac_exist_annotation[one_id] + reac_dict_kegg[one_id]))
        else:
          reac_exist_annotation[one_id] = list(set(reac_dict_kegg[one_id]))
      # Next, reaction components for each reaction
      reac_components = {val.getId():list(set([k.species for k in val.getListOfReactants()]+[k.species for k in val.getListOfProducts()])) \
                         for val in self.model.getListOfReactions()}
      reaction_tuple = (reac_components, reac_exist_annotation)
      self.reactions = ra.ReactionAnnotation(inp_tuple=reaction_tuple)


  def getSpeciesAnnotation(self, name_to_annotate):
    """
    Predict annotations of species using
    the provided IDs (argument).
    Can be a singuler (string) or a list of
    strings. 

    Parameters
    ----------
    name_to_annotate: str/list-str
        ID of species to annotate

    Returns
    -------
    result: Recommendation (namedtuple)

    """
    if isinstance(name_to_annotate, str):
      inp_list = [name_to_annotate]
    else:
      inp_list = name_to_annotate

    pred_result = self.species.predictAnnotationByName(inp_list)
    pred_score = self.species.evaluatePredictedSpeciesAnnotation(inp_list)
    urls = {k:['https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A'+val[6:] \
            for val in pred_result[k][cn.CHEBI]] \
            for k in inp_list}
    result = [Recommendation(k,
                             np.round(pred_score[k], 2),
                             pred_result[k][cn.MATCH_SCORE],
                             urls[k]) \
              for k in pred_score.keys()]
    return result

  def getReactionAnnotation(self, name_to_annotate):
    """
    Predict annotations of reactions using
    the provided IDs (argument). 
    Can be either singular (string) or plural

    Parameters
    ----------
    name_to_annotate: str/list-str
        ID of reactions to annotate

    Returns
    -------
    result: Recommendation (namedtuple)

    """
    if isinstance(name_to_annotate, str):
      inp_list = [name_to_annotate]
    else:
      inp_list = name_to_annotate
    # First, collect all species IDs to annotate
    specs_to_annotate = list(set(itertools.chain(*[self.reactions.reaction_components[val] \
                                                   for val in inp_list])))
    # For now, just predict all species and continue? 
    spec_results = self.getSpeciesAnnotation(specs_to_annotate)
    pred_formulas = self.species.formula
    # Use predicted species in formula
    pred_reaction = self.reactions.predictAnnotation(inp_spec_dict=pred_formulas,
                                                     inp_reac_list=inp_list)
    pred_score = self.reactions.evaluatePredictedReactionAnnotation(inp_list)
    urls = {k:['https://www.rhea-db.org/rhea/'+val[0][5:] \
            for val in pred_reaction[k]] \
            for k in inp_list}
    result = [Recommendation(k,
                             np.round(pred_score[k], 2),
                             pred_reaction[k],
                             urls[k]) \
              for k in pred_score.keys()]
    return result


  def updateSpeciesAnnotation(self, update_dict):
    """
    Update annotation of species; 

    Parameters
    ----------

    Returns
    -------


    """
    pass


  def updateReactionAnnotation(self, update_dict):
    """
    Update annotation of reactions; 

    Parameters
    ----------

    Returns
    -------


    """
    pass


  def reportAnnotation(self):
    """
    Create (and save) a report (or table)
    summarizing predicted anntoations.
    """









    