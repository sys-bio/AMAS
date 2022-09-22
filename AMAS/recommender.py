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
                                        ['id', 'credibility', 'candidates', 'urls'])


class Recommender(object):

  def __init__(self,
               libsbml_cl=None, 
               libsbml_fpath=None,
               model_specs=None):
    # First of all, collect model information from libsbml model
    # and send the informaton to create species/reaction annotations
    if libsbml_cl:
      spec_tuple, reac_tuple = self._parseSBML(libsbml_cl)
    elif libsbml_fpath:
      spec_tuple, reac_tuple = self._parseSBML(libsbml_fpath)
    elif model_specs:
      spec_tuple = model_specs[0]
      reac_tuple = model_specs[1]
    else:
      return None
    self.species = sa.SpeciesAnnotation(inp_tuple=spec_tuple)
    self.reactions = ra.ReactionAnnotation(inp_tuple=reac_tuple)

  ## Previous version; 
  # def getSpeciesAnnotation(self, name_to_annotate):
  #   """
  #   Predict annotations of species using
  #   the provided IDs (argument).
  #   Can be a singuler (string) or a list of
  #   strings. 

  #   Parameters
  #   ----------
  #   name_to_annotate: str/list-str
  #       ID of species to annotate

  #   Returns
  #   -------
  #   result: Recommendation (namedtuple)

  #   """
  #   if isinstance(name_to_annotate, str):
  #     inp_list = [name_to_annotate]
  #   else:
  #     inp_list = name_to_annotate

  #   pred_result = self.species.predictAnnotationByName(inp_list)
  #   pred_score = self.species.evaluatePredictedSpeciesAnnotation(inp_list)
  #   urls = {k:['https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A'+val[6:] \
  #           for val in pred_result[k][cn.CHEBI]] \
  #           for k in inp_list}
  #   result = [Recommendation(k,
  #                            np.round(pred_score[k], 2),
  #                            pred_result[k][cn.MATCH_SCORE],
  #                            urls[k]) \
  #             for k in pred_score.keys()]
  #   return result

  # New version after discussion with Joe & Steve et al. 
  def getSpeciesAnnotation(self, pred_str=None, pred_id=None):
    """
    Predict annotations of species using
    the provided string or ID.
    If pred_str is given, directly run the string;
    if pred_id is given, search for species ID and display name 
    to find useable name. 

    Parameters
    ----------
    pred_str: str
        Species name to predict annotation withh
    pred_id: str
        ID of species (search for name using it)

    Returns
    -------
    result: Recommendation (namedtuple)

    """
    if pred_str:
      name_to_use = pred_str
      given_id = pred_str
    elif pred_id:
      name_to_use = self.species.getNameToUse(inp_id=pred_id)
      given_id = pred_id
    pred_res = self.species.predictAnnotationByEditDistance(name_to_use)  
    pred_score = self.species.evaluatePredictedSpeciesAnnotation(pred_result=pred_res)
    urls = [cn.CHEBI_DEFAULT_URL + val[6:] for val in pred_res[cn.CHEBI]]
    result = Recommendation(given_id,
                            np.round(pred_score, 2),
                            pred_res[cn.MATCH_SCORE],
                            urls)
    return result

  def getSpeciesListAnnotation(self, pred_list, id=False):
    """
    Get annotation of multiple species,
    given as a list (or an iterable object).
    self.getSpeciesAnnotation is applied to
    each element. 

    Parameters
    ----------
    pred_list: str-list
        :List of predicrted annotation
    id: bool
        :Indicator whether given strings are ids or direct names

    Returns
    -------
    list-Recommendation (list-namedtuple)
    """
    if id:
      return [self.getSpeciesAnnotation(pred_id=val) \
              for val in pred_list]
    else:
      return [self.getSpeciesAnnotation(pred_str=val) \
              for val in pred_list]

  def getReactionAnnotation(self, pred_id):
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
    # For now, just predict all species and continue? 
    specs2predict = self.reactions.reaction_components[pred_id]
    spec_results = self.getSpeciesListAnnotation(specs2predict, id=True)
    # based on the function above; need to recreate it. 
    pred_formulas = dict()
    for one_recom in spec_results:
      cands = [val[0] for val in one_recom.candidates]
      forms = list(set([sa.ref_shortened_chebi_to_formula[k] \
               for k in cands if k in sa.ref_shortened_chebi_to_formula.keys()]))
      pred_formulas[one_recom.id] = forms
    #
    pred_reaction = self.reactions.predictAnnotation(inp_spec_dict=pred_formulas,
                                                     inp_reac_list=[pred_id])
    pred_score = self.reactions.evaluatePredictedReactionAnnotation([pred_id])
    urls = [cn.RHEA_DEFAULT_URL + val[0][5:] for val in pred_reaction[pred_id]]
    result = Recommendation(pred_id,
                            np.round(pred_score[pred_id], 2),
                            pred_reaction[pred_id],
                            urls)
    return result

  def getReactionListAnnotation(self, pred_list):
    """
    Get annotation of multiple reactions.
    Instead of applying getReactionAnnotation 
    for each reaction,
    it'll predict all component species first
    and proceed (this will reduce computational cost).

    Parameters
    ----------
    pred_list: str-list

    Returns
    -------
    list-Reccommendation (list-namedtuple)
    """
    # First, collect all species IDs to annotate
    specs_to_annotate = list(set(itertools.chain(*[self.reactions.reaction_components[val] \
                                                   for val in pred_list])))
    # For now, just predict all species and continue? 
    spec_results = self.getSpeciesListAnnotation(specs_to_annotate, id=True)
    # pred_formulas = self.species.formula
    pred_formulas = dict()
    for one_recom in spec_results:
      cands = [val[0] for val in one_recom.candidates]
      forms = list(set([sa.ref_shortened_chebi_to_formula[k] \
               for k in cands if k in sa.ref_shortened_chebi_to_formula.keys()]))
      pred_formulas[one_recom.id] = forms
    # Use predicted species in formula
    pred_reaction = self.reactions.predictAnnotation(inp_spec_dict=pred_formulas,
                                                     inp_reac_list=pred_list)
    pred_score = self.reactions.evaluatePredictedReactionAnnotation(pred_list)
    urls = {k:[cn.RHEA_DEFAULT_URL+val[0][5:] \
            for val in pred_reaction[k]] \
            for k in pred_list}
    result = [Recommendation(k,
                             np.round(pred_score[k], 2),
                             pred_reaction[k],
                             urls[k]) \
              for k in pred_score.keys()]
    return result
    # return [self.getReactionAnnotation(pred_id=val) \
    #         for val in pred_list]    


  def _parseSBML(self, sbml):
    """
    Parse SBML file and return 
    two tuples, for species and reactions 
    respecitvely,
    equivalent to model_specs.
    Can use either libsbml.Document class
    or a file path.

    Parameters
    ----------
    sbml: str(filepath)/libsbml.Document

    Returns
    -------
    species_tuple, reaction_tuple: 
        Two tuples to create species/reaction annotation classes
    """
    if isinstance(sbml, str):
      reader = libsbml.SBMLReader()
      document = reader.readSBML(sbml)
    elif isinstance(sbml, libsbml.SBMLDocument):
      document = sbml
    model = document.getModel()
    # Create species_annotation instance
    exist_spec_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                                 for val in model.getListOfSpecies()}
    exist_spec_annotation_filt = {val:exist_spec_annotation_raw[val] for val in exist_spec_annotation_raw.keys() \
                                  if exist_spec_annotation_raw[val] is not None}
    species_names = {val.getId():val.name for val in model.getListOfSpecies()}
    species_tuple = (species_names, exist_spec_annotation_filt)
    # Create reaction_annotation instance
    # Annotation of Rhea
    reac_dict_raw_rhea = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.RHEA) \
                         for r in model.getListOfReactions()}
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
                         for r in model.getListOfReactions()}
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
                       for val in model.getListOfReactions()}
    reaction_tuple = (reac_components, reac_exist_annotation)
    return species_tuple, reaction_tuple


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









    