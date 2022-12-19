# recommender.py
"""
Recomender for running annotation predictions.
This module is going to be directly used by the users.
"""

# import collections
import compress_pickle
import itertools
import libsbml
import numpy as np
import os

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra



class Recommender(object):

  def __init__(self,
               libsbml_cl=None, 
               libsbml_fpath=None,
               model_specs=None):
    """
    Parameters
    ----------
    libsbml_cl: libsbml.SBMLDocument
        A libsbml document class instance
    libsbml_fpath: str
        File path of an SBML model
    mdoel_specs: tuple/list
        Iterable object of two tuples including model specifications
    """
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
      spec_tuple = None
      reac_tuple = None
    self.species = sa.SpeciesAnnotation(inp_tuple=spec_tuple)
    self.reactions = ra.ReactionAnnotation(inp_tuple=reac_tuple)

  # New version after discussion with Joe & Steve et al. 
  def getSpeciesAnnotation(self,
                           pred_str=None,
                           pred_id=None,
                           update=True,
                           method='cdist'):
    """
    Predict annotations of species using
    the provided string or ID.
    If pred_str is given, directly run the string;
    if pred_id is given, search for species ID and display name 
    to find useable name. 

    Parameters
    ----------
    pred_str: str
        Species name to predict annotation with
    pred_id: str
        ID of species (search for name using it)
    update: bool
        If true, update existing species annotations
        (i.e., replace or create new entries)
        in self.species.candidates and self.species.formula
    methood: str
        One of ['cdist', 'edist']
        'cdist' represents Cosine Similarity
        'edist' represents Edit Distance.
        Default method id 'cdist'

    Returns
    -------
    Recommendation (namedtuple)

    """
    if method == 'edist':
      if pred_str:
        name_to_use = pred_str
        given_id = pred_str
      elif pred_id:
        name_to_use = self.species.getNameToUse(inp_id=pred_id)
        given_id = pred_id
      pred_res = self.species.predictAnnotationByEditDistance(name_to_use)  
    elif method == 'cdist':
      if pred_str: 
        given_id = pred_str
        pred_res = self.species.predictAnnotationByCosineSimilarity(inp_strs=[pred_str])[pred_str]
      elif pred_id:
        given_id = pred_id
        pred_res = self.species.predictAnnotationByCosineSimilarity(inp_ids=[pred_id])[pred_id]
    #
    pred_score = self.species.evaluatePredictedSpeciesAnnotation(pred_result=pred_res)
    urls = [cn.CHEBI_DEFAULT_URL + val[6:] for val in pred_res[cn.CHEBI]]
    result = cn.Recommendation(given_id,
                               np.round(pred_score, 2),
                               pred_res[cn.MATCH_SCORE],
                               urls)
    if update:
      _ = self.species.updateSpeciesWithRecommendation(result)
    return result

  def getSpeciesListAnnotation(self,
                               pred_strs=None,
                               pred_ids=None,
                               update=True,
                               method='cdist'):
    """
    Get annotation of multiple species,
    given as a list (or an iterable object).
    self.getSpeciesAnnotation is applied to
    each element. 

    Parameters
    ----------
    pred_strs: str-list
        :Species names to predict annotations with
    pred_ids: str-list
        :Species IDs to predict annotations with
         (model info should have been already loaded)
    update: bool
        :If true, update the current annotations
        (i.e., replace or create new entries)
        in self.species.candidates and self.species.formula
    methood: str
        One of ['cdist', 'edist']
        'cdist' represents Cosine Similarity
        'edist' represents Edit Distance.
        Default method id 'cdist'

    Returns
    -------
    list-Recommendation (list-namedtuple)
    """
    if method == 'edist':
      if pred_strs:
        return [self.getSpeciesAnnotation(pred_str=val, update=update, method='edist') \
                for val in pred_strs]
      elif pred_ids:
        return [self.getSpeciesAnnotation(pred_id=val, update=update, method='edist') \
                for val in pred_ids]
    elif method == 'cdist':
      if pred_strs: 
        pred_res = self.species.predictAnnotationByCosineSimilarity(inp_strs=pred_strs)
      elif pred_ids: 
        pred_res = self.species.predictAnnotationByCosineSimilarity(inp_ids=pred_ids)
      res_recom = []
      for one_k in pred_res.keys():
        pred_score = self.species.evaluatePredictedSpeciesAnnotation(pred_result=pred_res[one_k])
        urls = [cn.CHEBI_DEFAULT_URL + val[6:] for val in pred_res[one_k][cn.CHEBI]]
        result = cn.Recommendation(one_k,
                                   np.round(pred_score, 2),
                                   pred_res[one_k][cn.MATCH_SCORE],
                                   urls)
        res_recom.append(result)
        if update:
          _ = self.species.updateSpeciesWithRecommendation(result)
      return res_recom

  def getReactionAnnotation(self, pred_id,
                            use_exist_species_annotation=False,
                            use_species_formula=None,
                            update=True):
    """
    Predict annotations of reactions using
    the provided IDs (argument). 
    Can be either singular (string) or plural

    Parameters
    ----------
    pred_id: str
        A single ID of reaction to annotate
    # TODO:
    use_exist_speices_annotation: bool
        If True, use existing species annotation
    # TODO: 
    Try the use_species_formula, so if it is given,
    use information from recom.species.formula

    Returns
    -------
    Recommendation (namedtuple)
    """
    # For now, just predict all species and continue? 
    specs2predict = self.reactions.reaction_components[pred_id] 
    if use_exist_species_annotation:
      pred_formulas = {val:self.species.exist_annotation_formula[val] \
                       for val in specs2predict \
                       if val in self.species.exist_annotation_formula.keys()}
    else:
      pred_formulas = {}
    remaining_species = [val for val in specs2predict if val not in pred_formulas.keys()]

    if len(remaining_species) > 0:
      spec_results = self.getSpeciesListAnnotation(pred_ids=remaining_species,
                                                   update=True)
      for one_recom in spec_results:
        chebis = [val[0] for val in one_recom.candidates]
        forms = list(set([cn.REF_CHEBI2FORMULA[k] \
                 for k in chebis if k in cn.REF_CHEBI2FORMULA.keys()]))
        pred_formulas[one_recom.id] = forms
    pred_reaction = self.reactions.predictAnnotation(inp_spec_dict=pred_formulas,
                                                     inp_reac_list=[pred_id],
                                                     update=update)
    pred_score = self.reactions.evaluatePredictedReactionAnnotation(pred_reaction)
    urls = [cn.RHEA_DEFAULT_URL + val[0][5:] for val in pred_reaction[cn.MATCH_SCORE][pred_id]]
    result = cn.Recommendation(pred_id,
                               np.round(pred_score[pred_id], 2),
                               pred_reaction[cn.MATCH_SCORE][pred_id],
                               urls)
    return result

  def getReactionListAnnotation(self, pred_ids,
                                use_exist_species_annotation=False,
                                update=True):
    """
    Get annotation of multiple reactions.
    Instead of applying getReactionAnnotation 
    for each reaction,
    it'll predict all component species first
    and proceed (this will reduce computational cost).

    Parameters
    ----------
    pred_ids: str-list
        For now, it only accommodates calling by reaction IDs.

    Returns
    -------
    list-Reccommendation (list-namedtuple)
    """
    # First, collect all species IDs to annotate
    specs_to_annotate = list(set(itertools.chain(*[self.reactions.reaction_components[val] \
                                                   for val in pred_ids])))

    if use_exist_species_annotation:
      pred_formulas = {val:self.species.exist_annotation_formula[val] \
                       for val in specs_to_annotate \
                       if val in self.species.exist_annotation_formula.keys()}
    else:
      pred_formulas = {}
    remaining_species = [val for val in specs_to_annotate if val not in pred_formulas.keys()]

    if len(remaining_species) > 0:
      spec_results = self.getSpeciesListAnnotation(pred_ids=remaining_species,
                                                   update=True)
      for one_recom in spec_results:
        chebis = [val[0] for val in one_recom.candidates]
        forms = list(set([cn.REF_CHEBI2FORMULA[k] \
                 for k in chebis if k in cn.REF_CHEBI2FORMULA.keys()]))
        pred_formulas[one_recom.id] = forms
    # Use predicted species in formula
    pred_reaction = self.reactions.predictAnnotation(inp_spec_dict=pred_formulas,
                                                     inp_reac_list=pred_ids,
                                                     update=update)
    pred_score = self.reactions.evaluatePredictedReactionAnnotation(pred_reaction)
    urls = {k:[cn.RHEA_DEFAULT_URL+val[0][5:] \
            for val in pred_reaction[cn.MATCH_SCORE][k]] \
            for k in pred_ids}
    result = [cn.Recommendation(k,
                               np.round(pred_score[k], 2),
                               pred_reaction[cn.MATCH_SCORE][k],
                               urls[k]) \
              for k in pred_score.keys()]
    return result 

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
    (tuple, tuple): 
        Two tuples to create species/reaction annotation classes
        (species_tuple, reaction_tuple)
    """
    if isinstance(sbml, str):
      reader = libsbml.SBMLReader()
      document = reader.readSBML(sbml)
    elif isinstance(sbml, libsbml.SBMLDocument):
      document = sbml
    model = document.getModel()
    # Create species_annotation instance
    # exist_spec_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
    #                              for val in model.getListOfSpecies()}
    # exist_spec_annotation_filt = {val:exist_spec_annotation_raw[val] for val in exist_spec_annotation_raw.keys() \
    #                               if exist_spec_annotation_raw[val] is not None}
    exist_spec_annotation = tools.extractExistingSpeciesAnnotation(model)
    species_names = {val.getId():val.name for val in model.getListOfSpecies()}
    species_tuple = (species_names, exist_spec_annotation)
    # Create reaction_annotation instance
    # Annotation of Rhea
    # reac_dict_raw_rhea = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.RHEA) \
    #                      for r in model.getListOfReactions()}
    # reac_dict_raw_filt_rhea = {k:reac_dict_raw_rhea[k] \
    #                            for k in reac_dict_raw_rhea.keys() \
    #                            if reac_dict_raw_rhea[k] is not None}
    # reac_dict_format_rhea = {k:[cn.RHEA_HEADER+val for val in reac_dict_raw_filt_rhea[k]] \
    #                            for k in reac_dict_raw_filt_rhea.keys()}
    # reac_dict_rhea = dict()
    # for one_id in reac_dict_format_rhea.keys():
    #   one_itm = list(set([cn.REF_RHEA2BI[val] for val in reac_dict_format_rhea[one_id] \
    #              if val in cn.REF_RHEA2BI.keys()]))
    #   if len(one_itm) > 0:
    #     reac_dict_rhea[one_id] = one_itm
    # # Annotation of KEGG (mapped to corresponding Rhea BI term) 
    # reac_dict_raw_kegg = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.KEGG_REACTION) \
    #                      for r in model.getListOfReactions()}
    # reac_dict_raw_filt_kegg = {k:reac_dict_raw_kegg[k] \
    #                            for k in reac_dict_raw_kegg.keys() \
    #                            if reac_dict_raw_kegg[k] is not None}
    # reac_dict_kegg = {k:[cn.REF_KEGG2RHEA_BI[val] \
    #                      for val in reac_dict_raw_filt_kegg[k] if val in cn.REF_KEGG2RHEA_BI.keys()] \
    #                   for k in reac_dict_raw_filt_kegg.keys()}
    # reac_dict_filt_kegg = {k: reac_dict_kegg[k] for k in reac_dict_kegg.keys() \
    #                           if reac_dict_kegg[k]}
    # reac_exist_annotation = reac_dict_rhea
    # for one_id in reac_dict_filt_kegg.keys():
    #   if one_id in reac_exist_annotation.keys():
    #     reac_exist_annotation[one_id] = list(set(reac_exist_annotation[one_id] + reac_dict_filt_kegg[one_id]))
    #   else:
    #     reac_exist_annotation[one_id] = list(set(reac_dict_filt_kegg[one_id]))

    reac_exist_annotation = tools.extractExistingReactionAnnotation(inp_model=model)
    # Next, reaction components for each reaction
    reac_components = {val.getId():list(set([k.species for k in val.getListOfReactants()]+[k.species for k in val.getListOfProducts()])) \
                       for val in model.getListOfReactions()}
    reaction_tuple = (reac_components, reac_exist_annotation)
    return species_tuple, reaction_tuple


  def getSpeciesStatistics(self):
    """
    Get recall and precision 
    of species in a model, for both species and
    reactions.
    This method works only if there exists
    annotation in  model; otherwise
    None will be returned. 
    In the result, values will be 
    returned after rounding to the two decimal places. 

    Parameters
    ----------
    None 

    Returns
    -------
    None/dict
        Return None if there is nothing to evaluate
        (i.e., if there is no existing model annotation)
    """
    # get dictionary of formulas if they exist
    refs = {val:self.species.exist_annotation_formula[val] \
            for val in self.species.exist_annotation_formula.keys() \
            if self.species.exist_annotation_formula[val]}
    specs2eval = list(refs.keys())
    if len(specs2eval) == 0:
      return None
    specsdict2pred = {val:self.species.getNameToUse(val) for val in specs2eval}
    preds = {val:self.species.predictAnnotationByEditDistance(specsdict2pred[val])[cn.FORMULA] \
             for val in specsdict2pred}
    recall = tools.getRecall(ref=refs, pred=preds, mean=True)
    precision = tools.getPrecision(ref=refs, pred=preds, mean=True)
    return {cn.RECALL: np.round(recall, 2), cn.PRECISION: np.round(precision, 2)}


  def getReactionStatistics(self):
    """
    Get recall and precision 
    of reactions in a model, for both species and
    reactions.
    This method works only if there exists
    annotation in  model; otherwise
    None will be returned. 
    In the result, values will be 
    returned after rounding to the two decimal places. 

    Parameters
    ----------
    None 

    Returns
    -------
    None/dict
        Return None if there is nothing to evaluate
        (i.e., if there is no existing model annotation)
    """
    # For reactions, component species should be
    # predicted first. 
    refs = self.reactions.exist_annotation
    specs2pred = list(set(itertools.chain(*([self.reactions.reaction_components[val] for val in refs.keys()]))))
    specsdict2pred = {val:self.species.getNameToUse(val) for val in specs2pred}
    specs_predicted = {val:self.species.predictAnnotationByEditDistance(specsdict2pred[val])[cn.FORMULA] \
                       for val in specs2pred}
    preds = self.reactions.predictAnnotation(inp_spec_dict=specs_predicted,
                                             inp_reac_list=refs.keys(),
                                             update=True)[cn.CANDIDATES]
    recall = tools.getRecall(ref=refs, pred=preds, mean=True)
    precision = tools.getPrecision(ref=refs, pred=preds, mean=True)
    return {cn.RECALL: np.round(recall, 2), cn.PRECISION: np.round(precision, 2)}




    