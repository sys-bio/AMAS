"""Annotation for Species."""

from annotation_recommender import constants as cn
from annotation_recommender import tools

import editdistance
import libsbml
import numpy as np
import os
import pickle

# below might be in constants or main script
with open(os.path.join(cn.REF_DIR, 'chebi_shortened_formula_30apr2022.pickle'), 'rb') as f:
  ref_shortened_chebi_to_formula = pickle.load(f)
with open(os.path.join(cn.REF_DIR, 'chebi_synonyms.pickle'), 'rb') as f:
  chebi_synonyms = pickle.load(f)
#
chebi_low_synonyms = dict()
for one_k in chebi_synonyms.keys():
  chebi_low_synonyms[one_k] = list(set([val.lower() for val in chebi_synonyms[one_k]]))

species_rf = pickle.load(open(os.path.join(cn.REF_DIR, 'species_randomforestcv.sav'), 'rb'))


class SpeciesAnnotation(object):

  def __init__(self, libsbml_fpath=None):
    # self.exist_annotation stores existing CHEBI annotations in the model
    # If none exists, set None
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      exist_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                        for val in self.model.getListOfSpecies()}
      exist_annotation_filt = {val:exist_annotation_raw[val] for val in exist_annotation_raw.keys() \
                               if exist_annotation_raw[val] is not None}
      self.exist_annotation = {k:tools.transformCHEBIToFormula(exist_annotation_filt[k], ref_shortened_chebi_to_formula) \
                               for k in exist_annotation_filt.keys()}
    else:
      self.model = None
      self.exist_annotation = None
    # Below are predicted annotations;
    # Once created, each will be {species_ID: float/str-list}
    self.match_score = None
    self.candidates = None
    self.formula = None
      

  def predictAnnotationByEditDistance(self, inp_str):
    """
    Predict annotation using the argument string 
    and Levenshtein edit distance method. 

    Parameters
    ----------
    inp_str: str
        String to predict CHEBI annotation

    Returns
    -------
    one_result: dict
        {key: value(s)}
        'key' can be match_score, chebi, etc. 
    """
    one_result = dict()
    # For now, choose the terms that are included in the CHEBI-formula mapping reference
    dist_dict_min = {one_k:np.min([editdistance.eval(inp_str.lower(), val) for val in chebi_low_synonyms[one_k]]) \
                     for one_k in chebi_low_synonyms.keys() if one_k in ref_shortened_chebi_to_formula.keys()}
    min_min_dist = np.min([dist_dict_min[val] for val in dist_dict_min.keys()])
    one_match_score = 1 - min_min_dist/len(inp_str)
    one_result[cn.MATCH_SCORE] = one_match_score
    min_min_chebis = [one_k for one_k in dist_dict_min.keys() \
                      if dist_dict_min[one_k]==min_min_dist and one_k in ref_shortened_chebi_to_formula.keys()]
    # predicted formula of the species
    one_result[cn.CHEBI] = min_min_chebis
    min_min_formula = list(set([ref_shortened_chebi_to_formula[val] for val in min_min_chebis]))
    one_result[cn.FORMULA] = min_min_formula
    return one_result


  def predictAnnotationByName(self, inp_spec_list=None,
                              specnames_dict=None,
                              update=True):
    """
    Predict list of species annotations
    using species names/IDs.
    Rule is 1) use species name, 
    2) if not provided, use species ID.

    Alternatively, user can directly predict annotations
    without incurring libsbml model methods,
    by using specnames_dict.
  
    Parameters
    ----------
    inp_spec_list: str-list (or iterable list of strings)
        List of species IDs to extract names
    specnames_dict: dict
        {spec_id: spec_name(str)}

    Returns
    -------
    result:  
        Should be {species_id: 
                      {'chebi':[CHEBI terms]},
                      {'score': match_score},
                      {'formula': [chemical formulas in string]}
                      {'formula2chebi': [CHEHBI terms]}
                  }
        match_score is expected to be between 0.0-1.0
    """
    if specnames_dict is None:
      result = dict()
      # If no list if given, predict all elements' annotations
      if inp_spec_list is None:
        spec_list = [val.getId() for val in self.model.getListOfSpecies()]
      else:
        spec_list = inp_spec_list
      for one_spec_id in spec_list:
        one_spec_name = self.model.getSpecies(one_spec_id).name.lower()
        if len(one_spec_name) == 0:
          one_spec_name = one_spec_id.lower()
        result[one_spec_id] = self.predictAnnotationByEditDistance(inp_str=one_spec_name)
    else:
      result = {val:self.predictAnnotationByEditDistance(inp_str=specnames_dict[val]) \
                for val in specnames_dict.keys()}      
    #
    if update:
      self.match_score = {spec_id:result[spec_id][cn.MATCH_SCORE] for spec_id in result.keys()}
      self.candidates = {spec_id:result[spec_id][cn.CHEBI] for spec_id in result.keys()}
      self.formula = {spec_id:result[spec_id][cn.FORMULA] for spec_id in result.keys()}
    return result


  def getAccuracy(self,
                  ref_annotation=None,
                  pred_annotation=None):
    """
    Compute accuracy of species annotation.
    A list of annotations of 
    a single species (identified by each ID) 
    is considered accurate if it includes
    the corresponding value of ref_annotation.
    (More precisely, if there is at least one
    intersection).
  
    Parameters
    ----------
    ref_annotation: dict
        {species_id: [str-annotation]}
        if None, get self.exist_annotation
    pred_annotation: dict
        {species_id: [str-annotation]}
        if None, get self.candidates

    Returns
    -------
    : float
    """
    accuracy = []
    if ref_annotation is None:
      ref = self.exist_annotation
    else:
      ref = ref_annotation
    if pred_annotation is None:
      pred = self.formula
    else:
      pred = pred_annotation
    ref_keys = set(ref.keys())
    pred_keys = set(pred.keys())
    species_to_test = ref_keys.intersection(pred_keys)
    for one_k in species_to_test:
      if set(ref[one_k]).intersection(pred[one_k]):
        accuracy.append(True)
      else:
        accuracy.append(False)
    return np.mean(accuracy)

  def getNameToUse(self, inp_id):
    """
    Get name to use;
    If .name is not '', use it;
    otherwise use ID
  
    Parameters
    ----------
    inp_id: ID of model element
  
    Returns
    -------
    res_name: str
    """
    one_species = self.model.getSpecies(inp_id)
    species_name = one_species.name
    if len(species_name) > 0:
      res_name = species_name
    else:
      res_name = inp_id
    return res_name

  def evaluatePredictedSpeciesAnnotation(self, inp_list,
                                         fitted_model=species_rf):
    """
    Evaluate the quality of annotation;
    for each individual species.
    
    Parameters
    ---------
    inp_list: str-list?
        List of species IDs to evaluate

    Returns
    -------
    res: dict {species_id: level-of-species-prediction-being-correct}
        Information of whether confident or not
    """
    name_lengths = [len(self.getNameToUse(inp_id=val)) for val in inp_list]
    nums_candidates = [len(self.candidates[val]) for val in inp_list]
    match_scores = [self.match_score[val] for val in inp_list]
    data2prediction = list(zip(name_lengths, nums_candidates, match_scores))
    # # loaded_model is loaded fitted logistic regression CV model
    # pred_probs = [val[1] for val in fitted_model.predict(data2prediction)]
    # Collect probability to be correct
    res = {val[0]:val[1] for val in list(zip(inp_list, fitted_model.predict(data2prediction)))}
    return res


