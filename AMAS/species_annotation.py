# species_annotation.py
"""
<Annotation for Species>
species_annotation creates and predicts
annotation of libsbml species,
mainly using editdistance method.
(can be changed in future)
"""


from AMAS import constants as cn
from AMAS import tools

import compress_pickle
import editdistance
import libsbml
import numpy as np
import operator
import os
import pickle


with open(os.path.join(cn.REF_DIR, 'chebi_low_synonyms_comp.lzma'), 'rb') as f:
  CHEBI_LOW_SYNONYMS = compress_pickle.load(f)
# A trained random forest model 
SPECIES_RF = pickle.load(open(os.path.join(cn.REF_DIR, 'species_randomforestcv.sav'), 'rb'))


class SpeciesAnnotation(object):

  def __init__(self, libsbml_fpath=None,
               inp_tuple=None):

    """
    Parameters
    ----------
    libsbml_fpath: str
        File path of an SBMl (.xml) model

    inp_tuple: tuple 
        Tuple of model information,
        first element (index 0) is information on 
        species names,
        second element (index 1) is existing 
        ChEBI information.
        ({species_id: species_display_name},
         {species_id: ChEBI terms})
    """

    # self.exist_annotation stores existing CHEBI annotations in the model
    # If none exists, set None
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      self.names = {val.getId():val.name for val in self.model.getListOfSpecies()}
      exist_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                        for val in self.model.getListOfSpecies()}
      exist_annotation_chebi = {val:exist_annotation_raw[val] for val in exist_annotation_raw.keys() \
                               if exist_annotation_raw[val] is not None}
      self.exist_annotation = exist_annotation_chebi
      self.exist_annotation_formula = {k:tools.transformCHEBIToFormula(exist_annotation_chebi[k], cn.REF_CHEBI2FORMULA) \
                                       for k in exist_annotation_chebi.keys()}
    # inp_tuple: ({species_id:species_name}, {species_id: [CHEBI annotations]})
    elif inp_tuple is not None:
      self.model = None
      self.names = inp_tuple[0]
      self.exist_annotation = inp_tuple[1]
      self.exist_annotation_formula = {k:tools.transformCHEBIToFormula(inp_tuple[1][k], cn.REF_CHEBI2FORMULA) \
                                       for k in inp_tuple[1].keys()}
    else:
      self.model = None
      self.names = None
      self.exist_annotation = None
      self.exist_annotation_formula = None
    # Below are predicted annotations in dictionary format
    # Once created, each will be {species_ID: float/str-list}
    # self.match_score = dict()
    self.candidates = dict()
    self.formula = dict()
      
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
    dict
        {key: value(s)}
        'key' can be match_score, chebi, etc. 
    """
    one_result = dict()
    # For now, choose the terms that are included in the CHEBI-formula mapping reference
    dist_dict_min = {one_k:np.min([editdistance.eval(inp_str.lower(), val) for val in CHEBI_LOW_SYNONYMS[one_k]]) \
                     for one_k in CHEBI_LOW_SYNONYMS.keys() if one_k in cn.REF_CHEBI2FORMULA.keys()}
    min_min_dist = np.min([dist_dict_min[val] for val in dist_dict_min.keys()])
    min_min_chebis = [one_k for one_k in dist_dict_min.keys() \
                      if dist_dict_min[one_k]==min_min_dist and one_k in cn.REF_CHEBI2FORMULA.keys()]
    # Results are sorted based on match_score (average of 1 - (editdistance/len_synonyms)
    res_tuple = [(one_chebi,
                  np.round(np.max([1.0-editdistance.eval(inp_str.lower(), val)/len(val) \
                                    for val in CHEBI_LOW_SYNONYMS[one_chebi]]), 2)) \
                 for one_chebi in min_min_chebis] 
    res_tuple.sort(key=operator.itemgetter(1), reverse=True)
    #  CHEBI part is added, because we want a sorted list after computing res_tuple
    one_result[cn.NAME_USED] = inp_str
    one_result[cn.CHEBI] = [val[0] for val in res_tuple]
    one_result[cn.MATCH_SCORE] = res_tuple
    min_min_formula = list(set([cn.REF_CHEBI2FORMULA[val] for val in min_min_chebis]))
    one_result[cn.FORMULA] = min_min_formula
    return one_result


  # def predictAnnotationByName(self, inp_spec_list=None,
  #                             specnames_dict=None,
  #                             update=True):
  #   """
  #   Predict list of species annotations
  #   using species names/IDs.
  #   Rule is 1) use species name, 
  #   2) if not provided, use species ID.

  #   Alternatively, user can directly predict annotations
  #   without incurring libsbml model methods,
  #   by using specnames_dict.
  
  #   Parameters
  #   ----------
  #   inp_spec_list: str-list (or iterable list of strings)
  #       List of species IDs to extract names
  #   specnames_dict: dict
  #       {spec_id: spec_name(str)}

  #   Returns
  #   -------
  #   dict  
  #       Should be {species_id: 
  #                     {'chebi':[CHEBI terms]},
  #                     {'score': match_score},
  #                     {'formula': [chemical formulas in string]}
  #                     {'formula2chebi': [CHEHBI terms]}
  #                 }
  #       match_score is expected to be between 0.0-1.0
  #   """
  #   if specnames_dict is None:
  #     result = dict()
  #     # If no list if given, predict all elements' annotations
  #     if inp_spec_list is None:
  #       spec_list = list(self.names)
  #       # spec_list = [val.getId() for val in self.model.getListOfSpecies()]
  #     else:
  #       spec_list = inp_spec_list
  #     for one_spec_id in spec_list:
  #       one_spec_name = self.names[one_spec_id].lower()
  #       # one_spec_name = self.model.getSpecies(one_spec_id).name.lower()
  #       if len(one_spec_name) == 0:
  #         one_spec_name = one_spec_id.lower()
  #       result[one_spec_id] = self.predictAnnotationByEditDistance(inp_str=one_spec_name)
  #   else:
  #     result = {val:self.predictAnnotationByEditDistance(inp_str=specnames_dict[val]) \
  #               for val in specnames_dict.keys()}      
  #   # Might separate method or just leave it?
  #   if update:
  #     self.match_score.update({spec_id:result[spec_id][cn.MATCH_SCORE] for spec_id in result.keys()})
  #     self.candidates.update({spec_id:result[spec_id][cn.CHEBI] for spec_id in result.keys()})
  #     self.formula.update({spec_id:result[spec_id][cn.FORMULA] for spec_id in result.keys()})
  #   return result


  # def getAccuracy(self,
  #                 ref_annotation=None,
  #                 pred_annotation=None):
  #   """
  #   Compute accuracy of species annotation.
  #   A list of annotations of 
  #   a single species (identified by each ID) 
  #   is considered accurate if it includes
  #   the corresponding value of ref_annotation.
  #   (More precisely, if there is at least one
  #   intersection).
  
  #   Parameters
  #   ----------
  #   ref_annotation: dict
  #       {species_id: [str-annotation]}
  #       if None, get self.exist_annotation_formula
  #   pred_annotation: dict
  #       {species_id: [str-annotation]}
  #       if None, get self.candidates

  #   Returns
  #   -------
  #   float
  #   """
  #   accuracy = []
  #   if ref_annotation is None:
  #     ref = self.exist_annotation_formula
  #   else:
  #     ref = ref_annotation
  #   if pred_annotation is None:
  #     pred = self.formula
  #   else:
  #     pred = pred_annotation
  #   ref_keys = set(ref.keys())
  #   pred_keys = set(pred.keys())
  #   # select species that can be evaluated
  #   species_to_test = ref_keys.intersection(pred_keys)
  #   for one_k in species_to_test:
  #     if set(ref[one_k]).intersection(pred[one_k]):
  #       accuracy.append(True)
  #     else:
  #       accuracy.append(False)
  #   return np.mean(accuracy)


  # def getRecall(self,
  #              ref_annotation=None,
  #              pred_annotation=None,
  #              mean=True):
  #   """
  #   More precise version of 'accuracy',
  #   recall is the fraction of correct
  #   elements detected. Currently,
  #   it is calculated as the fraction of 
  #   correct chemical formula detected

  #   Parameters
  #   ----------
  #   ref_annotation: dict
  #       {species_id: [str-annotation, i.e., formula]}
  #       if None, get self.exist_annotation_formula
  #   pred_annotation: dict
  #       {species_id: [str-annotation, i.e., formula]}
  #       if None, get self.candidates  
  #   mean: bool
  #       If True, get model-level average
  #       If False, get value of each ID

  #   Returns
  #   -------
  #   float/dict{id: float}
  #       Depending on the 'mean' argument
  #   """
  #   recall = dict()
  #   if ref_annotation is None:
  #     ref = self.exist_annotation_formula
  #   else:
  #     ref = ref_annotation
  #   if pred_annotation is None:
  #     pred = self.formula
  #   else:
  #     pred = pred_annotation
  #   ref_keys = set(ref.keys())
  #   pred_keys = set(pred.keys())
  #   # select species that can be evaluated
  #   species_to_test = ref_keys.intersection(pred_keys)
  #   # go through each species
  #   for one_k in species_to_test:
  #     num_intersection = len(set(ref[one_k]).intersection(pred[one_k]))
  #     recall[one_k] = num_intersection / len(set(ref[one_k]))
  #   if mean:
  #     return np.mean([recall[val] for val in recall.keys()])
  #   else:
  #     return recall


  # def getPrecision(self,
  #                  ref_annotation=None,
  #                  pred_annotation=None,
  #                  mean=True):
  #   """
  #   A complementary term of 'recall'
  #   precision is the fraction of correct
  #   elements detected from all detected formulas. 
  
  #   Parameters
  #   ----------
  #   ref_annotation: dict
  #       {species_id: [str-annotation, i.e., formula]}
  #       if None, get self.exist_annotation_formula
  #   pred_annotation: dict
  #       {species_id: [str-annotation, i.e., formula]}
  #       if None, get self.candidates  
  #   mean: bool
  #       If True, get model-level average
  #       If False, get value of each ID

  #   Returns
  #   -------
  #   : float/dict{id: float}
  #       Depending on the 'mean' argument
  #   """
  #   precision = dict()
  #   if ref_annotation is None:
  #     ref = self.exist_annotation_formula
  #   else:
  #     ref = ref_annotation
  #   if pred_annotation is None:
  #     pred = self.formula
  #   else:
  #     pred = pred_annotation
  #   ref_keys = set(ref.keys())
  #   pred_keys = set(pred.keys())
  #   # select species that can be evaluated
  #   species_to_test = ref_keys.intersection(pred_keys)
  #   # go through each species
  #   for one_k in species_to_test:
  #     num_intersection = len(set(ref[one_k]).intersection(pred[one_k]))
  #     precision[one_k] = num_intersection / len(set(pred[one_k]))
  #   if mean:
  #     return np.mean([precision[val] for val in precision.keys()])
  #   else:
  #     return precision

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
    str
    """
    species_name = self.names[inp_id]
    if len(species_name) > 0:
      res_name = species_name
    else:
      res_name = inp_id
    return res_name

  def evaluatePredictedSpeciesAnnotation(self,
                                         pred_result=None,
                                         id_list=None,
                                         fitted_model=SPECIES_RF):
    """
    Evaluate the quality of annotation;
    for each individual species.
    
    Parameters
    ---------
    pred_result: dict 
        Result of prediction of one species.
        {'name_used':str, 'chebi':[str-chebi],
         'match_score': [chebi_tuples],
         'formula': [str-formula]}
    id_list: str-list
        List of species IDs to evaluate

    Returns
    -------
    dict {species_id: level-of-species-prediction-being-correct}
        Information of whether confident or not
    """
    name_length = len(pred_result)
    num_candidates = len(pred_result[cn.CHEBI])
    match_score = np.mean([val[1] for val in pred_result[cn.MATCH_SCORE]])
    return fitted_model.predict([[name_length, num_candidates, match_score]])[0]


  def updateSpeciesWithRecommendation(self, inp_recom):
    """
    Update species_annotation class using
    Recommendation namedtuple.
  
    self.candidates is a sorted list of tuples,
    (chebi_id: match_score)
    self.formula is a unsorted list of unique formulas
  
    Parameters
    ----------
    inp_recom: Recommendation
       A namedtuple. Created by recom.getSpeciesAnnotation
  
    Returns
    -------
    None
    """
    self.candidates.update({inp_recom.id: inp_recom.candidates})
    formulas2update = list(set([cn.REF_CHEBI2FORMULA[val[0]] for val in inp_recom.candidates]))
    self.formula.update({inp_recom.id: formulas2update})
    return None


  def updateSpeciesWithDict(self, inp_dict):
    """
    A direct way of updating species annotations,
    using ChEBI terms.
    As match scores are given
    when exact matches are found, 
    match scores were given as 1.0. 
  
    Parameters
    ----------
    inp_dict: dict
        {species_id: [chebi terms]}
  
    Returns
    -------
    None
    """
    info2upd_candidates = {k:[(val, 1.0) for val in inp_dict[k]] for k in inp_dict.keys()}
    info2upd_formula = {k:[cn.REF_CHEBI2FORMULA[chebi] for chebi in inp_dict[k]] for k in inp_dict.keys()}
    self.candidates.update(info2upd_candidates)
    self.formula.update(info2upd_formula)














