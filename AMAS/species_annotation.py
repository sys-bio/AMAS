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

import collections
import compress_pickle
import editdistance
import itertools
import libsbml
import numpy as np
import operator
import os
import pandas as pd
import pickle
import re


with open(os.path.join(cn.REF_DIR, 'chebi_low_synonyms_comp.lzma'), 'rb') as f:
  CHEBI_LOW_SYNONYMS = compress_pickle.load(f)
CHARCOUNT_COMB_DF = compress_pickle.load(os.path.join(cn.REF_DIR, 'charcount_df_scaled.lzma'),
                                         compression="lzma")
CHARCOUNT_DF = CHARCOUNT_COMB_DF.iloc[:, :-2]
CHEBI_DF = CHARCOUNT_COMB_DF.iloc[:, -2:]
# A trained random forest model 
SPECIES_RF = compress_pickle.load(os.path.join(cn.REF_DIR, 'species_rf_fitted.lzma'))


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
      self.exist_annotation = tools.extractExistingSpeciesAnnotation(self.model)
      exist_annotation_formula_raw = {k:tools.transformCHEBIToFormula(self.exist_annotation[k], cn.REF_CHEBI2FORMULA) \
                                      for k in self.exist_annotation.keys()}
      self.exist_annotation_formula = {val:exist_annotation_formula_raw[val] for val in exist_annotation_formula_raw.keys() \
                                       if exist_annotation_formula_raw[val]}

    # inp_tuple: ({species_id:species_name}, {species_id: [CHEBI annotations]})
    elif inp_tuple is not None:
      self.model = None
      self.names = inp_tuple[0]
      self.exist_annotation = inp_tuple[1]
      exist_annotation_formula_raw = {k:tools.transformCHEBIToFormula(inp_tuple[1][k], cn.REF_CHEBI2FORMULA) \
                                      for k in inp_tuple[1].keys()}
      self.exist_annotation_formula = {val:exist_annotation_formula_raw[val] for val in exist_annotation_formula_raw.keys() \
                                       if exist_annotation_formula_raw[val]}
    else:
      self.model = None
      self.names = None
      self.exist_annotation = None
      self.exist_annotation_formula = None
    # Below are predicted annotations in dictionary format
    # Once created, each will be {species_ID: float/str-list}
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
        {'name_used': str,
         'chebi': [list-ChEBI],
         'match_score': [(ChEBI, float)],
         'formula': [list-formula]} 
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
                                    for val in CHEBI_LOW_SYNONYMS[one_chebi]]), cn.ROUND_DIGITS)) \
                 for one_chebi in min_min_chebis] 
    res_tuple.sort(key=operator.itemgetter(1), reverse=True)
    #  CHEBI part is added, because we want a sorted list after computing res_tuple
    one_result[cn.NAME_USED] = inp_str
    one_result[cn.CHEBI] = [val[0] for val in res_tuple]
    one_result[cn.MATCH_SCORE] = res_tuple
    min_min_formula = list(set([cn.REF_CHEBI2FORMULA[val] for val in min_min_chebis]))
    one_result[cn.FORMULA] = min_min_formula
    return one_result

  # Methods to use Cosine Similarity
  def getCountOfIndividualCharacters(self, inp_str):
    """
    Get a list of characters
    between a-z and 0-9. 
  
    Parameters
    ----------
    inp_str: str
  
    Returns
    -------
    : list
    """
    return collections.Counter(itertools.chain(*re.findall('[a-z0-9]+', inp_str.lower())))

  def prepareCounterQuery(self,
                          specs,
                          ref_cols,
                          use_id=True):
    """
    Prepare a query vector, which will be used  
    as a vector for predictor variables.
    Input will be a list of
    IDs using which names_used will be determined. 
    In addition, querys will also be scaled
    by the length of each vector. 
  
    There is 'use_id' option, so
    if False, directly use the string
    instead of searching for used_name. 
  
    Parameters
    ----------
    list-str: specs
        IDs of species
    ref_cols: list-str
        Column names to use
    use_id: bool
        If False, directly use the string
        If True, use getNameToUse
      
    Returns
    -------
    : pandas.DataFrame
    : dict
    """
    name_used = dict()
    query_mat = pd.DataFrame(0, index=ref_cols, columns=specs)
    for one_spec in specs:
      if use_id:
        name2use = self.getNameToUse(one_spec)
        # characters are lowered in getCountOfIndividualCharacters()
        char_counts = self.getCountOfIndividualCharacters(name2use)
        name_used[one_spec] = name2use
      else:
        name2use = one_spec
        # characters are lowered in getCountOfIndividualCharacters()
        char_counts = self.getCountOfIndividualCharacters(name2use)
        name_used[one_spec] = name2use
      for one_char in char_counts:
        query_mat.loc[one_char, one_spec] = char_counts[one_char] 
    # Now, scale it using the vector distance
    div_row = query_mat.apply(lambda col : np.sqrt(np.sum([val**2 for val in col])), axis = 0)
    norm_query = query_mat.divide(div_row, axis=1)
    return norm_query, name_used


  def predictAnnotationByCosineSimilarity(self,
                                          inp_strs=None,
                                          inp_ids=None,
                                          ref_df=CHARCOUNT_DF,
                                          chebi_df=CHEBI_DF):
    """
    Predict annotation by taking cosine distance 
    of character count vectors.
  
    Parameters
    ----------
    inp_strs: list-str
        Strings that will directly used
        for prediction
    inp_ids: list-str
        IDs with which name2use will be
        determined
    ref_df: DataFrame
        Reference database
    chebi_df: DataFrame
        ChEBI information sharing the index with ref_df    

    Returnsa
    -------
    : dict/None
        {'name_used': str,
         'chebi': [list-ChEBI],
         'match_score': [(ChEBI, float)],
         'formula': [list-formula]} 
      if no name/ID is given, return None
    """
    if inp_ids:
      one_query, name_used = self.prepareCounterQuery(specs=inp_ids,
                                                      ref_cols=ref_df.columns,
                                                      use_id=True)
    elif inp_strs:
      one_query, name_used = self.prepareCounterQuery(specs=inp_strs,
                                                      ref_cols=ref_df.columns,
                                                      use_id=False)  
    else:
      return None
    multi_mat = ref_df.dot(one_query)
    max_val = multi_mat.max()
    result = dict()
    for one_spec in one_query.columns:
      one_res = dict()
      one_res[cn.NAME_USED] = name_used[one_spec]
      cand_index = multi_mat[abs(multi_mat[one_spec]-max_val[one_spec])<cn.TOLERANCE].index
      one_res[cn.CHEBI] = list(set(chebi_df.loc[cand_index, cn.CHEBI]))
      one_res[cn.MATCH_SCORE] = [(val, np.round(max_val[one_spec], cn.ROUND_DIGITS)) \
                                 for val in one_res[cn.CHEBI]]
      one_res[cn.FORMULA] = list(set([cn.REF_CHEBI2FORMULA[val] for val in one_res[cn.CHEBI] \
                                      if val in cn.REF_CHEBI2FORMULA.keys()])) 
      result[one_spec] = one_res
    return result

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
                                         fitted_model=SPECIES_RF):
    """
    Predict the probability of 
    the candidate set including the 'true'
    annotation. 
    
    Parameters
    ---------
    pred_result: dict 
        Result of prediction of one species.
        {'name_used':str, 'chebi':[str-chebi],
         'match_score': [chebi_tuples],
         'formula': [str-formula]}


    Returns
    -------
    : float
    """
    name_length = len(pred_result[cn.NAME_USED])
    num_candidates = len(pred_result[cn.CHEBI])
    match_score = np.mean([val[1] for val in pred_result[cn.MATCH_SCORE]])
    num_formulas = len(pred_result[cn.FORMULA])
    proba_correct = fitted_model.predict_proba([[name_length,
                                                 num_candidates,
                                                 match_score,
                                                 num_formulas]])[0][1]
    return proba_correct


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
    formulas2update = list(set([cn.REF_CHEBI2FORMULA[val[0]] \
                                for val in inp_recom.candidates \
                                if val[0] in cn.REF_CHEBI2FORMULA.keys()]))
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














