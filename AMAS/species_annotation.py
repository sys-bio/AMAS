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
import warnings


with open(os.path.join(cn.REF_DIR, 'chebi_low_synonyms_comp.lzma'), 'rb') as f:
  CHEBI_LOW_SYNONYMS = compress_pickle.load(f)
CHARCOUNT_COMB_DF = compress_pickle.load(os.path.join(cn.REF_DIR, 'charcount_df_scaled.lzma'),
                                         compression="lzma")
CHARCOUNT_DF = CHARCOUNT_COMB_DF.iloc[:, :-2]
CHEBI_DF = CHARCOUNT_COMB_DF.iloc[:, -2:]

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

  def getCScores(self,
                 inp_strs,
                 mssc,
                 cutoff,
                 ref_df=CHARCOUNT_DF,
                 chebi_df=CHEBI_DF):
    """
    Compute the eScores
    of query strings with
    all possible ChEBI terms. 
    A sorted list of tuples 
    (CHEBI:XXXXX, eScore)
    will be returned.
    Only unique strings 
    will be calculated to avoid
    cases such as {'a': 'a',
                   'a': 'b'}.
  
    Parameters
    ----------
    inp_strs: list-str
        List of strings
    mssc: match score selection criteria
        'top' will recommend candidates with
        the highest match score above cutoff
        'above' will recommend all candidates with
        match scores above cutoff
    cutoff: float
        Cutoff value; only candidates with match score
        at or above the cutoff will be recommended.
    ref_df: DataFrame
        Reference database
    chebi_df: DataFrame
        ChEBI information sharing the index with ref_df 
  
    Returns
    -------
    :dict
        {one_str: [(CHEBI:XXXXX, 1.0), ...]}
    """
    unq_strs = list(set(inp_strs))
    one_query, name_used = self.prepareCounterQuery(specs=unq_strs,
                                                    ref_cols=ref_df.columns,
                                                    use_id=False) 
    multi_mat = ref_df.dot(one_query)
    # updated code to avoid repeated prediction
    cscores = dict()
    multi_mat[cn.CHEBI] = chebi_df[cn.CHEBI]
    for spec in inp_strs:
      # Get max-value of each chebi term
      g_res = multi_mat.loc[:,[cn.CHEBI, spec]].groupby([cn.CHEBI]).max()[spec]
      spec_cscore = tools.applyMSSC(pred=zip(g_res.index, g_res),
                                    mssc=mssc,
                                    cutoff=cutoff)
      cscores[spec] = spec_cscore
    # cscores = dict()
    # for spec in unq_strs:
    #   spec_cscore = tools.applyMSSC(pred=zip(chebi_df[cn.CHEBI], multi_mat[spec]),
    #                                 mssc=mssc,
    #                                 cutoff=cutoff)
    #   spec_cscore.sort(key=operator.itemgetter(1), reverse=True)
    #   cscores[spec] = spec_cscore
    return cscores

  def getOneEScore(self, one_s, two_s):
    """
    Compute the eScore 
    of a pair of two strings using
    the formula below:
    1.0 - (editdistance(one_s, two_s) / max(len(one_s, two_s)))
  
    Values should be between 0.0 and 1.0.
  
    Parameters
    ----------
    one_s: str
    two_s: str
  
    Returns
    -------
    : float (0.0-1.0)
    """
    edist = editdistance.eval(one_s, two_s)/ max(len(one_s), len(two_s))
    escore = 1.0 - edist
    return escore

  def getEScores(self,
                 inp_strs,
                 mssc,
                 cutoff):
    """
    Compute the eScores
    of a list of query strings with
    all possible ChEBI terms. 
    A sorted list of tuples 
    (CHEBI:XXXXX, eScore)
    will be returned.
    Only unique strings
    will be calculated. 
  
    Parameters
    ----------
    inp_strs: str
        List of strings
    mssc: match score selection criteria
        'top' will recommend candidates with
        the highest match score above cutoff
        'above' will recommend all candidates with
        match scores above cutoff
    cutoff: float
        Cutoff value; only candidates with match score
        at or above the cutoff will be recommended.
  
    Returns
    -------
    :dict
        {one_str: [(CHEBI:XXXXX, 1.0), ...]}
    """
    unq_strs = list(set(inp_strs))
    escores = dict()
    for spec in unq_strs:
      spec_escore = [(one_k, np.max([self.getOneEScore(spec.lower(), val) \
                               for val in CHEBI_LOW_SYNONYMS[one_k]])) \
                     for one_k in CHEBI_LOW_SYNONYMS.keys() \
                     if one_k in cn.REF_CHEBI2FORMULA.keys()]
      mssc_escore = tools.applyMSSC(pred=spec_escore,
                                    mssc=mssc,
                                    cutoff=cutoff)
      mssc_escore.sort(key=operator.itemgetter(1), reverse=True)
      escores[spec] = mssc_escore
    return escores

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

  def updateSpeciesWithRecommendation(self, inp_recom):
    """
    Update species_annotation class using
    Recommendation namedtuple.
  
    self.candidates is a sorted list of tuples,
    (chebi_id: match_score)
    self.formula is a unsorted list of unique formulas
  
    Parameters
    ----------
    inp_recom: cn.Recommendation
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
    info2upd_formula = {k:[cn.REF_CHEBI2FORMULA[chebi] \
                       for chebi in inp_dict[k]] for k in inp_dict.keys()}
    self.candidates.update(info2upd_candidates)
    self.formula.update(info2upd_formula)




