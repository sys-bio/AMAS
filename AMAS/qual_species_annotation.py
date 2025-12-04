# gene_annotation.py
"""
<Annotation for Qualitative Species>
qual_species_annotation creates and predicts
annotation of qualitative species.
"""

from AMAS import constants as cn
from AMAS import tools

import collections
import itertools
import libsbml
import numpy as np
import pandas as pd
import operator
import re
import warnings


class QualSpeciesAnnotation(object):

  def __init__(self, libsbml_fpath=None,
               inp_tuple=None):

    """
    Parameters
    ----------
    libsbml_fpath: str
        File path of an SBMl (.xml/.sbml) model
        Assuming it use qual SBML extension.

    inp_tuple: tuple 
        Tuple of model information,
        first element (index 0) is information on 
        qualitative species names,
        second element (index 1) is existing 
        qualitative species information.
        ({species_id: species_display_name},
         {species_id: qualitative species terms})
    """
    # self.exist_annotation stores existing annotations in the model
    ## depends on the qualifier, by default it is NCBI gene
    # If none exists, set None

    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      if document.getModel().getPlugin("qual") is not None:
        self.model = document.getModel().getPlugin("qual")
        self.names = {val.getId():val.getName() for val in self.model.getListOfQualitativeSpecies()}
      else:
        self.model = None
        self.names = None
      self.exist_annotation = tools.extractExistingQualitativeSpeciesAnnotation(self.model, qualifier=cn.NCBI_GENE, description=True)
      
    # inp_tuple: ({gene_id:gene_name}, {gene_id: [NCBI gene annotations]})
    elif inp_tuple is not None:
      self.model = None
      self.names = inp_tuple[0]
      self.exist_annotation = inp_tuple[1]
      
    else:
      self.model = None
      self.names = None
      self.exist_annotation = None
      
    # Below are predicted annotations in dictionary format
    # Once created, each will be {gene_ID: float/str-list}
    self.candidates = dict()

  def getCScores(self,
                 tax_id,
                 inp_strs,
                 mssc,
                 cutoff):
    """
    Compute the cScores
    of query strings with
    all possible NCBI gene terms. 
    A sorted list of tuples 
    (ncbigene:XXXXX, cScore)
    will be returned.
    Only unique strings 
    will be calculated to avoid
    cases such as {'a': 'a',
                   'a': 'b'}.
  
    Parameters
    ----------
    tax_id: str
        Taxonomy ID of the species
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
  
    Returns
    -------
    :dict
        {one_str: [(XXXXX, 1.0), ...]}
    """
    # Filter the reference dataframe based on the taxonomy
    ref_df = tools.get_gene_char_count_df(tax_id).iloc[:, :-3]
    info_df = tools.get_gene_char_count_df(tax_id).iloc[:, -2:]

    # Prepare the query matrix
    unq_strs = list(set(inp_strs))
    one_query, name_used = self.prepareCounterQuery(species=unq_strs,
                                                    ref_cols=ref_df.columns,
                                                    use_id=False) 
    multi_mat = ref_df.dot(one_query)
    # updated code to avoid repeated prediction
    cscores = dict()
    multi_mat[cn.NCBI_GENE] = info_df[cn.NCBI_GENE]
    for species in inp_strs:    
      # Get max-value of each term
      g_res = multi_mat.loc[:,[cn.NCBI_GENE, species]].groupby([cn.NCBI_GENE]).max()[species]
      species_cscore = tools.applyMSSC(pred=zip(g_res.index, g_res),
                                    mssc=mssc,
                                    cutoff=cutoff)
      cscores[species] = species_cscore
    return cscores

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
                          species,
                          ref_cols,
                          use_id=True):
    """
    Prepare a query vector, which will be used  
    as a vector for predictor variables.
    Input will be a list of
    IDs using which names_used will be determined. 
    In addition, querys will also be scaled
    by the length of each vector. 
  
    Parameters
    ----------
    species: list-str
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
    query_mat = pd.DataFrame(0, index=ref_cols, columns=species)
    for one_species in species:
      if use_id:
        name2use = self.getNameToUse(one_species)
      else:
        name2use = one_species

      # Remove everything after the first underscore
      name2use = name2use.split('_', 1)[0]

      # characters are lowered in getCountOfIndividualCharacters()
      char_counts = self.getCountOfIndividualCharacters(name2use)
      name_used[one_species] = name2use
      for one_char in char_counts:
        query_mat.loc[one_char, one_species] = char_counts[one_char] 
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

  def updateQualSpeciesWithRecommendation(self, inp_recom):
    """
    Update qual_species_annotation class using
    Recommendation namedtuple.
  
    self.candidates is a sorted list of tuples,
    (qual_species: match_score)
  
    Parameters
    ----------
    inp_recom: cn.Recommendation
       A namedtuple. Created by recom.getQualSpeciesRecommendation
  
    Returns
    -------
    None
    """
    self.candidates.update({inp_recom.id: inp_recom.candidates})
    return None

  def updateQualSpeciesWithDict(self, inp_dict):
    """
    A direct way of updating qual species annotations,
    using qual species terms.
    As match scores are given
    when exact matches are found, 
    match scores were given as 1.0. 
  
    Parameters
    ----------
    inp_dict: dict
        {qual_species_id: [ncbigene terms]}
  
    Returns
    -------
    None
    """
    info2upd_candidates = {k:[(val, 1.0) for val in inp_dict[k]] for k in inp_dict.keys()}
    self.candidates.update(info2upd_candidates)

################### ES score (not used)###################
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

  def prepare_synonyms_dict(info_df):
      """
      Prepares a dictionary from the info_df DataFrame, grouping synonyms by GeneID.

      Parameters
      ----------
      info_df: pd.DataFrame
          DataFrame containing GeneID and synonym columns (both lowercased).

      Returns
      -------
      dict
          Dictionary with GeneID as keys and lists of synonyms as values.
      """
      # Ensure that GeneID and Synonyms are lowercase
      info_df = info_df.applymap(str.lower)
      
      # Group synonyms by GeneID
      synonyms_dict = info_df.groupby(cn.NCBI_GENE)['synonym'].apply(list).to_dict()
      return synonyms_dict

  def getEScores(self,
                 tax_id,
                 inp_strs,
                 mssc,
                 cutoff):
    """
    Compute the eScores
    of a list of query strings with
    all possible NCBI gene ids. 
    A sorted list of tuples 
    (XXXXX, eScore)
    will be returned.
    Only unique strings
    will be calculated. 
  
    Parameters
    ----------
    tax_id: str
        Taxonomy ID of the species
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
        {one_str: [(XXXXX, 1.0), ...]}
    """
    # Filter the reference dataframe based on the taxonomy
    ref_df = tools.get_gene_char_count_df(tax_id)
    synonyms_dict = self.prepare_synonyms_dict(ref_df)

    unq_strs = list(set(inp_strs))
    escores = dict()
    for spec in unq_strs:
      spec_escore = [(one_k, np.max([self.getOneEScore(spec.lower(), val) \
                               for val in synonyms_dict[one_k]])) \
                     for one_k in synonyms_dict.keys()]
      mssc_escore = tools.applyMSSC(pred=spec_escore,
                                    mssc=mssc,
                                    cutoff=cutoff)
      mssc_escore.sort(key=operator.itemgetter(1), reverse=True)
      escores[spec] = mssc_escore
    return escores


