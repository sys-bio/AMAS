# reaction_annotation.py
"""
<Annotation for Reaction>
reaction_annotation creates and predicts
annotation of libsbml reactions,
using existing  species annotations or,
annotations predicted via species_annotation.
"""

from AMAS import constants as cn
from AMAS import tools

import compress_pickle
import itertools
import libsbml
import numpy as np
import operator
import os
import pandas as pd
import pickle


with open(os.path.join(cn.REF_DIR, 'data2ref_mat.lzma'), 'rb') as handle:
  REF_DAT = compress_pickle.load(handle)
# might need to be deleted after trying Jaccard Index
REF_NONZERO_COLS = compress_pickle.load(os.path.join(cn.REF_DIR, 'ref_nonzero_cols.lzma'),
                        compression="lzma",
                        set_default_extension=False)

# first of list is list of columns
cols = REF_DAT[0]
# second, list of indices
inds = REF_DAT[1]
# third, list of index (column, [non-zero rows])
ref_mat_pairs = REF_DAT[2]
REF_MAT = pd.DataFrame(0, index=inds, columns=cols)
for val in ref_mat_pairs:
  REF_MAT.iloc[val[1], val[0]] = 1

REACTION_RF = compress_pickle.load(os.path.join(cn.REF_DIR, 'reactions_rf_fitted.lzma'))

class ReactionAnnotation(object):

  def __init__(self, libsbml_fpath=None, 
               inp_tuple=None):
  
    """
    Parameters
    ----------
    libsbml_fpath: str
        File path of an SBMl model (.xml)
    inp_tuple: tuple
        ({reaction_id: [unique components (that is, species) of that reaction]},
         {reaction_id: Rhea terms})
    """
    # self.exist_annotation stores 
    # existing KEGG Reaction or Rhea annotations in the model.
    # If none exists, set None.
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      self.exist_annotation = tools.extractExistingReactionAnnotation(inp_model=self.model)
      self.reaction_components = {val.getId():list(set([k.species for k in val.getListOfReactants()]+\
                                                       [k.species for k in val.getListOfProducts()])) \
                                  for val in self.model.getListOfReactions()}
    elif inp_tuple is not None:
      self.model = None
      self.reaction_components = inp_tuple[0]
      self.exist_annotation = inp_tuple[1]
    else:
      self.model = None
      self.reaction_components = None
      self.exist_annotation = None
    # Attributes after prediction
    self.candidates = None
    self.query_df = None

  def getReactionComponents(self,
                            inp_reaction):
    """
    Get component of reactions in species IDs
    (both reactants and products)
    of a reaction. 

    Parameters
    ----------
    reaction_id: str/libsbml.Reaction


    Returns
    -------
    None/str-list (list of species IDs)
    """
    if isinstance(inp_reaction, libsbml.Reaction):
      one_reaction = inp_reaction
    elif isinstance(inp_reaction, str):
      one_reaction = self.model.getReaction(inp_reaction)
    else:
      return None
    reactants = [val.species for val in one_reaction.getListOfReactants()]
    products = [val.species for val in one_reaction.getListOfProducts()]
    r_components = list(set(reactants + products))
    return r_components
  
  def getRScores(self,
                 spec_dict,
                 reacs,
                 mssc,
                 cutoff,
                 ref_mat=REF_MAT):
    """
    Get a sorted list of
    Rhea-rScore tuples.
    [(RHEA:XXXXX, 1.0), etc.]
  
    Parameters
    ----------
    dict: inp_spec_dict
        Dictionoary, {species id: formula(str-list)}
    reacs: str-list
        IDs of reactions to predict annotatinos.
    mssc: match score selection criteria
        'top' will recommend candidates with
        the highest match score above cutoff
        'above' will recommend all candidates with
        match scores above cutoff
    cutoff: float
        Cutoff value; only candidates with match score
        at or above the cutoff will be recommended.
    ref_mat: pd.DataFrame
        Reference matrix
      
    Returns
    -------
    :dict
        {one_str: [(Rhea:XXXXX, 1.0), ...]}
    """

    # def get_jaccard(inp1, inp2):
    #   """
    #   Get jaccard score of two lists; 
  
    #   Parameters
    #   ----------
    #   inp1: list-str
    #   inp2: list-str
  
    #   Returns
    #   -------
    #   : float
    #   """
    #   inters = set(inp1).intersection(inp2)
    #   uni = set(inp1).union(inp2)
    #   return len(inters) / len(uni)

    # r2comb_spec_formulas = dict()
    # for one_rid in reacs:
    #   r2comb_spec_formulas[one_rid] = list(set(itertools.chain(*[spec_dict[spec] \
    #                                    for spec in self.reaction_components[one_rid]])))
    # j_rscores = dict()
    # for rid in r2comb_spec_formulas.keys():
    #   score = REF_NONZERO_COLS.apply(lambda x: get_jaccard(x, r2comb_spec_formulas[rid]))
    #   reac_rscore = tools.applyMSSC(pred=zip(score.index, score),
    #                                 mssc=mssc,
    #                                 cutoff=cutoff)
    #   reac_rscore.sort(key=operator.itemgetter(1), reverse=True)
    #   j_rscores[rid] = reac_rscore
    # return j_rscores

    # BELOW IS THE ORIGINAL MINI-MAX VERSION
    # Get dictionary of reaction ID: species component
    r2pred_spec_formulas = dict()
    for one_rid in reacs:
      r2pred_spec_formulas[one_rid] = {spec:spec_dict[spec] \
                                       for spec in self.reaction_components[one_rid]}
    # prepare query df for prediction
    query_df = pd.DataFrame(0, 
                            index=ref_mat.columns,
                            columns=reacs)
    for one_rid in reacs:
      one_set_species = r2pred_spec_formulas[one_rid]
      # for each species element of the select reaction
      for spec_key in one_set_species.keys():
        one_spec = one_set_species[spec_key]
        # For each one_rid, set the values 1.0
        query_df.loc[[val for val in one_spec if val in query_df.index], one_rid] = 1
    multi_mat = ref_mat.dot(query_df)
    # new minimax of reference value
    max_multi_mat = np.max(multi_mat)
    query_colsum = pd.Series(0, index=max_multi_mat.index)
    for idx in query_colsum.index:
      query_colsum.at[idx] = np.min(np.sum(ref_mat.loc[multi_mat[multi_mat[idx]==max_multi_mat[idx]][idx].index,:],1))
    # divided 
    div_mat = multi_mat.divide(query_colsum, axis=1)
    rscores = dict()
    for reac in reacs:
      reac_rscore = tools.applyMSSC(pred=zip(div_mat.index, div_mat[reac]),
                                    mssc=mssc,
                                    cutoff=cutoff)
      reac_rscore.sort(key=operator.itemgetter(1), reverse=True)
      rscores[reac] = reac_rscore    
    return rscores

  def getRheaElementNum(self,
                        inp_rhea,
                        inp_df=REF_MAT):
    """
    Get Number of elements of
    the given rhea term.
    
    Parameters
    ----------
    inp_rhea: str
    
    Returns
    -------
    : int
    """
    return len(inp_df.loc[inp_rhea, :].to_numpy().nonzero()[0])


