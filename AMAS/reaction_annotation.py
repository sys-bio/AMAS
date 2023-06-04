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
import libsbml
import numpy as np
import operator
import os
import pandas as pd
import pickle


with open(os.path.join(cn.REF_DIR, 'data2ref_mat.lzma'), 'rb') as handle:
  REF_DAT = compress_pickle.load(handle)

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
    # get dictionary of reaction ID: species component
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
    # ref_rowsum = ref_mat.sum(1)
    query_colsum = query_df.sum(0)
    # divided by the number of elements of the QUERY
    div_mat = multi_mat.divide(query_colsum, axis=1)
    rscores = dict()
    for reac in reacs:
      reac_rscore = tools.applyMSSC(pred=list(zip(div_mat.index, div_mat[reac])),
                                    mssc=mssc,
                                    cutoff=cutoff)
      reac_rscore.sort(key=operator.itemgetter(1), reverse=True)
      rscores[reac] = reac_rscore    
    return rscores

  # # remove and replace
  # def predictAnnotation(self,
  #                       inp_spec_dict,
  #                       inp_reac_list=None,
  #                       inp_ref_mat=REF_MAT,
  #                       update=False):
  #   """
  #   Predict 1) reaction annotation candidates 
  #   and 2) match score of them
  #   using species dict (argument) etc.
  #   Also updates several class attributes as a result. 
  
  #   Parameters
  #   ----------
  #   dict: inp_spec_dict
  #       Dictionoary, {species id: formula(str-list)}
  #   inp_reac_list: str-list
  #       IDs of reactions to predict. If default, will do all reactions
  #   inp_ref_mat: pd.DataFrame
  #       Reference matrix
  #   update: bool
  #       Whether to save results as class attrributes or just return them.
  #       If True, only match_score is returned (and values are updated).
  #       If False, all relevant information is returned

  #   Returns
  #   -------
  #   : dict
  #       {'candidates': {reactionID: [candidates in RHEA]},
  #        'match_score': {reactionID: [(Rhea ID, match score: float between 0.0-1.0),]}
  #        'query_df': query_df}
  #   """
  #   # get libsbml.reaction and their IDs
  #   if inp_reac_list is not None:
  #     reaction_ids = inp_reac_list
  #   else:
  #     reaction_ids = list(self.reaction_components.keys())
  #   # get dictionary of reaction ID: species component
  #   r2pred_spec_formulas = dict()
  #   for one_rid in reaction_ids:
  #     r2pred_spec_formulas[one_rid] = {one_spec:inp_spec_dict[one_spec] \
  #                                      for one_spec in self.reaction_components[one_rid]}
  #   # prepare query df for prediction
  #   query_df = pd.DataFrame(0, 
  #                           index=inp_ref_mat.columns,
  #                           columns=reaction_ids)
  #   for one_rid in reaction_ids:
  #     one_set_species = r2pred_spec_formulas[one_rid]
  #     # for each species element of the select reaction
  #     for one_spec_key in one_set_species.keys():
  #       one_spec = one_set_species[one_spec_key]
  #       # For each one_rid, set the values 1.0
  #       query_df.loc[[val for val in one_spec if val in query_df.index], one_rid] = 1
  #   multi_mat = inp_ref_mat.dot(query_df)
  #   maxes = multi_mat.max()
  #   #
  #   # Collect candidates and calculate confidence score
  #   pred_cands = dict()
  #   pred_match_score = dict()
  #   for one_rid in maxes.index:
  #     one_multi = multi_mat.loc[:,one_rid]
  #     candidates = one_multi[one_multi==maxes[one_rid]].index
  #     # cand_data; (number of element matches, candidates)
  #     pred_cands[one_rid] = candidates
  #     # Now, match_scpre (calculated per each candidate) => replaced as a tuple
  #     match_score_per_cand = []
  #     for one_cand in candidates:
  #       if one_cand in cn.REF_RHEA2MASTER.keys():
  #         num_matches = maxes[one_rid]
  #         num_maxpos_matches = len(inp_ref_mat.loc[one_cand, :].to_numpy().nonzero()[0])
  #         match_score_per_cand.append((one_cand, np.round(num_matches/num_maxpos_matches, cn.ROUND_DIGITS)))
  #     match_score_per_cand.sort(key=operator.itemgetter(1), reverse=True)
  #     pred_match_score[one_rid] = match_score_per_cand
  #   if update:
  #     self.candidates = pred_match_score
  #     self.query_df = query_df
  #   #
  #   return {'candidates': pred_cands,
  #           'match_score': pred_match_score,
  #           'query_df': query_df}


  # # Develop a method to evaluate results using fitted model
  # def evaluatePredictedReactionAnnotation(self, pred_result,
  #                                         fitted_model=REACTION_RF):
  #   """
  #   Evaluate the quality of annotation;
  #   for each individual species.
  #   pred_result includes predicted results,
  #   a result of self.predictAnnotation().
  #   All informaton needed for prediction 
  #   is supposed to come from inp_dict, 
  #   not the stored reaction class itself.
  
  #   Parameters
  #   ---------
  #   pred_result: dict
  #       {'candidates': {reactionID: [candidates in RHEA]},
  #        'match_score': {reactionID: [(Rhea ID, match score: float between 0.0-1.0),]}
  #        'query_df': query_df}

  #   Returns
  #   -------  
  #   dict {reaction_id: probability-of-prediction-including-correct-value}
  #       Information of how algorithm is confident about the result
  #   """
  #   candidates_dict = pred_result[cn.CANDIDATES]
  #   match_score_dict = pred_result[cn.MATCH_SCORE]
  #   mean_rheas_num_dict = {one_k: np.mean([self.getRheaElementNum(val) \
  #                                          for val in candidates_dict[one_k]]) \
  #                          for one_k in candidates_dict.keys()}
  #   num_reac_comp_dict = {one_k: len(self.reaction_components[one_k]) \
  #                         for one_k in candidates_dict.keys()}
  #   num_candidates = {one_k: len(candidates_dict[one_k]) \
  #                         for one_k in candidates_dict.keys()}
  #   mean_match_scores = {one_k: np.mean([val[1] for val in match_score_dict[one_k]]) \
  #                         for one_k in candidates_dict.keys()}
  #   df2pred = pd.DataFrame([mean_rheas_num_dict,
  #                           num_reac_comp_dict,
  #                           num_candidates,
  #                           mean_match_scores]).T
  #   cred_pred = fitted_model.predict_proba(df2pred)
  #   prob_1_dict = {val: cred_pred[idx][1] for idx, val in enumerate(df2pred.index)}
  #   return prob_1_dict

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


