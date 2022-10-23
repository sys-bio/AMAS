"""Annotation for Reactions."""

from AMAS import constants as cn
from AMAS import tools

import compress_pickle
import libsbml
import numpy as np
import operator
import os
import pandas as pd
import pickle


with open(os.path.join(cn.REF_DIR, 'rhea_all2bi_comp.lzma'), 'rb') as f:
  ref_rhea2bi = compress_pickle.load(f)
with open(os.path.join(cn.REF_DIR, 'kegg2rhea_bi_comp.lzma'), 'rb') as handle:
  ref_kegg2rhea_bi = compress_pickle.load(handle)
with open(os.path.join(cn.REF_DIR, 'rhea2chebi_comp.lzma'), 'rb') as f:
  ref_rhea_to_chebi = compress_pickle.load(f)
# with open(os.path.join(cn.REF_DIR, 'chebi_shortened_formula_comp.lzma'), 'rb') as f:
#   ref_shortened_chebi_to_formula = compress_pickle.load(f)
with open(os.path.join(cn.REF_DIR, 'dat_ref_mat_comp.lzma'), 'rb') as handle:
  ref_dat = compress_pickle.load(handle)

# first of list is list of columns
cols = ref_dat[0]
# second, list of indices
inds = ref_dat[1]
# third, list of index (column, [non-zero rows])
ref_mat_pairs = ref_dat[2]
ref_mat = pd.DataFrame(0, index=inds, columns=cols)
for val in ref_mat_pairs:
  ref_mat.iloc[val[1], val[0]] = 1

reaction_rf = pickle.load(open(os.path.join(cn.REF_DIR, 'reaction_randomforestcv.sav'), 'rb'))

class ReactionAnnotation(object):

  def __init__(self, libsbml_fpath=None, 
               inp_tuple=None):
    # self.exist_annotation stores 
    # existing KEGG Reaction or Rhea annotations in the model.
    # If none exists, set None.
    if libsbml_fpath is not None:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
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
      exist_annotation = reac_dict_rhea
      for one_id in reac_dict_kegg.keys():
        if one_id in exist_annotation.keys():
          exist_annotation[one_id] = list(set(exist_annotation[one_id] + reac_dict_kegg[one_id]))
        else:
          exist_annotation[one_id] = list(set(reac_dict_kegg[one_id]))
      self.exist_annotation = exist_annotation
      self.reaction_components = {val.getId():list(set([k.species for k in val.getListOfReactants()]+[k.species for k in val.getListOfProducts()])) \
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
    self.match_score = None
    self.sum_match_score = None
    self.query_df = None
    # self.one_candidates = None

  # def getMatchScore(self, score_dict):
  #   """
  #   Calculate match scores using 
  #   the given dictionary
  
  #   Parameters
  #   ----------
  #   inp_match_score: dict {reaction_id: {candidate_id: float}}
  
  #   Returns
  #   -------
  #   match_score: float
  #   """
  #   return np.sum([np.max([score_dict[val][k] \
  #                          for k in score_dict[val].keys()]) \
  #                  for val in score_dict.keys()])

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
    r_components: None/str-list (list of species IDs)
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

  def predictAnnotation(self,
                        inp_spec_dict,
                        inp_reac_list=None,
                        inp_ref_mat=ref_mat,
                        update=False):
    """
    Predict 1) reaction annotation candidates 
    and 2) match score of them
    using species dict (argument) etc.
    Also updates several class attributes as a result. 
  
    Parameters
    ----------
    dict: inp_spec_dict
        Dictionoary, {species id: formula(str-list)}
    inp_reac_list: str-list
        IDs of reactions to predict. If default, will do all reactions
    inp_ref_mat: pd.DataFrame
        Reference matrix
    update: bool
        Whether to save results as class attrributes or just return them.
        If True, only match_score is returned (and values are updated).
        If False, all relevant information is returned

    Returns
    -------
    (If inplace=True)
    pred_match_score: dict
        Match score of each prediction
        {reaction ID: (Rhea ID, match score: float between 0.0-1.0)}
    (If inplace=False)
    res_dict: dict
    """
    # get libsbml.reaction and their IDs
    if inp_reac_list is not None:
      reaction_ids = inp_reac_list
    else:
      reaction_ids = list(self.reaction_components.keys())
    # get dictionary of reaction ID: species component
    r2pred_spec_formulas = dict()
    for one_rid in reaction_ids:
      r2pred_spec_formulas[one_rid] = {one_spec:inp_spec_dict[one_spec] \
                                       for one_spec in self.reaction_components[one_rid]}
    # prepare query df for prediction
    query_df = pd.DataFrame(0, 
                            index=inp_ref_mat.columns,
                            columns=reaction_ids)
    for one_rid in reaction_ids:
      one_set_species = r2pred_spec_formulas[one_rid]
      # for each species element of the select reaction
      for one_spec_key in one_set_species.keys():
        one_spec = one_set_species[one_spec_key]
        # For each one_rid, set the values 1.0
        query_df.loc[[val for val in one_spec if val in query_df.index], one_rid] = 1
    multi_mat = inp_ref_mat.dot(query_df)
    maxes = multi_mat.max()
    #
    # Collect candidates and calculate confidence score
    pred_cands = dict()
    pred_match_score = dict()
    for one_rid in maxes.index:
      one_multi = multi_mat.loc[:,one_rid]
      candidates = one_multi[one_multi==maxes[one_rid]].index
      # cand_data; (number of element matches, candidates)
      pred_cands[one_rid] = candidates
      # Now, match_scpre (calculated per each candidate) => replaced as a tuple
      match_score_per_cand = []
      for one_cand in candidates:
        if one_cand in ref_rhea2bi.keys():
          num_matches = maxes[one_rid]
          num_maxpos_matches = len(inp_ref_mat.loc[one_cand, :].to_numpy().nonzero()[0])
          match_score_per_cand.append((one_cand, num_matches/num_maxpos_matches))
      match_score_per_cand.sort(key=operator.itemgetter(1), reverse=True)
      pred_match_score[one_rid] = match_score_per_cand
    if update:
      self.candidates = pred_cands
      self.match_score = pred_match_score
      self.query_df = query_df
      # self.one_candidates = self.getBestOneCandidates(self.match_score)
      return pred_match_score
    else:
      return {'candidates': pred_cands,
              'match_score': pred_match_score,
              'query_df': query_df}
              # 'one_candidates': self.getBestOneCandidates(self.match_score)}

  # # Changed as match_score was changed to a tuple (rhea_term, match_score)
  # def getBestOneCandidates(self, inp_match_score=None):
  #   """
  #   Get a dictinoary of {reaction_id: [single candidate]}.
  #   If self.predictAnnotation should have been already run. 

  #   Parameters
  #   ----------
  #   inp_match_score: dict
  #       {reaction_id: {Rhea_id: match_score(i.e., float)}}

  #   Returns
  #   -------
  #   ranked_one_cands: dict
  #       {reaction_id: [one Rhea_id]}
  #   """
  #   if inp_match_score is None:
  #     match_score = self.match_score
  #   else:
  #     match_score = inp_match_score
  #   ranked_one_cands = dict()
  #   for one_k in match_score.keys():
  #     # already sorted; just need to get the first element. 
  #     ranked_one_cands[one_k] = [match_score[one_k][0][0]]
  #     # one_itm = pd.DataFrame.from_dict(match_score[one_k], orient='index', columns=['match_score'])
  #     # one_itm.sort_values(ascending=False, by='match_score', inplace=True)
  #     # ranked_one_cands[one_k] = [one_itm.index[0]]
  #   return ranked_one_cands

  # TODO: will be updated when iteration gets updated. 
  # def updateSpeciesByAReaction(self, 
  #                              inp_rid, inp_spec_dict,
  #                              inp_rhea, inp_ref_mat=ref_mat):
  #   """
  #   Update predicted species annotation
  #   using predicted rhea annotation candidate
  #   for a single reaction.
  #   Current version works for when there is just
  #   one candidate (i.e., one RHEA term).  
  
  #   Parameters
  #   ----------
  #   inp_rid: str
  #       Reactino ID to match the candidate
  #   inp_spec_dict: dict
  #       {species ID: chemical formula}
  #   inp_rhea: str
  #       A RHEA term predicted 
  #   inp_ref_mat: pandas.DataFrame
  
  #   Returns
  #   -------
  #   dict: {species ID: [list of CHEBI terms]}
  #         Suggested mapping from species ID to a chebi term
  #   """
  #   # Chebi terms (reference) associated with the given Rhea term
  #   # If there is no such Rhea term in the reference, return None
  #   if inp_rhea in ref_rhea2bi.keys():
  #     if ref_rhea2bi[inp_rhea] in ref_rhea_to_chebi.keys():
  #       rhea_term_to_chebi_elements = [val for val in ref_rhea_to_chebi[ref_rhea2bi[inp_rhea]] \
  #                                      if val in cn.ref_chebi2formula.keys()]
  #     else:
  #       return None
  #   else:
  #     return None
  #   # Similarly, species formula (reference) associated with the Rhea term. 
  #   one_r_elements_row = inp_ref_mat.loc[inp_rhea, :]
  #   one_r_elements = one_r_elements_row[one_r_elements_row!=0].index
  #   # Next, species formulas (query) associated with the predicted Rhea term.
  #   one_r_query_elements_row = self.query_df.loc[:, inp_rid]
  #   one_r_query_elements = one_r_query_elements_row[one_r_query_elements_row!=0].index
  #   #
  #   # Identify species not included vs. species included (in query), compared to ref_mat
  #   # species_not_included = [val for val in one_r_elements if val not in one_r_query_elements]
  #   formula_included = [val for val in one_r_elements if val in one_r_query_elements]
  #   # {predicted species ID: formula} for all elements in the reaction 
  #   all_species_in_a_reaction = self.getReactionComponents(inp_rid)
  #   spec2predicted_formula = {one_spec:inp_spec_dict[one_spec] \
  #                                      for one_spec in all_species_in_a_reaction}  
  #   #
  #   possibly_correct_species_ids = []
  #   chebi_term_already_used_in_prediction = []
  #   for one_incl_formula in formula_included:
  #     for one_specid in spec2predicted_formula.keys():
  #       # If one included formula can be found in one of the predicted spec_dict,
  #       # it will be possibly correct (we are looking for incorrect ones)
  #       if one_incl_formula in spec2predicted_formula[one_specid]:
  #         possibly_correct_species_ids.append(one_specid)
  #         for one_chebi_term in rhea_term_to_chebi_elements:
  #           if one_incl_formula == cn.ref_chebi2formula[one_chebi_term]:
  #             chebi_term_already_used_in_prediction.append(one_chebi_term)
  #   remaining_chebi = [val for val in rhea_term_to_chebi_elements \
  #                      if val not in chebi_term_already_used_in_prediction]
  #   remaining_specid = [val for val in all_species_in_a_reaction \
  #                       if val not in possibly_correct_species_ids]
  #   if len(remaining_specid) == 1 and len(remaining_chebi) == 1:    
  #     return {remaining_specid[0]: [remaining_chebi[0]]}
  #   else:
  #     None

  def getAccuracy(self,
                  ref_annotation=None,
                  pred_annotation=None):
    """
    Compute accuracy of reaction annotation.
    A list of annotations of 
    a single reaction (identified by each ID) 
    is considered accurate if it includes
    the corresponding value of ref_annotation.
    (More precisely, if there is at least one
    intersection).
    Once prediction is run, self.candidates
    can be used for pred_annotation. 
  
    Parameters
    ----------
    ref_annotation: dict
        {reaction_id: [str-annotation]}
        if None, get self.exist_annotation
    pred_annotation: dict
        {reaction_id: [str-annotation]}
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
      pred = self.candidates
    else:
      pred = pred_annotation
    ref_keys = set(ref.keys())
    pred_keys = set(pred.keys())
    reactions_to_test = ref_keys.intersection(pred_keys)
    for one_k in reactions_to_test:
      if set(ref[one_k]).intersection(pred[one_k]):
        accuracy.append(True)
      else:
        accuracy.append(False)
    return np.mean(accuracy)

  # Develop a method to evaluate results using fitted model
  def evaluatePredictedReactionAnnotation(self, inp_dict,
                                          fitted_model=reaction_rf):
    """
    Evaluate the quality of annotation;
    for each individual species.
  
    Parameters
    ---------
    inp_dict: dict
        {Reaction ID: number of candidatees}

    Returns
    -------  
    res: dict {reaction_id: probability-of-species-prediction-being-correct}
        Information of whether confident or not
    """
    # candidates_info = self.candidates
    cands_num_dict = {one_k: len(inp_dict[one_k]) for one_k in inp_dict.keys()}
    inp_list = list(inp_dict.keys())
    num_candidates = [cands_num_dict[val] for val in inp_list]
    # num_candidates = [len(candidates_info[val]) for val in inp_list]
    multi_mat = ref_mat.dot(self.query_df)
    maxes = multi_mat.max()
    max_match = [maxes[val] for val in inp_list]
    match_scores = self.match_score
    mean_match_score = [np.mean([val[1] for val in match_scores[k]]) for k in inp_list]
    med_match_score = [np.median([val[1] for val in match_scores[k]]) for k in inp_list]
    min_match_score = [np.min([val[1] for val in match_scores[k]]) for k in inp_list]
    max_match_score = [np.max([val[1] for val in match_scores[k]]) for k in inp_list]
    var_match_score = [np.var([val[1] for val in match_scores[k]]) for k in inp_list]
    data2prediction = list(zip(num_candidates,
                               max_match,
                               mean_match_score,
                               med_match_score,
                               min_match_score,
                               max_match_score,
                               var_match_score))
    pred_probs = fitted_model.predict(data2prediction)
    # Collect probability to be correct
    res = {val[0]:val[1] for val in list(zip(inp_list, pred_probs))}
    return res

