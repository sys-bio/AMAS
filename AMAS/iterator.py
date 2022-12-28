# iterator.py
"""
Iterator class will be used by 
the recommender; 
Recommener first predicts annotations
and (if it wants to use iteration)
sends it to the iterator,
whhich uses two dictionaries as arguments. 
Iterator tries and tests whether
iteration really improves the (match) scores,
and if so, informs recommender 
that it should (or can) update the
'current' annotations.

Currently, most up-to-date annotations
are stored in 
species.candidates & species.formula, (for speices)
and reactions.candidates (for reactions). 
"""

import copy
import itertools
import numpy as np
import os

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra


# Keys when evaluating match results. 
NEW_SCORE = 'new_score'
OLD_SCORE = 'old_score'
INCREASED = 'is_increased'
# Max limit for iteration
MAX_ITER = 3


class Iterator(object):

  def __init__(self,
               cur_spec_formula,
               reaction_cl,
               reactions_to_update=None):
    """
    Ideally, arguments should be directly from
    the relevant species.formula and reactions.candidates.
    Returns may be dictionaries 
    for species.candidates and reactions.candidates;
    (not 100% confirmed yet)

    Parameters
    ----------
    cur_spec_formula: dict
        {species_id: [predicted-formulas]}
        Current (most recent) annotation of species
    reaction_cl: AMAS.reaction_annotation.ReactionAnnotation
        reaction_annotation class instance with
        loaded SBML model information (reaction_components, etc.)
    reaction_to_update: list-str
        List of reactions to update; if None, use all reactions
        from self.reactions.candidates
    """
    self.orig_spec_formula = cur_spec_formula
    # Storing reaction candidates separately, 
    # as it may be different than self.reactions.candidates
    self.reactions = reaction_cl
    if reactions_to_update:
      self.r2upd = reactions_to_update
    else:
      self.r2upd = list(reaction_cl.candidates.keys())

  def getDictOfRheaComponentFormula(self, inp_rhea):
    """
    Get a dictionary {chebi_id: formula}
    from a given rhea term.
    Rhea term -> CheBI IDs -> Formulas
  
    Parameters
    ----------
    str: inp_rhea
        A Rhea identifier
  
    Returns
    -------
    : dict
        {chebi_id: formula-str}
    """
    chebis = cn.REF_RHEA2CHEBI[inp_rhea]
    return {val:cn.REF_CHEBI2FORMULA[val] for val in chebis \
            if val in cn.REF_CHEBI2FORMULA.keys()}


  def getDictMatchByItem(self,
                         chebi2ref_formula,
                         spec2pred_formula):
    """
    Get match between two keys,
    where there are exactly 
    one matching items.
    If all items are matched by 1-1
    (i.e., one species - one chebi),
    return the fully matched dictionary.
    (i.e., improve precision)
    If neither, return None.
    (i.e., nothing to update)
  
    Parameters
    ----------
    chebi2ref_formula: dict
        {chebi_term: a_species_formula(string)}
    spec2pred_formula: dict
        {species_id: [predicted_formulas]}
  
    Returns
    -------
    dict/None
        {species_id: [chebi_term]}
    """
    match_dict = {one_k:[spec_id for spec_id in spec2pred_formula.keys() \
                         if chebi2ref_formula[one_k] in spec2pred_formula[spec_id]
                        ] \
                  for one_k in chebi2ref_formula.keys()}
    unmatched_species = [val for val in spec2pred_formula.keys() \
                        if val not in list(itertools.chain(*match_dict.values()))]
    unmatched_chebi = [val for val in match_dict.keys() if not match_dict[val]]
    if len(unmatched_species) == 1 and len(unmatched_chebi) == 1:
      return {unmatched_species[0]: unmatched_chebi} 
    # reverse match_dict into the proper return format. 
    elif all([len(val[1])==1 for val in list(match_dict.items())]):
      return {match_dict[k][0]: [k] for k in match_dict.keys()}
    else:
      return None


  def getDictsToUpdate(self, reaction_id):
    """
    Using self.getDictMatchByItem(),
    get dictionaries to update

    Parameters
    ----------
    str: reaction_id

    Returns
    -------
    match_res: dict
        {species_id: [ChEBI terms]}

    match_res_formula: dict
        {species_id: [formula-str]}
    """
    one_rhea = self.reactions.candidates[reaction_id][0][0]
    # match_res will look like {species_id: [CHEBI term]}
    # filter to have only keys and items of one reaction
    filt_spec_formula = {k:self.orig_spec_formula[k] \
                           for k in self.reactions.reaction_components[reaction_id]}
    upd_spec_chebi = self.getDictMatchByItem(chebi2ref_formula=self.getDictOfRheaComponentFormula(one_rhea),
                                             spec2pred_formula=filt_spec_formula)
    if upd_spec_chebi:
      upd_spec_formula = {k:[cn.REF_CHEBI2FORMULA[chebi] \
                          for chebi in upd_spec_chebi[k]] for k in upd_spec_chebi.keys()}
    else:
      upd_spec_formula = None
    return upd_spec_chebi, upd_spec_formula


  def getUpdatedMatchScore(self, cur_spec_formulas, inp_spec2formula_dict):
    """
    Check whether it improves reaction measures; 
    if new value (sum of maximum match score per reaction)
    increased, return True; otherwise return False.
  
    Parameters
    ----------
    cur_spec_formulas: dict
        {'species_id': [formula-str]}
        Dictionary to be updated
      
    inp_spec_2formula_dict: dict
        {'species_id': [formula-str]}
        Dictionary to update
      
    Returns
    -------
    : dict
    """
    cur_spec_formulas.update(inp_spec2formula_dict)
    new_pred_res = self.reactions.predictAnnotation(inp_spec_dict = cur_spec_formulas,
                                                    inp_reac_list = list(self.r2upd))
    old_pred_res = self.reactions.predictAnnotation(inp_spec_dict = self.orig_spec_formula,
                                                    inp_reac_list = list(self.r2upd))
    # since candidates are already sorted, 
    # just check the match score (index '1') of the very first candidate tuple (index '0')
    new_pred_val = np.mean([new_pred_res[cn.MATCH_SCORE][k][0][1] \
                           for k in new_pred_res[cn.MATCH_SCORE].keys()])
    old_pred_val = np.mean([old_pred_res[cn.MATCH_SCORE][k][0][1] \
                           for k in old_pred_res[cn.MATCH_SCORE].keys()])
    return {NEW_SCORE: new_pred_val,
            OLD_SCORE: old_pred_val,
            INCREASED: new_pred_val>old_pred_val}

  def match(self):
    """
    Use self.runOneMatchCycle()
    and determine the final products to return.
    Will be used by the recommender or the user. 
    """
    all_upd_spec_chebi = dict()
    for _ in range(MAX_ITER):
      upd_spec_chebi = self.runOneMatchCycle()
      if upd_spec_chebi:
        all_upd_spec_chebi.update(upd_spec_chebi)
        # Update the formula attribute for the next iteration
        for one_k in upd_spec_chebi.keys():
          self.orig_spec_formula[one_k] = [cn.REF_CHEBI2FORMULA[val] \
                                           for val in upd_spec_chebi[one_k] \
                                           if val in cn.REF_CHEBI2FORMULA.keys()]
      else:
        break
    # Maybe run reaction once, and return final results :) 
    return all_upd_spec_chebi

  def runOneMatchCycle(self):
    """
    Using the methohds & information,
    determine species to update. 
    (Reaction will be updated in the following steps). 
    This method will directly used by
    the Recommender, or even the user. 

    Returns
    -------
    combine_upd_spec2chebi: dict
        {species_id: [ChEBI terms]}
    """
    combine_upd_spec2chebi = dict()
    # Use reactions existing in self.r2upd
    for one_reaction in self.r2upd:
      one_rhea_tup = self.reactions.candidates[one_reaction]
      one_rhea = one_rhea_tup[0][0]
      pred_spec_formulas = self.orig_spec_formula
      one_rhea2formula = self.getDictOfRheaComponentFormula(inp_rhea=one_rhea)
      upd_spec2chebi, upd_spec2formula = self.getDictsToUpdate(reaction_id=one_reaction)
      # Meaning, when examining match scores we only consider 
      # individual updates; not cumulated updtaes (so we don't use combine_spec2chhebi below)
      if upd_spec2formula: 
        upd_val = self.getUpdatedMatchScore(cur_spec_formulas = copy.deepcopy(self.orig_spec_formula),
                                            inp_spec2formula_dict = upd_spec2formula)

        if upd_val[INCREASED]:
          # update combine_upd_spec2chebi by combining the elements.
          for k in upd_spec2chebi.keys():
            if k in combine_upd_spec2chebi.keys():
              combine_upd_spec2chebi[k] = list(set(combine_upd_spec2chebi[k] + upd_spec2chebi[k]))
            else:
              combine_upd_spec2chebi[k] = upd_spec2chebi[k] 
    return combine_upd_spec2chebi









