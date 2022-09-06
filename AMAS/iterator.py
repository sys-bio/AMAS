# iterator.py
"""
Script for iteration;
can be refactored into class in future..
"""

import os
import pickle

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra

with open(os.path.join(cn.CHEBI_DIR, 'chebi_shortened_formula_30apr2022.pickle'), 'rb') as f:
  ref_shortened_chebi_to_formula = pickle.load(f)


def iterateAndGetUpdatedResults(spec_cl,
                                reac_cl,
                                num_iter=10,
                                show_message=False):
  """
  Using current species and reaction annotations
  (for species, get formulafied one),
  update both species and reaction annotations
  
  Parameters
  ----------
  spec_cl: species_annotation.SpeciesAnnotation class
  reac_cl: reaction_annotation.ReactionAnnotation class
  num_iter: int
  show_message: bool
      Prints out iteration number and match score
      
  Returns
  -------
  result: dict
      e.g., rep (int) 
          Last iteration position (actual run - 1) when algorithm quit
  """
  cur_candidates_dict = reac_cl.candidates
  cur_reac_match_score = reac_cl.sum_match_score
  cur_one_cands = reac_cl.one_candidates
  cur_spec_formula_dict = spec_cl.formula

  flag = False
  if show_message:
    print("Initial match score: %.02f" % cur_reac_match_score)
    print("*************************")
  for rep in range(0, num_iter):
    if flag:
      break
    if show_message:
      print("Iteration %d" % (rep+1))
    # updated chebi values
    all_upd_spec = dict()
    for one_k in cur_one_cands.keys():
      one_upd_spec = reac_cl.updateSpeciesByAReaction(inp_rid=one_k,
                                                      inp_spec_dict=cur_spec_formula_dict,
                                                      inp_rhea=cur_one_cands[one_k][0],
                                                      inp_ref_mat=ra.ref_mat)
      all_upd_spec = tools.updateDictKeyToList(all_upd_spec, one_upd_spec)

    # update species dictionary to use
    upd_spec_formula_dict = dict()
    for one_k in cur_spec_formula_dict.keys():
      if one_k in all_upd_spec.keys():
        upd_spec_formula_dict[one_k] = list(set([ref_shortened_chebi_to_formula[val] for val in all_upd_spec[one_k]]))
      else:
        upd_spec_formula_dict[one_k] = cur_spec_formula_dict[one_k]

    # Using upd_spec_formula_dict, predict reaction again 
    upd_reac_annotation = reac_cl.predictAnnotation(inp_spec_dict=upd_spec_formula_dict,
                                           inp_reac_list=None,
                                           inp_ref_mat=ra.ref_mat,
                                           update=False)
    upd_reac_match_score = upd_reac_annotation['sum_match_score']
    # Check wheter to continue;
    if upd_reac_match_score > cur_reac_match_score:
      cur_candidates_dict = upd_reac_annotation['candidates']
      cur_one_cands = upd_reac_annotation['one_candidates']
      cur_spec_formula_dict = upd_spec_formula_dict
      cur_reac_match_score = upd_reac_match_score
      if show_message:
        print("Updated match score: %.02f" % cur_reac_match_score)
        print("*************************")
    else:
      flag = True
      if show_message:
        print("Updated match score: %.02f" % cur_reac_match_score)
        print("Score not increasing. Quitting iteration...")
  if show_message:
    print("\nCalculation finished.")
  result = {'candidates': cur_candidates_dict,
            'spec_formula': cur_spec_formula_dict,
            'spec_chebi2update': all_upd_spec,
            'sum_match_score': cur_reac_match_score,
            'rep': rep}
  return result
















