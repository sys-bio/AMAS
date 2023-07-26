#!/usr/bin/env python

# recommend_annotation.py
"""
Predicts annotations of species and reactions using a local XML file
and the reaction ID. 
This is a combined version of recommend_species and recommend_reaction,
but is more convenient because user will just get the updated XML file or whole recommendations. 
Usage: python recommend_reaction.py files/BIOMD0000000190.xml --cutoff 0.6 --save csv --outfile res.csv 
"""

import argparse
import itertools
import libsbml
import numpy as np
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import iterator as it
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='Recommend annotations of an SBML model ' +\
                                               '(for both species and reactions) and save results.') 
  parser.add_argument('model', type=str, help='SBML model file (.xml).')
  # One or more reaction IDs can be given
  parser.add_argument('--cutoff', type=float, help='Match score cutoff.', nargs='?', default=0.0)
  parser.add_argument('--optimize', type=str, help='Whether to optimize or not. ' +\
                                                   'If y or yes is given, predictions will be ' +\
                                                   'optimized. N or no will not optimize predictions.',
                                              nargs='?',
                                              default='no')
  parser.add_argument('--mssc', type=str,
                                help='Match score selection criteria (MSSC). ' +\
                                     'Choose either "top" or "above". "top" recommends ' +\
                                     'the best candidates that are above the cutoff, ' +\
                                     'and "above" recommends all candidates that are above ' +\
                                     'the cutoff. Default is "top"',
                                nargs='?',
                                default='top')
  parser.add_argument('--save', type=str, 
                                help='Either "sbml" or "csv". ' +\
                                     'If "sbml" is chosen, model will be automatically ' +\
                                     'annotated with recommended candidates and saved. ' +\
                                     'If "csv" is chosen, recommendations will be saved ' +\
                                     'as a csv file. Default is "csv".',
                                nargs='?',
                                default='csv')
  parser.add_argument('--outfile', type=str, help='Path to save an output file.', nargs='?')
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  one_fpath = args.model
  cutoff = args.cutoff
  optim_raw = args.optimize
  if optim_raw.lower() in ['y', 'yes']:
    optim = True
  else:
    optim = False
  mssc = args.mssc.lower()
  save = args.save
  outfile = args.outfile
  #
  recom = recommender.Recommender(libsbml_fpath=one_fpath)
  specs = recom.getSpeciesIDs()
  print("...\nAnalyzing %d species...\n" % len(specs))
  reacts = recom.getReactionIDs()
  print("...\nAnalyzing %d reaction(s)...\n" % len(reacts))

  res_tab = recom.recommendAnnotation(mssc=mssc,
                                      cutoff=cutoff,
                                      optimize=optim,
                                      outtype='table')
  if save == 'csv':
    if outfile is None:
      outfile = os.path.join(os.getcwd(), 'recommendations.csv')
    recom.saveToCSV(res_tab, outfile)
  else:
    if outfile is None:
      outfile = os.path.join(os.getcwd(), 'updated_model.xml')
    res_sbml = recom.getSBMLDocument(sbml_document=recom.sbml_document,
                                     chosen=res_tab,
                                     auto_feedback=True)
    libsbml.writeSBMLToFile(res_sbml, outfile)   
  
  print("Recommendations saved as:\n%s\n" % os.path.abspath(outfile))



  # if save.lower() == 'sbml':
  #   if outfile is None:
  #     fin_path = os.path.join(os.getcwd(), 'updated_model.xml')
  #   else:
  #     fin_path = outfile
  #   recom.saveToSBML(fin_path)
  #   print("Annotations updated and saved as:\n%s\n" % os.path.abspath(fin_path))
  # elif save.lower() == 'csv':
  #   if outfile is None:
  #     fin_path = os.path.join(os.getcwd(), 'recommendations.csv')
  #   else:
  #     fin_path = outfile
  #   recom.saveToCSV(fin_path)



  #   res_spec = recom.getSpeciesListRecommendation(pred_ids=specs)

  # res_reac = recom.getReactionListRecommendation(pred_ids=reacts,
  #                                                spec_res=res_spec)
  # if optim:
  #   print("Optimizing predictions...\n")
  #   anot_iter = it.Iterator(cur_spec_formula=recom.species.formula,
  #                           reaction_cl=recom.reactions,
  #                           reactions_to_update=reacts)
  #   res_iter = anot_iter.match()
  #   recoms_tobe_added = []
  #   for one_spec in res_iter.keys():
  #     pred_reacs = [val.id for val in res_reac]
  #     reacs_using_one_spec = [val for val in pred_reacs \
  #                             if one_spec in recom.reactions.reaction_components[val]]
  #     filt_res_reac = [val for val in res_reac if val.id in reacs_using_one_spec]
  #     # match score of reactions using that species
  #     adj_match_score = np.mean(list(itertools.chain(*[[cand[1] for cand in val.candidates] \
  #                       for val in filt_res_reac])))
  #     cands = res_iter[one_spec]
  #     adj_formulas = list(set([cn.REF_CHEBI2FORMULA[k] \
  #                              for k in cands if k in cn.REF_CHEBI2FORMULA.keys()]))
  #     urls = [cn.CHEBI_DEFAULT_URL + val[6:] for val in cands]
  #     labels = [cn.REF_CHEBI2LABEL[val] for val in cands]
  #     adj_recom = cn.Recommendation(one_spec,
  #                                   [(val, adj_match_score) for val in cands],
  #                                   urls,
  #                                   labels)
  #     recoms_tobe_added.append(adj_recom)
  #   fin_spec_recom = recoms_tobe_added + \
  #                     [val for val in res_spec if val.id not in res_iter.keys()]
  #   fin_reac_recom = recom.getReactionListRecommendation(pred_ids=reacts,
  #                                                        spec_res=fin_spec_recom)
  # else:
  #   fin_spec_recom = res_spec
  #   fin_reac_recom = res_reac

  # RECOM_DICT = {'species': fin_spec_recom,
  #               'reaction': fin_reac_recom}
  # for one_k in RECOM_DICT:
  #   recom.current_type = one_k
  #   for one_recom in RECOM_DICT[one_k]:
  #     filt_df = recom.autoSelectAnnotation(df=recom.getDataFrameFromRecommendation(one_recom),
  #                                          cutoff=cutoff,
  #                                          mssc=mssc)
  #     recom.updateSelection(one_recom.id, filt_df)  
  # # save file
  # if save.lower() == 'sbml':
  #   if outfile is None:
  #     fin_path = os.path.join(os.getcwd(), 'updated_model.xml')
  #   else:
  #     fin_path = outfile
  #   recom.saveToSBML(fin_path)
  #   print("Annotations updated and saved as:\n%s\n" % os.path.abspath(fin_path))
  # elif save.lower() == 'csv':
  #   if outfile is None:
  #     fin_path = os.path.join(os.getcwd(), 'recommendations.csv')
  #   else:
  #     fin_path = outfile
  #   recom.saveToCSV(fin_path)
  #   print("Recommendations saved as:\n%s\n" % os.path.abspath(fin_path))


if __name__ == '__main__':
  main()