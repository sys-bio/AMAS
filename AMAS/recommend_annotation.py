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
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender

def main():
  parser = argparse.ArgumentParser(description='Recommend annotations of an SBML model (for both species and reactions) and save results') 
  parser.add_argument('model', type=str, help='SBML model file (.xml)')
  # One or more reaction IDs can be given
  parser.add_argument('--cutoff', type=float, help='minimum match score cutoff', nargs='?', default=0.0)
  parser.add_argument('--method', type=str,
                                  help='Choose either "top" or "above". "top" recommends ' +\
                                       'the best candidates that are above the cutoff, ' +\
                                       'and "above" recommends all candidates that are above ' +\
                                       'the cutoff. Default is "top"',
                                  nargs='?',
                                  default='top')
  parser.add_argument('--save', type=str, 
                                help='Choose either "sbml" or "csv". ' +\
                                     'If "sbml" is chosen, model will be automatically ' +\
                                     'annotated with recommended candidates and saved. ' +\
                                     'If "csv" is chosen, recommendations will be saved ' +\
                                     'as a csv file. Default is "csv"',
                                nargs='?',
                                default='csv')
  parser.add_argument('--outfile', type=str, help='path to save file.', nargs='?')
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  one_fpath = args.model
  cutoff = args.cutoff
  method = args.method
  save = args.save
  outfile = args.outfile
  try:
    recom = recommender.Recommender(libsbml_fpath=one_fpath)
    recom.current_type = 'species'
    specs = recom.getSpeciesIDs()
    print("...\nAnalyzing %d species...\n" % len(specs))
    res_spec = recom.getSpeciesListRecommendation(pred_ids=specs, get_df=True)
    for idx, one_df in enumerate(res_spec):
      filt_df = recom.autoSelectAnnotation(df=one_df,
                                           min_score=cutoff,
                                           method=method)
      recom.updateSelection(specs[idx], filt_df)

    recom.current_type = 'reaction'
    reacts = recom.getReactionIDs()
    print("...\nAnalyzing %d reaction(s)...\n" % len(reacts))
    res_reac = recom.getReactionListRecommendation(pred_ids=reacts, get_df=True)
    for idx, one_df in enumerate(res_reac):
      filt_df = recom.autoSelectAnnotation(df=one_df,
                                           min_score=cutoff,
                                           method=method)
      recom.updateSelection(reacts[idx], filt_df)
    # save file
    if save.lower() == 'sbml':
      if outfile is None:
        fin_path = os.path.join(os.getcwd(), 'updated_model.xml')
      else:
        fin_path = outfile
      recom.saveToSBML(fin_path)
      print("Annotations updated and saved as:\n%s\n" % os.path.abspath(fin_path))
    elif save.lower() == 'csv':
      if outfile is None:
        fin_path = os.path.join(os.getcwd(), 'recommendations.csv')
      else:
        fin_path = outfile
      recom.saveToCSV(fin_path)
      print("Recommendations saved as:\n%s\n" % os.path.abspath(fin_path))
  except:
    raise ValueError("Please check arguments.")


if __name__ == '__main__':
  main()