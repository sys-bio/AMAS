#!/usr/bin/env python

# recommend_reactions.py
"""
Predicts annotations of reaction(s) using a local XML file
and the reaction ID. 
Usage: python recommend_reaction.py files/BIOMD0000000190.xml --min_score 0.6 --outfile res.csv
"""

import argparse
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='Recommend reaction annotations of an SBML model and save results') 
  parser.add_argument('model', type=str, help='SBML model file (.xml)')
  # One or more reaction IDs can be given
  parser.add_argument('--reactions', type=str, help='ID(s) of reaction(s) to be recommended. ' +\
                                                   'If not provided, all reactions will be used', nargs='*')
  parser.add_argument('--min_len', type=int, help='Minimum number of reaction components (reactants and products) ' +\
                                                  'to be used for prediction. ' +\
                                                  'Reactions with fewer components than this value ' +\
                                                  'will be ignored. Default is zero.', nargs='?', default=0)
  parser.add_argument('--cutoff', type=float, help='Match score cutoff.', nargs='?', default=0.0)
  parser.add_argument('--method', type=str, 
                                  help='Choose either "top" or "above". "top" recommends ' +\
                                       'the best candidates that are above the cutoff, ' +\
                                       'and "above" recommends all candidates that are above ' +\
                                       'the cutoff. Default is "top".',
                                  nargs='?',
                                  default='top')
  parser.add_argument('--outfile', type=str, help='File path to save recommendation.', nargs='?',
                      default=os.path.join(os.getcwd(), 'reaction_rec.csv'))
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  one_fpath = args.model
  reacts = args.reactions
  min_len = args.min_len
  cutoff = args.cutoff
  method = args.method
  outfile = args.outfile

  recom = recommender.Recommender(libsbml_fpath=one_fpath)
  recom.current_type = 'reaction'
  # if nothing is given, predict all IDs
  if reacts is None:
    reacts = recom.getReactionIDs()
  print("...\nAnalyzing %d reaction(s)...\n" % len(reacts))
  # choosing reactions with at least 'min_len' components
  filt_reacts = [val for val in reacts \
                 if len(recom.reactions.reaction_components[val]) >= min_len]
  # stops if all elements were removed by filtering...
  if len(filt_reacts) == 0:
    print("No element found after the element filter.")
    return None
  res = recom.getReactionListRecommendation(pred_ids=filt_reacts, get_df=True)
  for idx, one_df in enumerate(res):
    filt_df = recom.autoSelectAnnotation(df=one_df,
                                         min_score=cutoff,
                                         method=method)
    recom.updateSelection(filt_reacts[idx], filt_df)
  # save file to csv
  recom.saveToCSV(outfile)
  print("Recommendations saved as:\n%s\n" % os.path.abspath(outfile))

if __name__ == '__main__':
  main()