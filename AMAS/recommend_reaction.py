# recommend_reaction.py
"""
Predicts annotations of reaction(s) using a local XML file
and the reaction ID. 
Usage: python recommend_reaction.py --model files/BIOMD0000000190.xml --min_score 0.6 --out_dir res
"""

import argparse
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='SBML file (.XML) and one or more reaction IDs in the model.') 
  parser.add_argument('model', type=str, help='SBML model file in the XML format')
  # One or more reaction IDs can be given
  parser.add_argument('--reaction', type=str, help='ID of reaction(s) in the model', nargs='*')
  parser.add_argument('--min_score', type=float, help='Minimum threshold', nargs='?', default=0.0)
  parser.add_argument('--method', type=str, help='Choose either "best" or "all". Default is "best".', nargs='?', default='best')
  parser.add_argument('--out_dir', type=str, help='Path of directory to save files', nargs='?', default=os.getcwd())
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  one_fpath = args.model
  reacts = args.reaction
  min_score = args.min_score
  method = args.method
  out_dir = args.out_dir
  try:
    recom = recommender.Recommender(libsbml_fpath=one_fpath)
    recom.current_type = 'reaction'
    # if nothing is given, predict all IDs
    if reacts is None:
      reacts = recom.getReactionIDs()
    print("\nAnalyzing %d reaction(s)...\n" % len(reacts))
    res = recom.getReactionListRecommendation(pred_ids=reacts, get_df=True)
    for idx, one_df in enumerate(res):
      filt_df = recom.autoSelectAnnotation(df=one_df,
                                           min_score=min_score,
                                           method=method)
      recom.updateSelection(reacts[idx], filt_df)
    # Create a new directory if it doesn't exist already
    if not os.path.exists(out_dir):
      os.mkdir(out_dir)
    recom.saveToCSV(os.path.join(out_dir, 'reaction_recommendation.csv'))
    recom.saveToSBML(os.path.join(out_dir, 'model_amas_reactions.xml'))
  except:
    raise ValueError("Please check arguments.")


if __name__ == '__main__':
  main()