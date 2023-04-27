#!/usr/bin/env python

# recommend_species.py
"""
Predicts annotations of species using a local XML file
and the species ID. 
Usage: python recommend_species.py files/BIOMD0000000190.xml --min_score 0.6 --outfile res.csv
"""

import argparse
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='Recommend species annotations of an SBML model and save results') 
  parser.add_argument('model', type=str, help='SBML model file (.xml)')
  # One or more species IDs can be given
  parser.add_argument('--species', type=str, help='ID(s) of species to be recommended. If not provided, all species will be used', nargs='*')
  parser.add_argument('--cutoff', type=float, help='minimum match score cutoff', nargs='?', default=0.0)
  parser.add_argument('--method', type=str,
                                  help='Choose either "top" or "above". "top" recommends ' +\
                                       'the best candidates that are above the cutoff, ' +\
                                       'and "above" recommends all candidates that are above ' +\
                                       'the cutoff. Default is "top"',
                                  nargs='?',
                                  default='top')
  parser.add_argument('--outfile', type=str, help='file path to save recommendation', nargs='?',
                      default=os.path.join(os.getcwd(), 'species_rec.csv'))
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  one_fpath = args.model
  specs = args.species
  cutoff = args.cutoff
  method = args.method
  outfile = args.outfile

  try:
    recom = recommender.Recommender(libsbml_fpath=one_fpath)
    recom.current_type = 'species'
    # if nothing is given, predict all IDs
    if specs is None:
      specs = recom.getSpeciesIDs()
    print("...\nAnalyzing %d species...\n" % len(specs))
    res = recom.getSpeciesListRecommendation(pred_ids=specs, get_df=True)
    for idx, one_df in enumerate(res):
      filt_df = recom.autoSelectAnnotation(df=one_df,
                                           min_score=cutoff,
                                           method=method)
      recom.updateSelection(specs[idx], filt_df)
    # save file to csv
    recom.saveToCSV(outfile)
    print("Recommendations saved as:\n%s\n" % os.path.abspath(outfile))
  except:
  	raise ValueError("Please check arguments.")


if __name__ == '__main__':
  main()