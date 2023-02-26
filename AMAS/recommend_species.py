# recommend_species.py
"""
Predicts annotations of species using a local XML file
and the species ID. 
Usage: python recommend_species <infile_path> <species_id_1> <species_id_2>.. etc.
"""

import argparse
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='SBML file (.XML) and one or more species IDs in the model.') 
  parser.add_argument('--model', type=str, help='SBML file in the XML format')
  # One or more species IDs can be given
  parser.add_argument('--species', type=str, help='ID of species in the model', nargs='+')
  parser.add_argument('--min_score', type=float, help='Minimum threshold')
  parser.add_argument('--out_dir', type=str, help='Path of directory to save files')
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.model)
  # check if all species are included in the species
  one_fpath = args.model
  specs = args.species
  min_score = args.min_score
  out_dir = args.out_dir
  try:
    recom = recommender.Recommender(libsbml_fpath=one_fpath)
    recom.current_type = 'species'
    res = recom.getSpeciesListRecommendation(specs, get_df=True)
    for idx, one_df in enumerate(res):
      filt_df = recom.autoSelectAnnotation(one_df, min_score)
      recom.updateSelection(specs[idx], filt_df)
    # store file to csv
    os.mkdir(out_dir)
    recom.saveToCSV(os.path.join(out_dir, 'recommendation.csv'))
    recom.saveToSBML(os.path.join(out_dir, 'sbml_model.xml'))
  except:
  	raise ValueError("Please check arguments.")


if __name__ == '__main__':
  main()