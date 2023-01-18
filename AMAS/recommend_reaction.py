# recommend_reaction.py
"""
Predicts annotations of reaction(s) using a local XML file
and the reaction ID. 
Usage: recommend_species <filepath> <reaction_id_1> <reaction_id_2>.. etc.
"""

import  argparse
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='SBML file (.XML) and one or more reaction IDs in the model.')
  parser.add_argument('file_path', type=str, help='SBML file in the XML format')
  # One or more species IDs can be given
  parser.add_argument('reaction_id', type=str, help='ID of reactions in the model', nargs='+')
  args = parser.parse_args()
  recom = recommender.Recommender(libsbml_fpath=args.file_path)
  # check if all species are included in the species
  reacts = args.reaction_id
  try:
    res_mkd = recom.getReactionListAnnotation(pred_ids=reacts, get_markdown=True)
    for one_mkd in res_mkd:
      print(one_mkd)
  except:
  	raise ValueError("Please check reaction IDs.")


if __name__ == '__main__':
  main()