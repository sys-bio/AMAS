#!/usr/bin/env python

# update_annotation.py
"""
Set annotation of a model file
Usage: python update_annotation.py res.csv files/BIOMD0000000190.xml BIOMD0000000190_upd.xml
"""

import argparse
import libsbml
import numpy as np
import os
from os.path import dirname, abspath
import pandas as pd
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import annotation_maker as am

def main():
  parser = argparse.ArgumentParser(description='Use user feedback file (.csv) to update model annotation')
  parser.add_argument('infile', type=str, help='Path of model file to update annotation')
  parser.add_argument('csv_select', type=str, help='CSV file with user choice')
  parser.add_argument('outfile', type=str, help='File path to save updated model file')
  # csv file with user choice 
  args = parser.parse_args()
  user_csv = pd.read_csv(args.csv_select)
  chosen = user_csv[user_csv['USE ANNOTATION']==1]
  outfile = args.outfile

  reader = libsbml.SBMLReader()
  document = reader.readSBML(args.infile)
  model = document.getModel()

  element_types = list(np.unique(chosen['type']))
  for one_type in element_types:
    maker = am.AnnotationMaker(one_type)
    df_type = chosen[chosen['type']==one_type]
    uids = list(np.unique(df_type['id']))
    meta_ids = {val:list(df_type[df_type['id']==val]['meta id'])[0] for val in uids}
    for one_id in uids:
      one_annotation = maker.getAnnotationString(list(df_type[df_type['id']==one_id]['annotation']),
                                                 meta_ids[one_id])
      if one_type == 'species':
        model.getSpecies(one_id).setAnnotation(one_annotation)
      elif one_type == 'reaction':
        model.getReaction(one_id).setAnnotation(one_annotation)
  libsbml.writeSBMLToFile(document, outfile)
  print("...\nUpdated model file saved as:\n%s\n" % os.path.abspath(outfile))


if __name__ == '__main__':
  main()