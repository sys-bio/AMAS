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
  parser = argparse.ArgumentParser(description='Update annotations of a model using user\'s feedback file (.csv)')
  parser.add_argument('infile', type=str, help='path of a model file (.xml) to update annotation')
  parser.add_argument('csv_select', type=str, help='feedback file (.csv) with user choice')
  parser.add_argument('outfile', type=str, help='file path to save model with updated annotations')
  # csv file with user choice
  args = parser.parse_args()
  user_csv = pd.read_csv(args.csv_select)
  # Only takes cells with values 'add' or 'delete'
  chosen = user_csv[(user_csv['UPDATE ANNOTATION']=='add') |\
                   (user_csv['UPDATE ANNOTATION']=='delete')]
  outfile = args.outfile

  reader = libsbml.SBMLReader()
  document = reader.readSBML(args.infile)
  model = document.getModel()
  ELEMENT_FUNC = {'species': model.getSpecies,
                  'reaction': model.getReaction}

  element_types = list(np.unique(chosen['type']))
  for one_type in element_types:
    maker = am.AnnotationMaker(one_type)
    ACTION_FUNC = {'delete': maker.deleteAnnotation,
                   'add': maker.addAnnotation}
    df_type = chosen[chosen['type']==one_type]
    uids = list(np.unique(df_type['id']))
    meta_ids = {val:list(df_type[df_type['id']==val]['meta id'])[0] for val in uids}
    for one_id in uids:
      orig_str = ELEMENT_FUNC[one_type](one_id).getAnnotationString()
      df_id = df_type[df_type['id']==one_id]
      dels = list(df_id[df_id[cn.DF_UPDATE_ANNOTATION_COL]=='delete'].loc[:, 'annotation'])
      adds = list(df_id[df_id[cn.DF_UPDATE_ANNOTATION_COL]=='add'].loc[:, 'annotation'])
      deled = maker.deleteAnnotation(dels, orig_str)
      added = maker.addAnnotation(adds, deled, meta_ids[one_id])
      ELEMENT_FUNC[one_type](one_id).setAnnotation(added)
  libsbml.writeSBMLToFile(document, outfile)
  print("...\nUpdated model file saved as:\n%s\n" % os.path.abspath(outfile))


if __name__ == '__main__':
  main()