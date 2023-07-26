#!/usr/bin/env python

# update_annotation.py
"""
Set annotation of a model file
Usage: python update_annotation.py res.csv files/BIOMD0000000190.xml BIOMD0000000190_upd.xml
"""

import argparse
import itertools
import libsbml
import numpy as np
import os
from os.path import dirname, abspath
import pandas as pd
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import annotation_mak er as am
from AMAS import tools

def main():
  parser = argparse.ArgumentParser(description='Update annotations of a model using user\'s feedback file (.csv)')
  parser.add_argument('infile', type=str, help='path of a model file (.xml) to update annotation')
  parser.add_argument('feedback', type=str, help='path of the file (.csv) containing user\'s feedback')
  parser.add_argument('outfile', type=str, help='file path to save model with updated annotations')
  # csv file with user choice
  args = parser.parse_args()
  user_csv = pd.read_csv(args.feedback)
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
    # going through one id at a time
    for one_id in uids:
      orig_str = ELEMENT_FUNC[one_type](one_id).getAnnotationString()
      df_id = df_type[df_type['id']==one_id]
      dels = list(df_id[df_id[cn.DF_UPDATE_ANNOTATION_COL]=='delete'].loc[:, 'annotation'])
      adds_raw = list(df_id[df_id[cn.DF_UPDATE_ANNOTATION_COL]=='add'].loc[:, 'annotation'])
      # existing annotations to be kept 
      keeps = list(df_id[df_id[cn.DF_UPDATE_ANNOTATION_COL]=='keep'].loc[:, 'annotation'])
      adds = list(set(adds_raw + keeps))
      # if type is 'reaction', need to map rhea terms back to ec/kegg terms to delete them. 
      if one_type == 'reaction':
        rhea_del_terms = list(set(itertools.chain(*[tools.getAssociatedTermsToRhea(val) for val in dels])))
        deled = maker.deleteAnnotation(rhea_del_terms, orig_str)
      elif one_type == 'species':
        deled = maker.deleteAnnotation(dels, orig_str)
      added = maker.addAnnotation(adds, deled, meta_ids[one_id])
      ELEMENT_FUNC[one_type](one_id).setAnnotation(added)
  libsbml.writeSBMLToFile(document, outfile)
  print("...\nUpdated model file saved as:\n%s\n" % os.path.abspath(outfile))


if __name__ == '__main__':
  main()