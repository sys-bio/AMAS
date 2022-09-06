# matcher.py
# Matcher class for running annotation predictions


import os
import pickle

from annotation_recommender import constants as cn
from annotation_recommender import tools
from annotation_recommender import species_annotation as sa
from annotation_recommender import reaction_annotation as ra




class Matcher(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier=cn.RHEA):
  	# creates both species & reaction annotation class instances
  	# 
    pass
    species_annotation = sa.SpeciesAnnotation(libsbml_fpath)
    reaction_anotation = ra.ReactionAnnotataion(libsbml_fpath)
    