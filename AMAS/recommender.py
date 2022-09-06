# recommender.py
# Recomender for running annotation predictions


import os
import pickle

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra




class Matcher(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier=cn.RHEA):
  	# creates both species & reaction annotation class instances
  	# 
    pass
    species_annotation = sa.SpeciesAnnotation(libsbml_fpath)
    reaction_anotation = ra.ReactionAnnotataion(libsbml_fpath)
    