# test_species_annotation.py
# Testing SpeciesAnnotation class


import libsbml
import os
import sys
import unittest

from annotation_recommender import species_annotation as sa
from annotation_recommender import reaction_annotation as ra
from annotation_recommender import constants as cn
from annotation_recommender import tools


E_COLI_PATH = os.path.join(os.getcwd(), 'e_coli_core.xml')
BIOMD_248_PATH = os.path.join(os.getcwd(), 'BIOMD0000000248.xml')
# IDs of species
M_FDP_C = 'M_fdp_c'
M_ATP_C = 'M_atp_c'
ONESET_SPECIES_IDS = [M_FDP_C, M_ATP_C]


#############################
# Tests
#############################
class TestSpeciesAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)
  
  def testPredictAnnotationByEditDistance(self):
    one_spec_name = self.spec_cl.model.getSpecies(M_FDP_C).name.lower()
    one_pred_spec = self.spec_cl.predictAnnotationByEditDistance(inp_str=one_spec_name)
    self.assertEqual(one_pred_spec[cn.MATCH_SCORE], 1.0)
    self.assertTrue('CHEBI:16905' in one_pred_spec[cn.CHEBI])
    self.assertTrue('CHEBI:49299' in one_pred_spec[cn.CHEBI])
    self.assertEqual(one_pred_spec[cn.FORMULA],  ['C6O12P2'])

  def testPredictAnnotationByName(self):
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_FDP_C])
    self.assertEqual(one_pred_spec[M_FDP_C][cn.MATCH_SCORE], 1.0)
    self.assertTrue('CHEBI:16905' in one_pred_spec[M_FDP_C][cn.CHEBI])
    self.assertTrue('CHEBI:49299' in one_pred_spec[M_FDP_C][cn.CHEBI])
    self.assertEqual(one_pred_spec[M_FDP_C][cn.FORMULA],  ['C6O12P2'])


  def testGetAccuracy(self):
    # test with a dummy name
    dummy_ref1 = {'a': ['ABC', 'BCD'],
                  'b': ['DEF']}
    pref_ref1 = {'a': ['ABC'],
                 'b': ['AAA']}
    accuracy1 = self.spec_cl.getAccuracy(ref_annotation=dummy_ref1,
                                            pred_annotation=pref_ref1)
    self.assertEqual(accuracy1, 0.5)
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_FDP_C])
    accuracy2 = self.spec_cl.getAccuracy()
    self.assertEqual(accuracy2, 1.0)
    