# test_species_annotation.py
# Testing SpeciesAnnotation class


import libsbml
import os
import sys
import unittest

from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra
from AMAS import constants as cn
from AMAS import tools


E_COLI_PATH = os.path.join(cn.TEST_DIR, 'e_coli_core.xml')
BIOMD_248_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000248.xml')
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
    self.assertTrue(('CHEBI:16905', 1.0) in one_pred_spec[cn.MATCH_SCORE])
    self.assertTrue(('CHEBI:49299', 1.0) in one_pred_spec[cn.MATCH_SCORE])
    self.assertTrue('CHEBI:16905' in one_pred_spec[cn.CHEBI])
    self.assertTrue('CHEBI:49299' in one_pred_spec[cn.CHEBI])
    self.assertEqual(one_pred_spec[cn.FORMULA],  ['C6O12P2'])

  def testPredictAnnotationByName(self):
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_FDP_C])
    self.assertTrue(('CHEBI:16905', 1.0) in one_pred_spec[M_FDP_C][cn.MATCH_SCORE])
    self.assertTrue(('CHEBI:49299', 1.0) in one_pred_spec[M_FDP_C][cn.MATCH_SCORE])
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

  def testGetNameToUse(self):
    self.assertEqual(self.spec_cl.getNameToUse('M_glc__D_e'), 'D-Glucose')

  def testEvaluatePredictedSpeciesAnnotation(self):
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=ONESET_SPECIES_IDS)
    scores = self.spec_cl.evaluatePredictedSpeciesAnnotation(ONESET_SPECIES_IDS)
    self.assertEqual(scores[M_FDP_C], 1.0)
    self.assertTrue(scores[M_ATP_C] < 0.04)  
    self.assertTrue(scores[M_ATP_C] > 0.038)    
  












