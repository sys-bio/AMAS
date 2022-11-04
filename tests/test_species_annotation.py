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
M_AMP_C = 'M_amp_c'
ONESET_SPECIES_IDS = [M_FDP_C, M_ATP_C]
ONE_CHEBI = 'CHEBI:15414'
DUMMY_RECOMMENDATION = cn.Recommendation('SAM',
                                         1.0,
                                         [('CHEBI:15414', 1.0), ('CHEBI:59789', 1.0)],
                                         ['https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15414',
                                         'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A59789'])
DUMMY_ID = 'SAM'
# Dummy data for calculating accuracy, recalll & precision
DUMMY_REF = {'a': ['ABC', 'BCD'],
              'b': ['DEF']}
DUMMY_PRED = {'a': ['ABC'],
             'b': ['AAA']}


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
    accuracy1 = self.spec_cl.getAccuracy(ref_annotation=DUMMY_REF,
                                         pred_annotation=DUMMY_PRED)
    self.assertEqual(accuracy1, 0.5)
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_FDP_C])
    accuracy2 = self.spec_cl.getAccuracy()
    self.assertEqual(accuracy2, 1.0)


  def testGetRecall(self):
    recall1 = self.spec_cl.getRecall(ref_annotation=DUMMY_REF,
                                     pred_annotation=DUMMY_PRED,
                                     mean=True)
    self.assertEqual(recall1, 0.25)
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_AMP_C])
    one_res_formula = {M_AMP_C: one_pred_spec[M_AMP_C][cn.FORMULA]}
    recall2 = self.spec_cl.getRecall(pred_annotation=one_res_formula,
                                     mean=True)
    self.assertEqual(recall2, 1.0)


  def testGetPrecision(self):
    precision1 = self.spec_cl.getPrecision(ref_annotation=DUMMY_REF,
                                           pred_annotation=DUMMY_PRED,
                                           mean=True)
    self.assertEqual(precision1, 0.5)
    one_pred_spec = self.spec_cl.predictAnnotationByName(inp_spec_list=[M_AMP_C])
    one_res_formula = {M_AMP_C: one_pred_spec[M_AMP_C][cn.FORMULA]}
    precision2 = self.spec_cl.getPrecision(pred_annotation=one_res_formula,
                                           mean=True)
    self.assertEqual(precision2, 0.1)


  def testGetNameToUse(self):
    self.assertEqual(self.spec_cl.getNameToUse('M_glc__D_e'), 'D-Glucose')

  def testEvaluatePredictedSpeciesAnnotation(self):
    fdp_pred_spec = self.spec_cl.predictAnnotationByEditDistance(inp_str=M_FDP_C)
    fdp_score = self.spec_cl.evaluatePredictedSpeciesAnnotation(pred_result=fdp_pred_spec)
    self.assertTrue(fdp_score < 0.15)
    self.assertTrue(fdp_score > 0.14)
    atp_pred_spec = self.spec_cl.predictAnnotationByEditDistance(inp_str=M_ATP_C)
    atp_score = self.spec_cl.evaluatePredictedSpeciesAnnotation(pred_result=atp_pred_spec)
    self.assertTrue(atp_score < 0.08)  
    self.assertTrue(atp_score > 0.07)    
  
  def testUpdateSpeciesWithRecommendation(self):
    one_upd = self.spec_cl.updateSpeciesWithRecommendation(DUMMY_RECOMMENDATION)
    self.assertEqual(one_upd, None)
    self.assertTrue((ONE_CHEBI, 1.0) in self.spec_cl.candidates[DUMMY_ID])
    one_formula = cn.REF_CHEBI2FORMULA[ONE_CHEBI]
    self.assertTrue(one_formula in self.spec_cl.formula[DUMMY_ID])










