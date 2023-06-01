# test_species_annotation.py
# Testing SpeciesAnnotation class


import libsbml
import numpy as np
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
M_GLUCOSE = 'M_glc__D_e'
D_GLUCOSE = 'D-Glucose'

ONESET_SPECIES_IDS = [M_FDP_C, M_ATP_C]
ONE_CHEBI = 'CHEBI:4167'
ATP_CHEBI = 'CHEBI:30616'
GLUCOSE_CHEBI = 'CHEBI:17634'
ATP_FORMULA = 'C10N5O13P3'
DUMMY_RECOMMENDATION = cn.Recommendation('M_glc__D_e',
                                         # 0.967,
                                         [('CHEBI:4167', 1.0), ('CHEBI:17634', 1.0), ('CHEBI:42758', 1.0)],
                                         ['https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A4167',
                                          'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A17634',
                                          'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A42758'],
                                         ['D-glucopyranose', 'D-glucose', 'aldehydo-D-glucose'])
DUMMY_ID = 'M_glc__D_e'
# Dummy data for calculating accuracy, recalll & precision
DUMMY_REF = {'a': ['ABC', 'BCD'],
              'b': ['DEF']}
DUMMY_PRED = {'a': ['ABC'],
             'b': ['AAA']}
DUMMY_DICT2UPDATE = {M_ATP_C: [ATP_CHEBI]}


#############################
# Tests
#############################
class TestSpeciesAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)

  def testGetCScores(self):
    res = self.spec_cl.getCScores(['hydrogen'])['hydrogen']
    chebis = [k[0] for k in res]
    vals = [k[1] for k in res]
    # checking the list is correctly sorted (max -> min)
    self.assertEqual(vals[0], np.max(vals))
    self.assertEqual(vals[-1], np.min(vals))
    # cheking hydrogen is indeed at the top
    self.assertTrue('CHEBI:18276' in chebis[:5])
    self.assertTrue('CHEBI:49637' in chebis[:5])

  def testGetOneEScore(self):
    res = self.spec_cl.getOneEScore('a', 'ab')
    self.assertEqual(res, 0.5)

  def testGetEScores(self):
    res = self.spec_cl.getEScores(['hydrogen'])['hydrogen']
    chebis = [k[0] for k in res]
    vals = [k[1] for k in res]
    # checking the list is correctly sorted (max -> min)
    self.assertEqual(vals[0], np.max(vals))
    self.assertEqual(vals[-1], np.min(vals))
    # cheking hydrogen is indeed at the top
    self.assertTrue('CHEBI:18276' in chebis[:5])
    self.assertTrue('CHEBI:49637' in chebis[:5])

  # def testPredictAnnotationByEditDistance(self):
  #   one_spec_name = self.spec_cl.model.getSpecies(M_FDP_C).name.lower()
  #   one_pred_spec = self.spec_cl.predictAnnotationByEditDistance(inp_str=one_spec_name)
  #   self.assertTrue(('CHEBI:16905', 1.0) in one_pred_spec[cn.MATCH_SCORE])
  #   self.assertTrue(('CHEBI:49299', 1.0) in one_pred_spec[cn.MATCH_SCORE])
  #   self.assertTrue('CHEBI:16905' in one_pred_spec[cn.CHEBI])
  #   self.assertTrue('CHEBI:49299' in one_pred_spec[cn.CHEBI])
  #   self.assertEqual(one_pred_spec[cn.FORMULA],  ['C6O12P2'])

  def testGetCountOfIndividualCharacters(self):
    one_res = self.spec_cl.getCountOfIndividualCharacters(DUMMY_ID)
    self.assertEqual(one_res['m'], 1)
    self.assertEqual(one_res['g'], 1)
    self.assertEqual(one_res['c'], 1)

  def testPrepareCounterQuery(self):
    one_query, one_name = self.spec_cl.prepareCounterQuery([M_GLUCOSE],
                                                           sa.CHARCOUNT_DF.columns)
    self.assertEqual(one_name[M_GLUCOSE], D_GLUCOSE)
    one_val = np.round(one_query.loc['g', M_GLUCOSE], 2)
    self.assertEqual(one_val, 0.35)


  def testPredictAnnotationByCosineSimilarity(self):
    one_res = self.spec_cl.predictAnnotationByCosineSimilarity(inp_ids=[M_GLUCOSE])
    self.assertEqual(one_res[M_GLUCOSE][cn.NAME_USED], D_GLUCOSE)
    self.assertTrue(GLUCOSE_CHEBI in one_res[M_GLUCOSE][cn.CHEBI])

  def testGetNameToUse(self):
    self.assertEqual(self.spec_cl.getNameToUse(M_GLUCOSE), D_GLUCOSE)

  # def testEvaluatePredictedSpeciesAnnotation(self):
  #   fdp_pred_spec = self.spec_cl.getEScores(inp_strs=[M_FDP_C])[M_FDP_C][:17]
  #   fdp_score = self.spec_cl.evaluatePredictedSpeciesAnnotation(pred=fdp_pred_spec,
  #                                                               name_used=self.spec_cl.getNameToUse(M_FDP_C))
  #   self.assertTrue(fdp_score < 0.906)
  #   self.assertTrue(fdp_score > 0.905)
  #   atp_pred_spec = self.spec_cl.getEScores(inp_strs=[M_ATP_C])[M_ATP_C][:2]
  #   atp_score = self.spec_cl.evaluatePredictedSpeciesAnnotation(pred=atp_pred_spec,
  #                                                               name_used=self.spec_cl.getNameToUse(M_ATP_C))
  #   self.assertTrue(atp_score < 0.973)     
  #   self.assertTrue(atp_score > 0.972)    
  
  def testUpdateSpeciesWithRecommendation(self):
    one_upd = self.spec_cl.updateSpeciesWithRecommendation(DUMMY_RECOMMENDATION)
    self.assertEqual(one_upd, None)
    self.assertTrue((ONE_CHEBI, 1.0) in self.spec_cl.candidates[DUMMY_ID])
    one_formula = cn.REF_CHEBI2FORMULA[ONE_CHEBI]
    self.assertTrue(one_formula in self.spec_cl.formula[DUMMY_ID])


  def testUpdateSpeciesWithDict(self):
    dummy_spec = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)
    dummy_spec.updateSpeciesWithDict(inp_dict=DUMMY_DICT2UPDATE)
    self.assertEqual(dummy_spec.candidates[M_ATP_C], [(ATP_CHEBI, 1.0)])
    self.assertEqual(dummy_spec.formula[M_ATP_C], [ATP_FORMULA])






