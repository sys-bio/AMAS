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
    res = self.spec_cl.getCScores(inp_strs=['hydrogen'],
                                  mssc='top',
                                  cutoff=0.0)['hydrogen']
    chebis = [k[0] for k in res]
    vals = [k[1] for k in res]
    # checking the list is correctly sorted (max -> min)
    self.assertTrue(abs(vals[0]-np.max(vals)) < cn.TOLERANCE)
    self.assertTrue(abs(vals[-1]-np.max(vals)) < cn.TOLERANCE)
    # cheking hydrogen is indeed at the top
    self.assertTrue('CHEBI:18276' in chebis[:5])
    self.assertTrue('CHEBI:49637' in chebis[:5])

  def testGetOneEScore(self):
    res = self.spec_cl.getOneEScore('a', 'ab')
    self.assertEqual(res, 0.5)

  def testGetEScores(self):
    res = self.spec_cl.getEScores(inp_strs=['hydrogen'],
                                  mssc='top',
                                  cutoff=0.0)['hydrogen']
    chebis = [k[0] for k in res]
    vals = [k[1] for k in res]
    # checking the list is correctly sorted (max -> min)
    self.assertTrue(abs(vals[0]-np.max(vals)) < cn.TOLERANCE)
    self.assertTrue(abs(vals[-1]-np.max(vals)) < cn.TOLERANCE)
    # cheking hydrogen is indeed at the top
    self.assertTrue('CHEBI:18276' in chebis[:5])
    self.assertTrue('CHEBI:49637' in chebis[:5])

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

  def testGetNameToUse(self):
    self.assertEqual(self.spec_cl.getNameToUse(M_GLUCOSE), D_GLUCOSE)  
  
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






