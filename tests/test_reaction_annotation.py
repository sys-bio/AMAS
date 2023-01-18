# test_reaction_annotation.py
# Test for ReactionAnnotation class

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
# ID of a reaction
R_PFK = 'R_PFK'
ATP = 'M_atp_c'
COMPONENTS = {'M_fdp_c', 'M_adp_c', 'M_atp_c', 'M_f6p_c', 'M_h_c'}
ONE_CANDIDATE = 'RHEA:12420'
ONE_CHEBI = 'CHEBI:30616'

# Dummy data for calculating accuracy, recalll & precision
DUMMY_REF = {'a': ['ABC', 'BCD'],
              'b': ['DEF']}
DUMMY_PRED = {'a': ['ABC'],
             'b': ['AAA']}


#############################
# Tests
#############################
class TestReactionAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)
    self.reac_cl = ra.ReactionAnnotation(libsbml_fpath = E_COLI_PATH)
    self.pred_species = self.spec_cl.predictAnnotationByCosineSimilarity(inp_ids=list(COMPONENTS))
    self.spec_formula_dict = {val: self.pred_species[val][cn.FORMULA] for val in list(COMPONENTS)}  
    self.pred_reaction = self.reac_cl.predictAnnotation(inp_spec_dict=self.spec_formula_dict,
                                                        inp_reac_list=[R_PFK],
                                                        update=True)

  def testGetReactionComponents(self):
    # When argument is string, i.e., reaction ID
    one_comps = self.reac_cl.getReactionComponents(R_PFK)
    self.assertEqual(COMPONENTS, set(one_comps))
    # When argument is a libsbml.Reaction instance
    one_r = self.reac_cl.model.getReaction(R_PFK)
    two_comps = self.reac_cl.getReactionComponents(one_r)
    self.assertEqual(COMPONENTS, set(two_comps))

  def testPredictAnnotation(self):
    match_scores_dict = self.pred_reaction[cn.MATCH_SCORE]
    self.assertTrue(ONE_CANDIDATE in [val[0] for val in match_scores_dict[R_PFK]])
    self.assertTrue(0.6 in [val[1] for val in match_scores_dict[R_PFK]])

  def testEvaluatePredictedReactionAnnotation(self):
    one_eval = self.reac_cl.evaluatePredictedReactionAnnotation(pred_result=self.pred_reaction)
    self.assertEqual(np.round(one_eval[R_PFK], cn.ROUND_DIGITS), 0.817)


  def testGetRheaElementNum(self):
    num_elements = self.reac_cl.getRheaElementNum(inp_rhea=ONE_CANDIDATE)
    self.assertEqual(num_elements, 5)









