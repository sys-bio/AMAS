# test_reaction_annotation.py
# Test for ReactionAnnotation class

import libsbml
import os
import sys
import unittest

from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra
from AMAS import constants as cn
from AMAS import tools


E_COLI_PATH = os.path.join(os.getcwd(), 'e_coli_core.xml')
BIOMD_248_PATH = os.path.join(os.getcwd(), 'BIOMD0000000248.xml')
# Below ??
# M_FDP_C = 'M_fdp_c'
# M_ATP_C = 'M_atp_c'
# ONESET_SPECIES_IDS = [M_FDP_C, M_ATP_C]

# ID of a reaction
R_PFK = 'R_PFK'
ATP = 'M_atp_c'
COMPONENTS = {'M_fdp_c', 'M_adp_c', 'M_atp_c', 'M_f6p_c', 'M_h_c'}
ONE_CANDIDATE = 'RHEA:12423'
ONE_CHEBI = 'CHEBI:30616'


#############################
# Tests
#############################
class TestReactionAnnotation(unittest.TestCase):

  def setUp(self):
    self.spec_cl = sa.SpeciesAnnotation(libsbml_fpath = E_COLI_PATH)
    self.reac_cl = ra.ReactionAnnotation(libsbml_fpath = E_COLI_PATH)
    pred_species = self.spec_cl.predictAnnotationByName(inp_spec_list=COMPONENTS)
    pred_reaction = self.reac_cl.predictAnnotation(inp_spec_dict=self.spec_cl.formula,
                                                   inp_reac_list=[R_PFK])
  ### Was used for iteration algorithm
  # def testGetMatchScore(self):
  #   one_dict = {'R1': {'M1': 0.7}}
  #   one_match_score = self.reac_cl.getMatchScore(score_dict=one_dict)
  #   self.assertEqual(one_match_score, 0.7)

  def testGetReactionComponents(self):
    # When argument is string, i.e., reaction ID
    one_comps = self.reac_cl.getReactionComponents(R_PFK)
    self.assertEqual(COMPONENTS, set(one_comps))
    # When argument is a libsbml.Reaction instance
    one_r = self.reac_cl.model.getReaction(R_PFK)
    two_comps = self.reac_cl.getReactionComponents(one_r)
    self.assertEqual(COMPONENTS, set(two_comps))

  def testPredictAnnotation(self):
    self.assertTrue(ONE_CANDIDATE in [val[0] for val in self.reac_cl.match_score[R_PFK]])
    self.assertTrue(0.8 in [val[1] for val in self.reac_cl.match_score[R_PFK]])

  def testGetBestOneCandidates(self):
    # When argument is directly given
    one_match_score = {'R1': [('RHEA:1', 1.0), ('RHEA:2', 0.5)]}
    self.assertEqual(self.reac_cl.getBestOneCandidates(one_match_score)['R1'],
                     ['RHEA:1'])
    # When argument is not given
    self.assertEqual(self.reac_cl.getBestOneCandidates()[R_PFK],
                     [ONE_CANDIDATE])

  ### Was used for iteration algorithm
  # def testUpdateSpeciesByAReaction(self):
  #   chebi2update = self.reac_cl.updateSpeciesByAReaction(inp_rid=R_PFK,
  #                                                        inp_spec_dict=self.spec_cl.formula,
  #                                                        inp_rhea=ONE_CANDIDATE)
  #   self.assertEqual(chebi2update[ATP], [ONE_CHEBI])

  def testGetAccuracy(self):
    pred = {R_PFK: ['RHEA:16112']}
    self.assertEqual(self.reac_cl.getAccuracy(pred_annotation=pred), 1.0)







