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

  def testGetReactionComponents(self):
    # When argument is string, i.e., reaction ID
    one_comps = self.reac_cl.getReactionComponents(R_PFK)
    self.assertEqual(COMPONENTS, set(one_comps))
    # When argument is a libsbml.Reaction instance
    one_r = self.reac_cl.model.getReaction(R_PFK)
    two_comps = self.reac_cl.getReactionComponents(one_r)
    self.assertEqual(COMPONENTS, set(two_comps))

  def testGetRScores(self):
    specs = {'M_f6p_c': ['C6O9P'],
             'M_fdp_c': ['C6O12P2'],
             'M_atp_c': ['C30N4O29P3'],
             'M_h_c': ['H', 'C6N3O2', 'C6N3O', '[3He]', 'C12N6O3'],
             'M_adp_c': ['C30O8P', 'C39O8P']}
    res = self.reac_cl.getRScores(spec_dict=specs,
                                  reacs=[R_PFK],
                                  mssc='top',
                                  cutoff=0.0)[R_PFK]
    res_vals = [val[1] for val in res]
    self.assertEqual(res_vals[0], np.max(res_vals))
    self.assertEqual(res_vals[-1], np.min(res_vals))

  def testGetRheaElementNum(self):
    num_elements = self.reac_cl.getRheaElementNum(inp_rhea=ONE_CANDIDATE)
    self.assertEqual(num_elements, 5)









