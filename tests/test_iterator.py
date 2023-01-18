# test_iterator.py
# unittest for iterator class

import copy
import libsbml
import os
import sys
import unittest

from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra
from AMAS import constants as cn
from AMAS import iterator as it
from AMAS import tools


E_COLI_PATH = os.path.join(cn.TEST_DIR, 'e_coli_core.xml')
INIT_SPEC_FORMULA = {'M_accoa_c': ['C23N7O17P3S'],
                     'M_for_c': ['CO2'], 
                     'M_adp_c': ['C115N8O85', 'C21ClN3O2', 'C19N2O2S',
                                 'C17Cl2F3N7O2S', 'C16NO6', 'C19O2',
                                 'C26N7O2S', 'C19O9P', 'C28N6OS',
                                 'C16F3IN2O4', 'C35N4O4', 'C18N4O11',
                                 'C26FIN5O4', 'C20O4', 'C12',
                                 'C10N5O10P2', 'C8NO6', 'C29N6O4S',
                                 'C18ClN2O6S2', 'C27O5', 'C18O2',
                                 'C9N4O5', 'C20N2O5S', 'C101N7O75',
                                 'C14N2O'],
                     'M_pyr_c': ['C3O3'],
                     'M_coa_c': ['C21N7O16P3S'],
                     'M_fdp_c': ['C6O12P2'],
                     'M_f6p_c': ['C6O9P'],
                     'M_atp_c': ['C20O4', 'C18ClN2O6S2'],
                     'M_h_c': ['H']}

REACTION_CANDIDATES = {'R_PFK': [('RHEA:12420', 0.8), ('RHEA:13377', 0.8),
                                 ('RHEA:14213', 0.8), ('RHEA:15653', 0.8),
                                 ('RHEA:16109', 0.8), ('RHEA:20105', 0.8)],
                       'R_PFL': [('RHEA:11844', 1.0), ('RHEA:17425', 0.8),
                                 ('RHEA:22988', 0.8), ('RHEA:22992', 0.8),
                                 ('RHEA:28042', 0.8), ('RHEA:12765', 0.667),
                                 ('RHEA:21912', 0.667), ('RHEA:44140', 0.667)]}

R_PFK = 'R_PFK'
R_PFL = 'R_PFL'
REACTIONS = [R_PFK, R_PFL]


ONE_RHEA = 'RHEA:12420'
ONE_CHEBI = 'CHEBI:15378'
MOLECULE_H = 'H'
TWO_CHEBI = 'CHEBI:58695'
MOLECULE_C6O9P = 'C6O9P'

SPECIES_ATP = 'M_atp_c'
CHEBI_ATP = 'CHEBI:30616'
FORMULA_ATP = 'C10N5O13P3'

ONE_RES_CHEBI = {'M_atp_c': ['CHEBI:30616']}
ONE_SPEC2FORMULA = {'M_atp_c': ['C10N5O13P3']}


#############################
# Tests
#############################
class TestIterator(unittest.TestCase):

  def setUp(self):
    reac_cl = ra.ReactionAnnotation(libsbml_fpath = E_COLI_PATH)
    reac_cl.candidates = REACTION_CANDIDATES
    self.anot_iter = it.Iterator(cur_spec_formula=copy.deepcopy(INIT_SPEC_FORMULA),
                                 reaction_cl=reac_cl,
                                 reactions_to_update=REACTIONS)

  def testGetDictOfRheaComponentFormula(self):
    comp_formula = self.anot_iter.getDictOfRheaComponentFormula(inp_rhea=ONE_RHEA)
    self.assertEqual(comp_formula[ONE_CHEBI], MOLECULE_H)
    self.assertEqual(comp_formula[TWO_CHEBI], MOLECULE_C6O9P)

  def testGetDictMatchByItem(self):
    rhea_comps = self.anot_iter.getDictOfRheaComponentFormula(ONE_RHEA)
    filt_spec_formula = {k:self.anot_iter.orig_spec_formula[k] \
                           for k in self.anot_iter.reactions.reaction_components[R_PFK]}
    upd_spec_chebi = self.anot_iter.getDictMatchByItem(chebi2ref_formula=rhea_comps,
                                                       spec2pred_formula=filt_spec_formula)
    self.assertEqual(upd_spec_chebi, ONE_RES_CHEBI)


  def testGetDictsToUpdate(self):
    spec2chebi, spec2formula = self.anot_iter.getDictsToUpdate(reaction_id=R_PFK)
    self.assertTrue(CHEBI_ATP in spec2chebi[SPECIES_ATP])
    self.assertTrue(FORMULA_ATP in spec2formula[SPECIES_ATP])

  def testGetUpdatedMatchScore(self):
    res = self.anot_iter.getUpdatedMatchScore(cur_spec_formulas=copy.deepcopy(INIT_SPEC_FORMULA),
                                              inp_spec2formula_dict=ONE_SPEC2FORMULA)
    self.assertEqual(res[it.NEW_SCORE], 1.0)
    self.assertEqual(res[it.OLD_SCORE], 0.9)
    self.assertTrue(res[it.INCREASED])

  def testMatch(self):
    res_match = self.anot_iter.match()
    self.assertEqual(res_match, ONE_RES_CHEBI)

  def testRunOneMatchCycle(self):
    match_res = self.anot_iter.runOneMatchCycle()
    self.assertEqual(match_res, ONE_RES_CHEBI)















