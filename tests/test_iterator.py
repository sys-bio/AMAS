# test_iterator.py
# unittest for iterator class

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
R_PFK = 'R_PFK'
R_PFL = 'R_PFL'
reactions = [R_PFK, R_PFL]
ONE_RHEA = 'RHEA:12423'
ONE_CHEBI = 'CHEBI:15378'
MOLECULE_H = 'H'


#############################
# Tests
#############################
class TestIterator(unittest.TestCase):

  def setUp(self):
    reac_cl = ra.ReactionAnnotation(libsbml_fpath = E_COLI_PATH)
    self.anot_iter = it.Iterator(cur_spec_formula=INIT_SPEC_FORMULA,
                                 reaction_cl=reac_cl,
                                 reactions_to_update = reactions)

  def testGetDictOfRheaComponentFormula(self):
    comp_formula = self.anot_iter.getDictOfRheaComponentFormula(inp_rhea=ONE_RHEA)
    self.assertEqual(comp_formula[ONE_CHEBI], MOLECULE_H)



















