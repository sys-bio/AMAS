# test_tools.py


import libsbml
import os
import pickle
import sys
import unittest

from AMAS import constants as cn
from AMAS import tools

with open(os.path.join(cn.CHEBI_DIR, 'chebi_shortened_formula_30apr2022.pickle'), 'rb') as f:
  ref_chebi2formula = pickle.load(f)
BIOMD_248_PATH = os.path.join(os.getcwd(), 'BIOMD0000000248.xml')

#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def setUp(self):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(BIOMD_248_PATH)
    self.model = document.getModel()
    self.annotation = self.model.getReaction('CreatineKinase').getAnnotationString()
  
  def testGetOntologyFromString(self):
    filt_annotation = tools.getOntologyFromString(string_annotation=self.annotation)
    self.assertTrue(('ec-code', '2.7.3.2') in filt_annotation)
    self.assertTrue(('kegg.reaction', 'R01881') in filt_annotation)
    self.assertTrue(('go', 'GO:0004111') in filt_annotation)


  def testGetQualifierFromString(self):
    ec_annotation = tools.getQualifierFromString(input_str=self.annotation,
                                                 qualifier='ec-code')
    self.assertEqual(ec_annotation, ['2.7.3.2'])
    kegg_annotation = tools.getQualifierFromString(input_str=self.annotation,
                                                   qualifier='kegg.reaction')
    self.assertEqual(kegg_annotation, ['R01881'])
    go_annotation = tools.getQualifierFromString(input_str=self.annotation,
                                                 qualifier='go')
    self.assertEqual(go_annotation, ['GO:0004111'])

  def testTransformCHEBIToFormula(self):
    res = tools.transformCHEBIToFormula(inp_list=['CHEBI:18357', 'CHEBI:10'],
                                        ref_to_formula_dict=ref_chebi2formula)
    self.assertTrue('C8NO3' in res)
    self.assertTrue('C36N2O6' in res)

  def testUpdateDictKeyToList(self):
    empty_dict = dict()
    oneitm_dict = {'a': 'x'}
    twoitms_dict = {'b': ['x', 'y']}
    empty_res = tools.updateDictKeyToList(inp_orig_dict=empty_dict,
                                          inp_new_dict=empty_dict)
    self.assertEqual(empty_res, dict())
    oneitm_res = tools.updateDictKeyToList(inp_orig_dict=empty_dict,
                                           inp_new_dict=oneitm_dict)
    self.assertEqual(oneitm_res, {'a': ['x']})
    twoitm_res = tools.updateDictKeyToList(inp_orig_dict=oneitm_res,
                                           inp_new_dict=twoitms_dict)
    self.assertEqual(set(twoitm_res.keys()), {'a','b'})
    self.assertEqual(twoitm_res['b'], ['x', 'y'])





