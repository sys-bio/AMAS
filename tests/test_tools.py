# test_tools.py


import libsbml
import os
import compress_pickle
import sys
import unittest

from AMAS import constants as cn
from AMAS import tools

with open(os.path.join(cn.REF_DIR, 'chebi_shortened_formula_comp.lzma'), 'rb') as f:
  ref_shortened_chebi_to_formula = compress_pickle.load(f)
BIOMD_248_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000248.xml')

# dictionary to test getRecall and getPrecision
DUMMY_ID = 'SAM'
# Dummy data for calculating accuracy, recalll & precision
DUMMY_REF = {'a': ['ABC', 'BCD'],
              'b': ['DEF']}
DUMMY_PRED = {'a': ['ABC'],
             'b': ['AAA']}
DUMMY_A = 'a'
DUMMY_B = 'b'

ATP = 'ATP'
ATP_CHEBI = ['CHEBI:15422']
REACTION_CREATINEKINASE = 'CreatineKinase'
CREATINEKINASE_ANNOTATION = ['RHEA:17160']

#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def setUp(self):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(BIOMD_248_PATH)
    self.model = document.getModel()
    self.annotation = self.model.getReaction(REACTION_CREATINEKINASE).getAnnotationString()

  def testExtractExistingSpeciesAnnotation(self):
    spec_annotation = tools.extractExistingSpeciesAnnotation(inp_model=self.model)
    self.assertTrue(spec_annotation[ATP], ATP_CHEBI)

  def testExtractExistingReactionAnnotation(self):
    one_dict_rhea = tools.extractExistingReactionAnnotation(inp_model=self.model)
    self.assertEqual(one_dict_rhea[REACTION_CREATINEKINASE], CREATINEKINASE_ANNOTATION)

  def textExtractRheaFromAnnotationString(self):
    one_str = model.getReaction(REACTION_CREATINEKINASE).getAnnotationString()
    one_list_rhea = tools.extractExistingRheaFromAnnotationString(inp_str=one_str)
    self.assertEqual(one_list_rhea, CREATINEKINASE_ANNOTATION)

  def testGetOntologyFromString(self):
    filt_annotation = tools.getOntologyFromString(string_annotation=self.annotation)
    self.assertTrue((cn.EC, '2.7.3.2') in filt_annotation)
    self.assertTrue((cn.KEGG_REACTION, 'R01881') in filt_annotation)
    self.assertTrue((cn.GO, 'GO:0004111') in filt_annotation)

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
                                        ref_to_formula_dict=ref_shortened_chebi_to_formula)
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

  def testGetRecall(self):
    recall1 = tools.getRecall(ref=DUMMY_REF,
                              pred=DUMMY_PRED,
                              mean=True)
    self.assertEqual(recall1, 0.25)
    recall2 = tools.getRecall(ref=DUMMY_REF,
                              pred=DUMMY_PRED,
                              mean=False)
    self.assertEqual(recall2[DUMMY_A], 0.5)
    self.assertEqual(recall2[DUMMY_B], 0.0)


  def testGetPrecision(self):
    precision1 = tools.getPrecision(ref=DUMMY_REF,
                                    pred=DUMMY_PRED,
                                    mean=True)
    self.assertEqual(precision1, 0.5)
    precision2 = tools.getPrecision(ref=DUMMY_REF,
                                    pred=DUMMY_PRED,
                                    mean=False)
    self.assertEqual(precision2[DUMMY_A], 1.0)
    self.assertEqual(precision2[DUMMY_B], 0.0)




