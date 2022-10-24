# test_recommender.py
# unittest for AMAS.recommender

import libsbml
import os
import sys
import unittest


from AMAS import constants as cn
from AMAS import recommender
from AMAS import species_annotation as sa
from AMAS import tools

BIOMD_190_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000190.xml')
ONE_SPEC_CAND = ('CHEBI:15414', 1.0)
ONE_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15414'
TWO_SPEC_CAND = ('CHEBI:15729', 1.0)
TWO_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15729'

ONE_REAC_CAND = ('RHEA:28830', 1.0)
ONE_REAC_URL = 'https://www.rhea-db.org/rhea/28830'

SPECIES_SAM = 'SAM'
SPECIES_SAM_NAME = 'S-adenosyl-L-methionine'
SPECIES_ORN = 'ORN'
REACTION_ODC = 'ODC'
REACTION_SAMDC = 'SAMdc'

ONE_CHEBI = 'CHEBI:15414'


#############################
# Tests
#############################
class TestRecommender(unittest.TestCase):
  def setUp(self):
    self.recom = recommender.Recommender(libsbml_fpath=BIOMD_190_PATH)

  def testGetSpeciesAnnotation(self):
    one_res = self.recom.getSpeciesAnnotation(pred_id=SPECIES_SAM, update=False)
    # one_res = species_annotation
    self.assertEqual(one_res.id, SPECIES_SAM)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(ONE_SPEC_CAND in one_res.candidates)
    self.assertTrue(ONE_SPEC_URL in one_res.urls)
    self.assertEqual(self.recom.species.candidates, {})
    self.assertEqual(self.recom.species.formula, {})
    two_res = self.recom.getSpeciesAnnotation(pred_str=SPECIES_SAM_NAME, update=True)
    self.assertEqual(two_res.id, SPECIES_SAM_NAME)
    self.assertEqual(two_res.credibility, 1.0)
    self.assertTrue(ONE_SPEC_CAND in two_res.candidates)
    self.assertTrue(ONE_SPEC_URL in two_res.urls)
    self.assertTrue((ONE_CHEBI, 1.0) in self.recom.species.candidates[SPECIES_SAM_NAME])
    one_formula = cn.ref_chebi2formula[ONE_CHEBI]
    self.assertTrue(one_formula in self.recom.species.formula[SPECIES_SAM_NAME])    

  def testGetSpeciesListAnnotation(self):
    specs = self.recom.getSpeciesListAnnotation(pred_ids=[SPECIES_SAM, SPECIES_ORN],
                                                update=False)
    one_res = specs[1]
    self.assertEqual(one_res.id, SPECIES_ORN)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(TWO_SPEC_CAND in one_res.candidates)
    self.assertTrue(TWO_SPEC_URL in one_res.urls)
    self.assertEqual(self.recom.species.candidates, {})
    self.assertEqual(self.recom.species.formula, {})
    two_specs = self.recom.getSpeciesListAnnotation(pred_ids=[SPECIES_SAM, SPECIES_ORN],
                                                    update=True)
    self.assertTrue((ONE_CHEBI, 1.0) in self.recom.species.candidates[SPECIES_SAM])
    one_formula = cn.ref_chebi2formula[ONE_CHEBI]
    self.assertTrue(one_formula in self.recom.species.formula[SPECIES_SAM])      

  def testGetReactionAnnotation(self):
    one_res = self.recom.getReactionAnnotation(REACTION_ODC)
    self.assertEqual(one_res.id, REACTION_ODC)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(ONE_REAC_CAND in one_res.candidates)
    self.assertTrue(ONE_REAC_URL in one_res.urls)

  def testGetReactionListAnnotation(self):
    reacs = self.recom.getReactionListAnnotation(pred_list=[REACTION_ODC, REACTION_SAMDC])
    one_res = reacs[0]
    self.assertEqual(one_res.id, REACTION_ODC)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(ONE_REAC_CAND in one_res.candidates)
    self.assertTrue(ONE_REAC_URL in one_res.urls)

  def testParseSBML(self):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(BIOMD_190_PATH)
    dummy_recom = recommender.Recommender(document)
    # checking if model was loaded successfully
    self.assertEqual(len(dummy_recom.species.names), 11)
    self.assertEqual(len(dummy_recom.reactions.reaction_components), 13)
    self.assertTrue(SPECIES_SAM in dummy_recom.species.names.keys())
    self.assertTrue(REACTION_ODC in dummy_recom.reactions.reaction_components.keys())
    self.assertEqual(len(dummy_recom.species.exist_annotation), 0)
    self.assertEqual(len(dummy_recom.reactions.exist_annotation), 9)


