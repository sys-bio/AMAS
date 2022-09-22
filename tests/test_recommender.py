# test_recommender.py
# unittest for AMAS.recommender

import libsbml
import os
import sys
import unittest


from AMAS import constants as cn
from AMAS import recommender
from AMAS import tools

BIOMD_190_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000190.xml')
ONE_SPEC_CAND = ('CHEBI:15414', 1.0)
ONE_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15414'
TWO_SPEC_CAND = ('CHEBI:15729', 1.0)
TWO_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15729'

ONE_REAC_CAND = ('RHEA:28830', 1.0)
ONE_REAC_URL = 'https://www.rhea-db.org/rhea/28830'

SPECIES_SAM = 'SAM'
SPECIES_ORN = 'ORN'
REACTION_ODC = 'ODC'
REACTION_SAMDC = 'SAMdc'


#############################
# Tests
#############################
class TestRecommender(unittest.TestCase):
  def setUp(self):
    self.recom = recommender.Recommender(libsbml_fpath=BIOMD_190_PATH)

  def testGetSpeciesAnnotation(self):
    species_annotation = self.recom.getSpeciesAnnotation(pred_id=SPECIES_SAM)
    one_res = species_annotation
    self.assertEqual(one_res.id, SPECIES_SAM)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(ONE_SPEC_CAND in one_res.candidates)
    self.assertTrue(ONE_SPEC_URL in one_res.urls)

  def testGetSpeciesListAnnotation(self):
    specs = self.recom.getSpeciesListAnnotation(pred_list=[SPECIES_SAM, SPECIES_ORN],
                                                id=True)
    one_res = specs[1]
    self.assertEqual(one_res.id, SPECIES_ORN)
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(TWO_SPEC_CAND in one_res.candidates)
    self.assertTrue(TWO_SPEC_URL in one_res.urls)

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


