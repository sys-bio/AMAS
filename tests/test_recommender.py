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
ONE_REAC_CAND = ('RHEA:28830', 1.0)
ONE_REAC_URL = 'https://www.rhea-db.org/rhea/28830'

#############################
# Tests
#############################
class TestRecommender(unittest.TestCase):
  def setUp(self):
    self.recom = recommender.Recommender(libsbml_fpath=BIOMD_190_PATH)

  def testGetSpeciesAnnotation(self):
    species_annotation = self.recom.getSpeciesAnnotation('SAM')
    one_res = species_annotation[0]
    self.assertEqual(one_res.id, 'SAM')
    self.assertEqual(one_res.credibility, 1.0)
    self.assertTrue(ONE_SPEC_CAND in one_res.candidates)
    self.assertTrue(ONE_SPEC_URL in one_res.urls)

  def testGetReactionAnnotation(self):
    reaction_annotation = self.recom.getReactionAnnotation('ODC')
    one_res = reaction_annotation[0]
    self.assertEqual(one_res.id, 'ODC')
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
    self.assertTrue('SAM' in dummy_recom.species.names.keys())
    self.assertTrue('ODC' in dummy_recom.reactions.reaction_components.keys())
    self.assertEqual(len(dummy_recom.species.exist_annotation), 0)
    self.assertEqual(len(dummy_recom.reactions.exist_annotation), 9)


