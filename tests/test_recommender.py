# test_recommender.py
# unittest for AMAS.recommender

import libsbml
import numpy as np
import os
import sys
import unittest


from AMAS import constants as cn
from AMAS import recommender
from AMAS import species_annotation as sa
from AMAS import tools

BIOMD_190_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000190.xml')
BIOMD_634_PATH = os.path.join(cn.TEST_DIR, 'BIOMD0000000634.xml')
E_COLI_PATH = os.path.join(cn.TEST_DIR, 'e_coli_core.xml')
ONE_SPEC_CAND = ('CHEBI:15414', 1.0)
ONE_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15414'
TWO_SPEC_CAND = ('CHEBI:15729', 1.0)
TWO_SPEC_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A15729'

ONE_REAC_CAND = ('RHEA:28827', 1.0)
ONE_REAC_URL = 'https://www.rhea-db.org/rhea/28827'

SPECIES_SAM = 'SAM'
SPECIES_SAM_NAME = 'S-adenosyl-L-methionine'
SPECIES_ORN = 'ORN'
SPECIES_ATP = 'ATP'
REACTION_ODC = 'ODC'
REACTION_SAMDC = 'SAMdc'
REACTION_SPMS = 'SpmS'
R_PFK = 'R_PFK'
R_PFL = 'R_PFL'
ECOLI_REACTIONS = [R_PFK, R_PFL]
ECOLI_ATP = 'M_atp_c'
ECOLI_RHEA = 'RHEA:12420'

ONE_CHEBI = 'CHEBI:15414'
ATP_CHEBI = 'CHEBI:30616'
FORMULA_ATP = 'C10N5O13P3'


RESULT_RECOM = cn.Recommendation('R_PFK', 0.817,
                                 [('RHEA:12420', 0.6), ('RHEA:13377', 0.6)],
                                 ['https://www.rhea-db.org/rhea/12420', 'https://www.rhea-db.org/rhea/13377'],
                                 ['tagatose-6-phosphate kinase activity', 'phosphoglucokinase activity'])
RESULT_MARKDOWN = '                      R_PFK (credibility score: 0.817)                      \n' + \
                  '+----+--------------+---------------+--------------------------------------+\n' + \
                  '|    | annotation   |   match_score | label                                |\n' + \
                  '+====+==============+===============+======================================+\n' + \
                  '|  1 | RHEA:12420   |         0.600 | tagatose-6-phosphate kinase activity |\n' + \
                  '+----+--------------+---------------+--------------------------------------+\n' + \
                  '|  2 | RHEA:13377   |         0.600 | phosphoglucokinase activity          |\n' + \
                  '+----+--------------+---------------+--------------------------------------+'

#############################
# Tests
#############################
class TestRecommender(unittest.TestCase):
  def setUp(self):
    self.recom = recommender.Recommender(libsbml_fpath=BIOMD_190_PATH)

  def testFilterRecommendationByThreshold(self):
    recom = recommender.Recommender(libsbml_fpath=E_COLI_PATH)
    one_recom = recom.getReactionRecommendation(pred_id=R_PFK)
    two_recom = recom.getReactionRecommendation(pred_id=R_PFL)
    self.assertEqual(None, recom.filterRecommendationByThreshold(inp_recom=one_recom, inp_thresh=0.8))
    filt_two_recom = recom.filterRecommendationByThreshold(inp_recom=two_recom, inp_thresh=0.8)
    self.assertEqual(len(two_recom.candidates), 8)
    self.assertEqual(len(filt_two_recom.candidates), 5)

  def testGetDataFrameFromRecommendation(self):
    df = self.recom.getDataFrameFromRecommendation(rec=RESULT_RECOM,
                                                   show_url=False)
    self.assertEqual(set(df.index), {1,2})


  def testGetMarkdownFromRecommendation(self):
    res = self.recom.getMarkdownFromRecommendation(rec=RESULT_RECOM,
                                                   show_url=False)
    self.assertEqual(res, RESULT_MARKDOWN)

  def testGetSpeciesRecommendation(self):
    one_res = self.recom.getSpeciesRecommendation(pred_id=SPECIES_SAM,
                                                  update=False,
                                                  method='edist')
    self.assertEqual(one_res.id, SPECIES_SAM)
    self.assertEqual(one_res.credibility, 0.975)
    self.assertTrue(ONE_SPEC_CAND in one_res.candidates)
    self.assertTrue(ONE_SPEC_URL in one_res.urls)
    self.assertEqual(self.recom.species.candidates, {})
    self.assertEqual(self.recom.species.formula, {})
    two_res = self.recom.getSpeciesRecommendation(pred_str=SPECIES_SAM_NAME,
                                              update=True,
                                              method='cdist')
    self.assertEqual(two_res.id, SPECIES_SAM_NAME)
    self.assertEqual(two_res.credibility, 0.975)
    self.assertTrue(ONE_SPEC_CAND in two_res.candidates)
    self.assertTrue(ONE_SPEC_URL in two_res.urls)
    self.assertTrue((ONE_CHEBI, 1.0) in self.recom.species.candidates[SPECIES_SAM_NAME])
    one_formula = cn.REF_CHEBI2FORMULA[ONE_CHEBI]
    self.assertTrue(one_formula in self.recom.species.formula[SPECIES_SAM_NAME])    

  def testGetSpeciesIDs(self):
    one_res = self.recom.getSpeciesIDs(pattern="*CoA")
    self.assertTrue('AcCoA' in one_res)
    self.assertTrue('CoA' in one_res)
    self.assertEqual(len(one_res), 2)
    none_res = self.recom.getSpeciesIDs(pattern="AAA")
    self.assertEqual(none_res, None)

  def testGetSpeciesListRecommendation(self):
    specs = self.recom.getSpeciesListRecommendation(pred_ids=[SPECIES_SAM, SPECIES_ORN],
                                                update=False, method='edist')
    one_res = specs[1]
    self.assertEqual(one_res.id, SPECIES_ORN)
    self.assertEqual(one_res.credibility, 0.972)
    self.assertTrue(TWO_SPEC_CAND in one_res.candidates)
    self.assertTrue(TWO_SPEC_URL in one_res.urls)
    self.assertEqual(self.recom.species.candidates, {})
    self.assertEqual(self.recom.species.formula, {})
    two_specs = self.recom.getSpeciesListRecommendation(pred_ids=[SPECIES_SAM, SPECIES_ORN],
                                                    update=True, method='cdist')
    self.assertTrue((ONE_CHEBI, 1.0) in self.recom.species.candidates[SPECIES_SAM])
    one_formula = cn.REF_CHEBI2FORMULA[ONE_CHEBI]
    self.assertTrue(one_formula in self.recom.species.formula[SPECIES_SAM])      

  def testGetReactionRecommendation(self):
    one_res = self.recom.getReactionRecommendation(REACTION_ODC)
    self.assertEqual(one_res.id, REACTION_ODC)
    self.assertEqual(one_res.credibility, 0.817)
    self.assertTrue(ONE_REAC_CAND in one_res.candidates)
    self.assertTrue(ONE_REAC_URL in one_res.urls)

  def testGetReactionIDs(self):
    one_res = self.recom.getReactionIDs(pattern='*CoA', by_species=True)
    self.assertEqual(len(one_res), 4)
    self.assertTrue('SSAT_for_S' in one_res)
    two_res = self.recom.getReactionIDs(pattern='*CoA', by_species=False)
    self.assertEqual(len(two_res), 2)
    self.assertTrue('VCoA' in two_res)

  def testGetReactionListRecommendation(self):
    reacs = self.recom.getReactionListRecommendation(pred_ids=[REACTION_ODC, REACTION_SAMDC])
    one_res = reacs[0]
    self.assertEqual(one_res.id, REACTION_ODC)
    self.assertEqual(one_res.credibility, 0.817)
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

  def testGetSpeciesStatistics(self):
    recom2 = recommender.Recommender(libsbml_fpath=BIOMD_634_PATH)
    spec_stats1 = recom2.getSpeciesStatistics(model_mean=True)
    self.assertEqual(spec_stats1[cn.RECALL], 1.000)
    self.assertEqual(spec_stats1[cn.PRECISION], 0.125)
    spec_stats2 = recom2.getSpeciesStatistics(model_mean=False)
    self.assertEqual(spec_stats2[cn.RECALL][SPECIES_ATP], 1.000)
    self.assertEqual(spec_stats2[cn.PRECISION][SPECIES_ATP], 0.200)

  def testGetReactionStatistics(self):
    reac_stats1 = self.recom.getReactionStatistics(model_mean=True)
    self.assertEqual(reac_stats1[cn.RECALL], 0.694)
    self.assertEqual(reac_stats1[cn.PRECISION], 0.631)
    reac_stats2 = self.recom.getReactionStatistics(model_mean=False)
    self.assertEqual(reac_stats2[cn.RECALL][REACTION_SPMS], 1.000)
    self.assertEqual(reac_stats2[cn.PRECISION][REACTION_SPMS], 0.333)

  def testUpdateAnnotationsByIteration(self):
    recom = recommender.Recommender(libsbml_fpath=E_COLI_PATH)
    _ = recom.getReactionListRecommendation(pred_ids=ECOLI_REACTIONS, spec_method='edist')
    self.assertEqual(recom.species.candidates[ECOLI_ATP][0][0], 'CHEBI:182955')
    self.assertEqual(recom.species.candidates[ECOLI_ATP][0][1], 0.231)
    self.assertTrue('C20O4' in recom.species.formula[ECOLI_ATP])
    self.assertEqual(recom.reactions.candidates[R_PFK][0][0], ECOLI_RHEA)
    self.assertEqual(recom.reactions.candidates[R_PFK][0][1], 0.8)
    recom.updateAnnotationsByIteration()
    self.assertEqual(recom.species.candidates[ECOLI_ATP][0][0], ATP_CHEBI)
    self.assertEqual(recom.species.candidates[ECOLI_ATP][0][1], 0.231)
    self.assertTrue(FORMULA_ATP in recom.species.formula[ECOLI_ATP])
    self.assertEqual(recom.reactions.candidates[R_PFK][0][0], ECOLI_RHEA)
    self.assertEqual(recom.reactions.candidates[R_PFK][0][1], 1.0)







