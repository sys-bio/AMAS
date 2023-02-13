# test_annotation_maker.py
# unittest for AMAS.annotation_maker

# import libsbml
# import numpy as np
import os
# import sys
import unittest


from AMAS import annotation_maker as am
from AMAS import constants as cn
# from AMAS import recommender
# from AMAS import species_annotation as sa
# from AMAS import tools

ONE_CHEBI = 'CHEBI:15414'
ONE_CHEBI_SCORE = 1.0
TWO_CHEBI = 'CHEBI:59789'
TWO_CHEBI_SCORE = 1.0
CANDIDATES = [(ONE_CHEBI, ONE_CHEBI_SCORE),
              (TWO_CHEBI, TWO_CHEBI_SCORE)]
CHEBI = 'chebi'

ONE_ANNOTATION_ITEM = '<rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:15414"/>'
ONE_SCORE_ITEM = '<rdf:li rdf:resource="http://amas/match_score/by_name/0.2"/>'
ONE_TAG = ['<annotation>', '</annotation>']
TWO_TAG = ['<rdf:RDF>', '</rdf:RDF>']

with open(os.path.join(cn.TEST_DIR, 'full_annotation_example.txt')) as file:
    lines = [line.rstrip() for line in file]
FULL_ANNOTATION = '\n'.join(lines)
with open(os.path.join(cn.TEST_DIR, 'annotation_block_example.txt')) as file:
    ONE_BLOCK = [line.rstrip() for line in file]

ONE_INSERTED = ['<annotation>', '  <rdf:RDF>', '  </rdf:RDF>', '</annotation>']
TWO_INSERTED = ['rdf:RDF']

#############################
# Tests
#############################
class TestAnnotationMaker(unittest.TestCase):
  def setUp(self):
    self.cands = CANDIDATES
    # '00001' indicates a meta ID
    self.maker = am.AnnotationMaker('species', '00001')

  def testGetIndent(self):
    one_indent = self.maker.getIndent(num_indents=1)
    self.assertEqual(one_indent, '  ')
    two_indent = self.maker.getIndent()
    self.assertEqual(two_indent, '')

  def testCreateAnnotationBlock(self):
    one_block = self.maker.createAnnotationBlock(identifier=ONE_CHEBI,
                                                 match_score=ONE_CHEBI_SCORE,
                                                 knowledge_resource=None,
                                                 score_type=None,
                                                 nested_prefix='bqbiol:hasProperty')
    still_good = True
    while len(one_block)>0 and len(ONE_BLOCK)>0 and still_good:
      one_item = one_block.pop()
      two_item = ONE_BLOCK.pop()
      self.assertEqual(one_item, two_item)
      if one_item != two_item:
      	still_good = False

  def testCreateAnnotationContainer(self):
  	one_container = self.maker._createAnnotationContainer(['annotation'])
  	self.assertEqual(one_container, ONE_TAG)

  def testCreateAnnotationItem(self):
    one_item = self.maker.createAnnotationItem(CHEBI,
                                               ONE_CHEBI)
    self.assertEqual(one_item, ONE_ANNOTATION_ITEM)

  def testCreateScoreItem(self):
    one_score = self.maker.createScoreItem(0.2,
                                           self.maker.score_type)
    self.assertEqual(one_score, ONE_SCORE_ITEM)

  def testCreateTag(self):
    one_item = self.maker.createAnnotationItem(CHEBI,
                                               ONE_CHEBI)
    self.assertEqual(one_item, ONE_ANNOTATION_ITEM)
    one_tag = self.maker.createTag('annotation', 0)
    self.assertEqual(one_tag, ONE_TAG)

  def testGetAnnotationString(self):
    one_str = self.maker.getAnnotationString(CANDIDATES)
    self.assertEqual(one_str, FULL_ANNOTATION)

  def testInsertEntry(self):
    one_ent = self.maker.insertEntry(inp_str='rdf:RDF',
                                     inp_list=ONE_TAG,
                                     insert_loc=None,
                                     is_prefix=True)
    self.assertEqual(one_ent, ONE_INSERTED)
    two_ent = self.maker.insertEntry(inp_str='rdf:RDF',
                                     inp_list=[],
                                     insert_loc=None,
                                     is_prefix=False)
    self.assertEqual(two_ent, TWO_INSERTED)

  def testInsertList(self):
    one_ins = self.maker.insertList(insert_to=ONE_TAG,
                                    insert_from=TWO_TAG,
                                    start_loc=None)
    self.assertEqual(one_ins, ONE_INSERTED)
    two_ins = self.maker.insertList(insert_to=ONE_TAG,
                                    insert_from=['rdf:RDF'],
                                    start_loc=None)
    self.assertEqual(two_ins,
                     ['<annotation>', '  rdf:RDF', '</annotation>'])













