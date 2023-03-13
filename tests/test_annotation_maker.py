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
ONE_CHEBI_SCORE = 1.00
TWO_CHEBI = 'CHEBI:59789'
TWO_CHEBI_SCORE = 1.00
CANDIDATES = [ONE_CHEBI, TWO_CHEBI]
CHEBI = 'chebi'

ONE_ANNOTATION_ITEM = '<rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:15414"/>'
ONE_CONTAINER = ['<annotation>', '</annotation>']
TWO_CONTAINER = ['<rdf:RDF>', '</rdf:RDF>']

with open(os.path.join(cn.TEST_DIR, 'full_annotation_example.txt')) as file:
    lines = [line.rstrip() for line in file]
FULL_ANNOTATION = '\n'.join(lines)

with open(os.path.join(cn.TEST_DIR, 'existing_annotation_example.txt')) as file:
    lines = [line.rstrip() for line in file]
EXISTING_ANNOTATION = '\n'.join(lines)

with open(os.path.join(cn.TEST_DIR, 'augmented_annotation_example.txt')) as file:
    lines = [line.rstrip() for line in file]
AUGMENTED_ANNOTATION = '\n'.join(lines)
CHEBI_TO_AUGMENT = ['CHEBI:15414', 'CHEBI:00000']

ONE_INSERTED = ['<annotation>', '  <rdf:RDF>', '  </rdf:RDF>', '</annotation>']
TWO_INSERTED = ['rdf:RDF']

#############################
# Tests
#############################
class TestAnnotationMaker(unittest.TestCase):
  def setUp(self):
    self.cands = CANDIDATES
    self.maker = am.AnnotationMaker('species')

  def testGetIndent(self):
    one_indent = self.maker.getIndent(num_indents=1)
    self.assertEqual(one_indent, '  ')
    two_indent = self.maker.getIndent()
    self.assertEqual(two_indent, '')

  def testCreateAnnotationContainer(self):
  	one_container = self.maker.createAnnotationContainer(['annotation'])
  	self.assertEqual(one_container, ONE_CONTAINER)

  def testCreateAnnotationItem(self):
    one_item = self.maker.createAnnotationItem(CHEBI,
                                               ONE_CHEBI)
    self.assertEqual(one_item, ONE_ANNOTATION_ITEM)

  def testCreateTag(self):
    one_item = self.maker.createAnnotationItem(CHEBI,
                                               ONE_CHEBI)
    self.assertEqual(one_item, ONE_ANNOTATION_ITEM)
    one_tag = self.maker.createTag('annotation')
    self.assertEqual(one_tag, ONE_CONTAINER)

  def testGetAnnotationString(self):
    one_str = self.maker.getAnnotationString(CANDIDATES, meta_id='metaid_00001')
    self.assertEqual(one_str, FULL_ANNOTATION)

  def testInsertList(self):
    one_ins = self.maker.insertList(insert_to=ONE_CONTAINER,
                                    insert_from=TWO_CONTAINER,
                                    start_loc=None)
    self.assertEqual(one_ins, ONE_INSERTED)
    two_ins = self.maker.insertList(insert_to=ONE_CONTAINER,
                                    insert_from=['rdf:RDF'],
                                    start_loc=None)
    self.assertEqual(two_ins,
                     ['<annotation>', '  rdf:RDF', '</annotation>'])

  def testDivideExistingAnnotation(self):
    res_dict = self.maker.divideExistingAnnotation(EXISTING_ANNOTATION)
    container = res_dict['container']
    items = res_dict['items']
    self.assertEqual(container[0], '<annotation>')
    self.assertEqual(items[0], '<rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:15414"/>')

  def testAugmentAnnotation(self):
    augmented_annotation = self.maker.augmentAnnotation(CHEBI_TO_AUGMENT, EXISTING_ANNOTATION)
    self.assertEqual(augmented_annotation, AUGMENTED_ANNOTATION)

