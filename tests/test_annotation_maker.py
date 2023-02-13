# test_annotation_maker.py
# unittest for AMAS.annotation_maker

# import libsbml
# import numpy as np
# import os
# import sys
import unittest


from AMAS import annotation_maker as am
# from AMAS import constants as cn
# from AMAS import recommender
# from AMAS import species_annotation as sa
# from AMAS import tools

ONE_CHEBI = 'CHEBI:15414'
TWO_CHEBI = 'CHEBI:59789'
CANDIDATES = [(ONE_CHEBI, 1.0), (TWO_CHEBI, 1.0)]
CHEBI = 'chebi'

ONE_ANNOTATION_ITEM = '<rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:15414"/>'

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

  def testCreateAnnotationLine(self):
  	one_item = self.maker.createAnnotationItem(CHEBI,
  		                                       ONE_CHEBI)
  	self.assertEqual(one_item, ONE_ANNOTATION_ITEM)










