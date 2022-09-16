# test_recommender.py
# unittest for AMAS.recommender

import libsbml
import os
import sys
import unittest


from AMAS import constants as cn
from AMAS import recommender
from AMAS import tools

BIOMD_190_PATH = os.path.join(os.getcwd(), 'BIOMD0000000190.xml')

#############################
# Tests
#############################
class TestRecommender(unittest.TestCase):
  def setUp(self):
    self.recom = recommender.Recommender(libsbml_fpath=BIOMD_190_PATH)