# constants.py
"""
Constants for modlues
"""

import collections
import compress_pickle
import os

# Folder for reference data
CUR_DIR = os.path.dirname(os.path.realpath(__file__))
REF_DIR = os.path.join(CUR_DIR, 'files')
TEST_DIR = os.path.join(CUR_DIR, os.pardir, 'tests')

# Strings used in the modules
CANDIDATES = 'candidates'
CHEBI = "chebi"
RHEA = "rhea"
KEGG_REACTION = "kegg.reaction"
MATCH_SCORE = "match_score"
NAME_USED = "name_used"
FORMULA = "formula"
QUERY_DF = 'query_df'

# Default URLs for CHEBI/Rhea
CHEBI_DEFAULT_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A'
RHEA_DEFAULT_URL = 'https://www.rhea-db.org/rhea/'

# Output; namedtuple 'Recommendation'
Recommendation = collections.namedtuple('Recommendation',
                                        ['id', 'credibility', 'candidates', 'urls'])

with open(os.path.join(REF_DIR, 'chebi_shortened_formula_comp.lzma'), 'rb') as f:
  ref_chebi2formula = compress_pickle.load(f)