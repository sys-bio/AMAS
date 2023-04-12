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
CHEBI = 'chebi'
OBO_CHEBI = 'obo.chebi'
EC = 'ec-code'
EC_HEADER = 'EC:'
KEGG_HEADER = 'KEGG:'
GO = 'go'
RHEA = "rhea"
RHEA_HEADER = 'RHEA:'
KEGG_REACTION = "kegg.reaction"
MATCH_SCORE = "match_score"
NAME_USED = "name_used"
FORMULA = "formula"
QUERY_DF = 'query_df'
RECALL = 'recall'
PRECISION = 'precision'

# For resulting DataFrame
DF_MATCH_SCORE_COL = 'match score'
DF_UPDATE_ANNOTATION_COL = 'UPDATE ANNOTATION'

# Tolerance to determine identical numerical values
TOLERANCE = 0.00001
# Digits to be rounded up
ROUND_DIGITS = 3

# Default URLs for CHEBI/Rhea
CHEBI_DEFAULT_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A'
RHEA_DEFAULT_URL = 'https://www.rhea-db.org/rhea/'

# Output; namedtuple 'Recommendation'
Recommendation = collections.namedtuple('Recommendation',
                                        ['id', 'credibility', 'candidates', 'urls', 'labels'])

with open(os.path.join(REF_DIR, 'chebi_shortened_formula_comp.lzma'), 'rb') as f:
  REF_CHEBI2FORMULA = compress_pickle.load(f)
with open(os.path.join(REF_DIR, 'chebi2label.lzma'), 'rb') as f:
  REF_CHEBI2LABEL = compress_pickle.load(f)
with open(os.path.join(REF_DIR, 'ec2mrhea.lzma'), 'rb') as handle:
  REF_EC2RHEA = compress_pickle.load(handle)
with open(os.path.join(REF_DIR, 'kegg2mrhea.lzma'), 'rb') as handle:
  REF_KEGG2RHEA = compress_pickle.load(handle)
with open(os.path.join(REF_DIR, 'rhea_all2master.lzma'), 'rb') as f:
  REF_RHEA2MASTER = compress_pickle.load(f)
with open(os.path.join(REF_DIR, 'mrhea2chebi_prime.lzma'), 'rb') as f:
  REF_RHEA2CHEBI = compress_pickle.load(f)
with open(os.path.join(REF_DIR, 'rhea2label.lzma'), 'rb') as f:
  REF_RHEA2LABEL = compress_pickle.load(f)  