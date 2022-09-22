# constants.py
"""
Constants for modlues
"""

import os

# BASE_DIR = '/Users/woosubs/Desktop/AutomateAnnotation/'
# DATA_DIR = os.path.join(BASE_DIR, "DATA")
# ALGO_DIR = os.path.join(DATA_DIR, "algo")
# CHEBI_DIR = os.path.join(DATA_DIR, "chebi")
# RHEA_DIR = os.path.join(DATA_DIR, "rhea")


# Folder for reference data
# CUR_DIR = os.getcwd()
CUR_DIR = os.path.dirname(os.path.realpath(__file__))
REF_DIR = os.path.join(CUR_DIR, os.pardir, 'files')
TEST_DIR = os.path.join(CUR_DIR, os.pardir, 'tests')

# Strings used in the modules
CHEBI = "chebi"
RHEA = "rhea"
KEGG_REACTION = "kegg.reaction"
MATCH_SCORE = "match_score"
NAME_USED = "name_used"
FORMULA = "formula"

# Default URLs for CHEBI/Rhea
CHEBI_DEFAULT_URL = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A'
RHEA_DEFAULT_URL = 'https://www.rhea-db.org/rhea/'