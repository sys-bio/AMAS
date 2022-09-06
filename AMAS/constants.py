# constants.py
"""
Constants for modlues
"""

import os

BASE_DIR = '/Users/woosubs/Desktop/AutomateAnnotation/'
DATA_DIR = os.path.join(BASE_DIR, "DATA")
ALGO_DIR = os.path.join(DATA_DIR, "algo")
CHEBI_DIR = os.path.join(DATA_DIR, "chebi")
RHEA_DIR = os.path.join(DATA_DIR, "rhea")


# Folder for reference data
CUR_DIR = os.getcwd()
REF_DIR = os.path.join(CUR_DIR, os.pardir, 'data')

# Strings used in the modules
CHEBI = "chebi"
RHEA = "rhea"
KEGG_REACTION = "kegg.reaction"
MATCH_SCORE = "match_score"
FORMULA = "formula"