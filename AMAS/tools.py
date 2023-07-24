# tools.py

from AMAS import constants as cn

import itertools
import numpy as np
import re

def applyMSSC(pred,
              mssc,
              cutoff):
  """
  Apply MSSC to a predicted results. 
  If cutoff is too high, 
  return an empty list.

  Parameters
  ----------
  pred: list-tuple
      [(CHEBI:XXXXX, 1.0), etc.]
  mssc: string
  cutoff: float

  Returns
  -------
  filt: list-tuple
      [(CHEBI:XXXXX, 1.0), etc.]
  """
  filt_pred = [val for val in pred if val[1]>=cutoff]
  if not filt_pred:
    return []
  if mssc == 'top':
    max_val = np.max([val[1] for val in filt_pred])
    res_pred = [val for val in filt_pred if val[1]==max_val]
  elif mssc == 'above':
    res_pred = filt_pred
  return res_pred

def extractExistingSpeciesAnnotation(inp_model, qualifier=cn.CHEBI):
  """
  Get existing annotation of species
  that contains ChEBI terms

  Parameters
  ---------
  qualifier: str
      'chebi' or 'obo.chebi'?
  """
  exist_raw = {val.getId():getQualifierFromString(val.getAnnotationString(), [cn.CHEBI, cn.OBO_CHEBI]) \
               for val in inp_model.getListOfSpecies()}
  exist_filt = {val:exist_raw[val] for val in exist_raw.keys() \
                if exist_raw[val]}
  return exist_filt

def extractExistingReactionAnnotation(inp_model):
  """
  Get existing annotation of reactions in Rhea
  in Bi-directional format (RHEA:10003, etc.)
  This will extract annotations from three
  knowledge resources:
  1. Rhea (mapped to BI-format)
  2. KEGG (kegg.reaction mapped to Rhea-BI)
  3. EC-Number (or ec-code, mapped to list of Rhea-BIs)
  
  Once they are mapped to a list of Rhea terms,
  a list of unique Rhea-Bi terms will be filtered
  and be assigned to exist_annotation of the 
  reaction_annotation class instance.
  
  Parameters
  ----------
  inp_model: libsbml.Model
  
  Returns
  -------
  dict
  """
  exist_raw = {val.getId():extractRheaFromAnnotationString(val.getAnnotationString()) \
               for val in inp_model.getListOfReactions()}
  exist_filt = {val:exist_raw[val] for val in exist_raw.keys() \
                if exist_raw[val]}
  return exist_filt

def extractRheaFromAnnotationString(inp_str):
  """
  Extract Rhea from existing annotation string,
  by directly extracting Rhea,
  and converting from KEGGG and EC-Code. 
  
  Parameters
  ----------
  inp_str: str
  
  Returns
  -------
  list-str
  """
  exist_rheas = [cn.RHEA_HEADER+val for val in getQualifierFromString(inp_str, cn.RHEA)]
  map_rhea_bis = [cn.REF_RHEA2MASTER[val] for val in exist_rheas if val in cn.REF_RHEA2MASTER.keys()]

  exist_keggs = [cn.KEGG_HEADER+val for val in getQualifierFromString(inp_str, cn.KEGG_REACTION)]
  map_kegg2rhea = list(itertools.chain(*[cn.REF_KEGG2RHEA[val] \
                                         for val in exist_keggs if val in cn.REF_KEGG2RHEA.keys()]))

  exist_ecs = [cn.EC_HEADER+val for val in getQualifierFromString(inp_str, cn.EC)]
  map_ec2rhea = list(itertools.chain(*[cn.REF_EC2RHEA[val] \
                                      for val in exist_ecs if val in cn.REF_EC2RHEA.keys()]))

  return list(set(map_rhea_bis + map_kegg2rhea + map_ec2rhea))


def formatRhea(one_rhea):
  """
  Format rhea values; 
  if 'RHEA:' is not in the name,
  add it; if not, ignore it
  
  Parameters
  ----------
  str: one_rhea
  
  Returns
  -------
  :str
  """
  if one_rhea[:4].lower() == 'rhea':
    str_to_add = one_rhea[5:] 
  else:
    str_to_add = one_rhea
  return cn.RHEA_HEADER + str_to_add


def getOntologyFromString(string_annotation,
                          bqbiol_qualifiers=['is', 'isVersionOf']):
  """
  Parse string and return string annotation,
  marked as <bqbiol:is> or <bqbiol:isVersionOf>.
  If neither exists, return None.

  Parameters
  ----------
  string_annotation: str
  bqbiol_qualifiers: str-list
      Use 'is' and 'isVersionOf' by default
  

  Returns
  -------
  list-tuple (ontology type, ontology id)
       Return [] if none is provided
  """
  combined_str = ''
  for one_qualifier in bqbiol_qualifiers:
    one_match = '<bqbiol:' + one_qualifier + \
                '[^a-zA-Z].*?<\/bqbiol:' + \
                one_qualifier + '>'
    one_matched = re.findall(one_match,
                  string_annotation,
                  flags=re.DOTALL)
    if len(one_matched)>0:
      matched_filt = [s.replace("      ", "") for s in one_matched]
      one_str = '\n'.join(matched_filt) 
    else:
      one_str = ''
    combined_str = combined_str + one_str
  identifiers_list = re.findall('identifiers\.org/.*/', combined_str)
  result_identifiers = [(r.split('/')[1],r.split('/')[2].replace('\"', '')) \
                        for r in identifiers_list]
  return result_identifiers


def getQualifierFromString(input_str, qualifier):
  """
  Parses string and returns an identifier. 
  If not, return None.
  Qualifier is allowed to be
  either a string or a list of string. 

  Parameters
  ----------
  str/list-str: (list of) string_annotation

  Returns
  -------
  str (ontology Id)
      Returns an empty list if none is provided
  """
  ontologies = getOntologyFromString(input_str)
  # To make sure it works, make it lower
  if isinstance(qualifier, str):
    qualifier_list = [val for val in ontologies if val[0].lower()==qualifier.lower()]
  elif isinstance(qualifier, list):
    lower_qualifiers = [q.lower() for q in qualifier]
    qualifier_list = [val for val in ontologies \
                      if val[0].lower() in lower_qualifiers]
  if qualifier_list:
    return [val[1] for val in qualifier_list]
  else:
    return []


def getPrecision(ref, pred, mean=True):
  """
  (A model element-agnostic
  version of the method.)
  A complementary term of 'recall',
  precision is the fraction of correct
  elements detected from all detected elements. 

  Parameters
  ----------
  ref: dict
      {id: [str-annotation, e,g., formula/Rhea]}
  pred: dict
      {id: [str-annotation, e,g., formula/Rhea]} 
  mean: bool
      If True, get model-level average
      If False, get value of each ID

  Returns
  -------
  : float/dict{id: float}
      Depending on the 'mean' argument
  """
  ref_keys = set(ref.keys())
  pred_keys = set(pred.keys())
  precision = dict()
  # select species that can be evaluated
  species_to_test = ref_keys.intersection(pred_keys)
  # go through each species
  for one_k in species_to_test:
    num_intersection = len(set(ref[one_k]).intersection(pred[one_k]))
    num_predicted = len(set(pred[one_k]))
    precision[one_k] = num_intersection / num_predicted
  # return value is rounded up to the three decimal places
  if mean:
    return np.round(np.mean([precision[val] for val in precision.keys()]), cn.ROUND_DIGITS)
  else:
    return {val:np.round(precision[val],cn.ROUND_DIGITS) for val in precision.keys()}


def getRecall(ref, pred, mean=True):
  """
  (A model element-agnostic
  version of the method.)
  A precise version of 'accuracy',
  recall is the fraction of correct
  elements detected.
  Arguments are given as dictionaries. 

  Parameters
  ----------
  ref: dict
      {id: [str-annotation, e,g., formula/Rhea]}
      Annotations from reference. Considered 'correct'
  pred: dict
      {id: [str-annotation, e,g., formula/Rhea]}
      Annotations to be evaluated. 
  mean: bool
      If True, get the average across the keys.
      If False, get value of each key.

  Returns
  -------
  float/dict{id: float}
      Depending on the 'mean' argument
  """
  ref_keys = set(ref.keys())
  pred_keys = set(pred.keys())
  # select species that can be evaluated
  species_to_test = ref_keys.intersection(pred_keys)
  recall = dict()
  # go through each species
  for one_k in species_to_test:
    num_intersection = len(set(ref[one_k]).intersection(pred[one_k]))
    recall[one_k] = num_intersection / len(set(ref[one_k]))
  if mean:
    return np.round(np.mean([recall[val] for val in recall.keys()]), cn.ROUND_DIGITS)
  else:
    return {val:np.round(recall[val],cn.ROUND_DIGITS) for val in recall.keys()}


def transformCHEBIToFormula(inp_list, ref_to_formula_dict):
  """
  transform input list of CHEBI terms
  to list of annotations. 
  
  Parameters
  ----------
  inp_list: str-list
  
  Returns
  -------
  res: str-list
  """
  inp_formulas = [ref_to_formula_dict[val] for val in inp_list \
                  if val in ref_to_formula_dict.keys()]
  res = list(set([val for val in inp_formulas if val is not None]))
  return res


def updateDictKeyToList(inp_orig_dict, inp_new_dict):
  """
  Update inp_orig_dict using inp_up_dict.
  If key of inp_up_dict is already in inp_orig_dict,
  simply append the item list, 
  otherwise create a new list with a single item. 
  
  Parameters
  ----------
  inp_orig_dict: dict
      {key: [items]}
  inp_new_dict: dict
      {key: [items]} / {key: item}
      
  Returns
  -------
  res_dict: dict
      {key: [list of items]}
  """
  res_dict = inp_orig_dict.copy()
  # If nothing to update; return original dictionary
  if inp_new_dict is None:
    return res_dict
  for one_k in inp_new_dict.keys():
    # make item to a list, it is already not
    if isinstance(inp_new_dict[one_k], list):
      itm2add = inp_new_dict[one_k]
    else:
      itm2add = [inp_new_dict[one_k]]
    if one_k in res_dict.keys():
      res_dict[one_k] = list(set(res_dict[one_k] + itm2add))
    else:
      res_dict[one_k] = itm2add
  return res_dict

def getAssociatedTermsToRhea(inp_rhea):
  """
  Get a list of associated terms 
  of a rhea term. 
  The resulting list will contain 
  the original rhea term, 
  associated EC & KEGG numbers. 
  
  Parameters
  ----------
  inp_rhea: str
  
  Returns
  -------
  : list-str
  """
  if inp_rhea in cn.REF_RHEA2ECKEGG.keys():
    return cn.REF_RHEA2ECKEGG[inp_rhea] + [inp_rhea]
  else:
    return [inp_rhea]

  