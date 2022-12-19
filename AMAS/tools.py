# tools.py

from AMAS import constants as cn

import itertools
import numpy as np
import re


def extractExistingSpeciesAnnotation(inp_model, qualifier=cn.CHEBI):
  """
  Get existing annotation of species
  that contains ChEBI terms
  """
  exist_raw = {val.getId():getQualifierFromString(val.getAnnotationString(), qualifier) \
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
  map_rhea_bis = [cn.REF_RHEA2BI[val] for val in exist_rheas if val in cn.REF_RHEA2BI.keys()]

  exist_keggs = [val for val in getQualifierFromString(inp_str, cn.KEGG_REACTION)]
  map_kegg2rhea = [cn.REF_KEGG2RHEA_BI[val] for val in exist_keggs if val in cn.REF_KEGG2RHEA_BI.keys()]

  exist_ecs = [cn.EC_HEADER+val for val in getQualifierFromString(inp_str, cn.EC)]
  map_ec2rhea = list(itertools.chain(*[cn.REF_EC2RHEA_BI[val] \
                                      for val in exist_ecs if val in cn.REF_EC2RHEA_BI.keys()]))

  return list(set(map_rhea_bis + map_kegg2rhea + map_ec2rhea))


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
  If not, return None

  Parameters
  ----------
  str: string_annotation

  Returns
  -------
  str (ontology Id)
      Returns an empty list if none is provided
  """
  ontologies = getOntologyFromString(input_str)
  # To make sure it works, make it lower
  qualifier_list = [val for val in ontologies if val[0]==qualifier.lower()]
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
    precision[one_k] = num_intersection / len(set(pred[one_k]))
  if mean:
    return np.mean([precision[val] for val in precision.keys()])
  else:
    return precision


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
    return np.mean([recall[val] for val in recall.keys()])
  else:
    return recall


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





  