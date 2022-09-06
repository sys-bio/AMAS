# tools.py

import re


# def getOntologyFromString(string_annotation):
#   """
#   Parse string and return string annotation,
#   marked as <bqbiol:is> or <bqbiol:isVersionOf>.
#   If neither exists, return None.

#   Parameters
#   ----------
#   str: string_annotation

#   Returns
#   -------
#   list-tuple (ontology type, ontology id)
#        Return [] if none is provided
#   """
#   # first, extracts strings tagged as bqbiol:is or bqbiol:isVersionOf.
#   is_str = ''
#   isVersionOf_str = ''
#   is_str_match = re.findall('<bqbiol:is[^a-zA-Z].*?<\/bqbiol:is>',
#                             string_annotation,
#                             flags=re.DOTALL)
#   if len(is_str_match)>0:
#     is_str_match_filt = [s.replace("      ", "") for s in is_str_match]
#     is_str = '\n'.join(is_str_match_filt)  
#   #
#   is_VersionOf_str_match = re.findall('<bqbiol:isVersionOf[^a-zA-Z].*?<\/bqbiol:isVersionOf>',
#                                       string_annotation,
#                                       flags=re.DOTALL)  
#   #
#   if len(is_VersionOf_str_match) > 0:
#     is_VersionOf_str_match_filt = [s.replace("      ", "") for s in is_VersionOf_str_match]
#     isVersionOf_str = '\n'.join(is_VersionOf_str_match_filt) 
#   #
#   combined_str = is_str + isVersionOf_str
#   if combined_str == '':
#     return []
#   identifiers_list = re.findall('identifiers\.org/.*/', combined_str)
#   return [(r.split('/')[1],r.split('/')[2].replace('\"', '')) \
#           for r in identifiers_list]


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
      Return None if none is provided
  """
  ontologies = getOntologyFromString(input_str)
  # To make sure it works, make it lower
  qualifier_list = [val for val in ontologies if val[0]==qualifier.lower()]
  if qualifier_list:
    return [val[1] for val in qualifier_list]
  else:
    return None


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





  