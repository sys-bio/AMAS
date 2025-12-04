# tools.py

from AMAS import constants as cn

import itertools
import numpy as np
import re
import os
import compress_pickle
import requests
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
  exist_raw = {val.getId():getOneOntologyFromString(val.getAnnotationString(), [cn.CHEBI, cn.OBO_CHEBI]) \
               for val in inp_model.getListOfSpecies()}
#   print(f"exist_raw: {exist_raw}")
  exist_filt = {val:exist_raw[val] for val in exist_raw.keys() \
                if exist_raw[val]}
#   print(f"exist_filt: {exist_filt}")
  # remove duplicates
  exist_filt = {val:list(set(exist_filt[val])) for val in exist_filt.keys()}
  return exist_filt

def extractExistingGeneAnnotation(inp_model, qualifier=cn.NCBI_GENE):
  """
  Get existing annotation of genes
  that contains NCBI gene terms

  Parameters
  ---------
  qualifier: str
      'ncbi.gene'
  """
  exist_raw = {val.getIdAttribute():getOneOntologyFromString(val.getAnnotationString(), [cn.NCBI_GENE]) \
               for val in inp_model.getListOfGeneProducts()}
  exist_filt = {val:exist_raw[val] for val in exist_raw.keys() \
                if exist_raw[val]}
  # remove duplicates
  exist_filt = {val:list(set(exist_filt[val])) for val in exist_filt.keys()}
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
  # remove duplicates
  exist_filt = {val:list(set(exist_filt[val])) for val in exist_filt.keys()}
  return exist_filt

def extractExistingQualitativeSpeciesAnnotation(inp_model, qualifier=cn.NCBI_GENE, description=False):
  """
  Get existing annotation of qualitative species
  that contains gene terms (by default, NCBI gene)

  Parameters
  ---------
  qualifier: str
      'ncbigene'
  description: bool
      If True, also extract the description using <bqbiol:isDescribedBy>
  """
  exist_raw = {val.getId():getOneOntologyFromString(val.getAnnotationString(), [qualifier], description) \
               for val in inp_model.getListOfQualitativeSpecies()}
  exist_filt = {val:exist_raw[val] for val in exist_raw.keys() \
                if exist_raw[val]}
  # remove duplicates
  exist_filt = {val:list(set(exist_filt[val])) for val in exist_filt.keys()}
  return exist_filt

def extractExistingTransitionAnnotation(inp_model, qualifier=cn.PUBMED, description=False):
  # TODO: need to be tested
  """
  Get existing annotation of transitions
  that contains pubmed terms

  Parameters
  ---------
  qualifier: str
      'pubmed'
  """
  exist_raw = {val.getId():getOneOntologyFromString(val.getAnnotationString(), [qualifier], description) \
               for val in inp_model.getListOfTransitions()}
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
  exist_rheas = [formatRhea(val) for val in getOneOntologyFromString(inp_str, cn.RHEA)]
  map_rhea_bis = [cn.REF_RHEA2MASTER[val] for val in exist_rheas if val in cn.REF_RHEA2MASTER.keys()]

  exist_keggs = [cn.KEGG_HEADER+val for val in getOneOntologyFromString(inp_str, cn.KEGG_REACTION)]
  map_kegg2rhea = list(itertools.chain(*[cn.REF_KEGG2RHEA[val] \
                                         for val in exist_keggs if val in cn.REF_KEGG2RHEA.keys()]))

  exist_ecs = [cn.EC_HEADER+val for val in getOneOntologyFromString(inp_str, cn.EC)]
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

def remove_duplicate_namespaces(xml_string):
    """
    Remove duplicate namespace declarations from the RDF_TAG in the XML string.

    Parameters
    ----------
    xml_string: str
        The XML string to clean.

    Returns
    -------
    str:
        The cleaned XML string with duplicate attributes removed.
    """
    # Regex to locate the <rdf:RDF> opening tag and extract its attributes
    rdf_opening_pattern = r"(<rdf:RDF\b)([^>]*)(>)"
    match = re.search(rdf_opening_pattern, xml_string)

    if match:
        opening_tag_start = match.group(1)  # '<rdf:RDF'
        attributes = match.group(2)  # All attributes
        closing_tag = match.group(3)  # '>'

        # Split attributes and retain only unique ones
        attr_list = attributes.split()
        unique_attrs = {}
        for attr in attr_list:
            key, value = attr.split("=", 1)
            if key not in unique_attrs:
                unique_attrs[key] = value

        # Reconstruct the cleaned <rdf:RDF> opening tag
        cleaned_attributes = " ".join([f'{key}={value}' for key, value in unique_attrs.items()])
        cleaned_rdf_opening = f"{opening_tag_start} {cleaned_attributes}{closing_tag}"

        # Replace the original <rdf:RDF> opening tag with the cleaned one
        xml_string = xml_string.replace(match.group(0), cleaned_rdf_opening)

    return xml_string

def divideExistingAnnotation(inp_str, qualifier):
    """
    Parse the annotation string to extract items in multiple <rdf:Bag> elements 
    under the specified qualifier (e.g., <bqbiol:is>), and keep the rest in the container.
    Ensures namespaces and structure are preserved when creating an empty qualifier block.

    Parameters
    ----------
    inp_str: str
        The full annotation string.
    qualifier: str
        The qualifier to target (e.g., bqbiol:is, bqbiol:isDescribedBy).

    Returns
    -------
    dict:
        - 'container': The annotation string with the qualifier blocks replaced by a single empty block.
        - 'items': A list of <rdf:li> elements found in all <rdf:Bag> containers of the specified qualifier.
        - If no qualifier blocks are found, returns None.
    """

    # Regex to match all blocks for the specified qualifier, including attributes
    qualifier_pattern = rf"(<{qualifier}\b[^>]*?>\s*<rdf:Bag>.*?</rdf:Bag>\s*</{qualifier}>)"
    qualifier_matches = re.findall(qualifier_pattern, inp_str, re.DOTALL)

    if not qualifier_matches:
        return None  # Return None if no blocks for the qualifier are found

    # Collect all <rdf:li> elements from each matched block
    rdf_li_pattern = r"<rdf:li[^>]*\/>"
    items = []
    for block in qualifier_matches:
        items.extend(re.findall(rdf_li_pattern, block))

    # Extract the opening tag with attributes from the first match
    match_prefix = re.match(rf"<{qualifier}.*?>", qualifier_matches[0])
    if match_prefix:
        qualifier_opening = match_prefix.group()  # Capture opening tag with attributes
    else:
        qualifier_opening = f"<{qualifier}>"

    # Construct a single empty qualifier block using the preserved attributes
    empty_qualifier_block = (
        f"      {qualifier_opening}\n"
        f"        <rdf:Bag>\n"
        f"        </rdf:Bag>\n"
        f"      </{qualifier}>"
    )
    # Remove all original qualifier blocks from the container
    stripped_annotation = re.sub(qualifier_pattern, "", inp_str, flags=re.DOTALL).strip()

    # Reinsert the empty qualifier block with namespaces added
    if cn.RDF_TAG not in stripped_annotation:
        container = stripped_annotation.replace(
            "<rdf:RDF",
            f"<{cn.RDF_TAG}",
        )
    else:
        container = stripped_annotation

    container = container.replace(
        "</rdf:Description>",
        f"{empty_qualifier_block}\n    </rdf:Description>"
    )
    # Remove repeated blank lines
    container = re.sub(r"\n\s*\n", "\n", container)
    container = remove_duplicate_namespaces(container)
    return {"container": container, "items": items}


def insertItemsBackToContainer(container, items, qualifier):
    """
    Insert <rdf:li> items back into the <rdf:Bag> of the specified qualifier in the container.
    If no block for the specified qualifier exists, create a new one with attributes preserved.

    Parameters
    ----------
    container: str
        The annotation string with the specified qualifier and its <rdf:Bag>.
    items: list
        List of <rdf:li> items to insert back into the <rdf:Bag> of the qualifier.
    qualifier: str
        The qualifier to target (e.g., bqbiol:is, bqbiol:isDescribedBy).

    Returns
    -------
    str:
        The updated annotation string with the <rdf:li> items inserted into the correct <rdf:Bag>.
    """
    # Regex to locate the <rdf:Bag> inside the specified qualifier
    bag_pattern = rf"(<{qualifier}\b[^>]*?>\s*<rdf:Bag>).*?</rdf:Bag>\s*</{qualifier}>"
    match = re.search(bag_pattern, container, re.DOTALL)

    # Construct the new <rdf:Bag> content with the items
    items_str = "\n".join([f"          {item}" for item in items])

    if match:
        # If qualifier block exists, replace its <rdf:Bag> content
        qualifier_opening = match.group(1)  # Capture opening tag with attributes
        updated_bag = f"{qualifier_opening}\n{items_str}\n        </rdf:Bag>"
        # Replace the old <rdf:Bag> block with the updated one
        updated_container = re.sub(bag_pattern, updated_bag + f"\n      </{qualifier}>", container, flags=re.DOTALL)
    else:
        # If qualifier block does not exist, create a new one
        qualifier_opening = f"<{qualifier} xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\">"
        updated_bag = f"      {qualifier_opening}\n        <rdf:Bag>\n{items_str}\n        </rdf:Bag>\n      </{qualifier}>"
        if cn.RDF_TAG not in container:
          updated_container = container.replace(
              "<rdf:RDF",
              f"<{cn.RDF_TAG}")
        updated_container = updated_container.replace(
            "</rdf:Description>",
            f"{updated_bag}\n    </rdf:Description>"
        )
    updated_container = remove_duplicate_namespaces(updated_container)

    return updated_container


def extract_ontology_from_items(items_list):
    """
    Extract ontology from items and return a flat list of tuples.
    Each tuple contains (ontology type, ontology id).

    Parameters
    ----------
    items_list : list
        A list of string items containing ontology annotations.

    Returns
    -------
    list of tuples
        A flat list of (ontology type, ontology id).
    """
    result_identifiers = []
    for item in items_list:
        # Extract identifiers from "urn:miriam" URIs
        miriam_identifiers = re.findall(r'urn:miriam:([^"]+)"', item)
        for r in miriam_identifiers:
            ontology_type, ontology_id = r.split(":", 1)
            result_identifiers.append((ontology_type, ontology_id))

        # Extract identifiers from "identifiers.org" URIs
        identifiers_list = re.findall(r'identifiers\.org/([^/]+)/([^/"]+)', item)
        for ontology_type, ontology_id in identifiers_list:
            # Clean ontology_id to remove any unwanted characters
            ontology_id = ontology_id.replace('"', '')
            result_identifiers.append((ontology_type, ontology_id))

    return result_identifiers


def getOntologyFromString(string_annotation, description=False):
    """
    Parse string and return string annotation,
    marked as <bqbiol:is> or <bqbiol:isVersionOf>;
    (and extract the description using <bqbiol:isDescribedBy>)
    If neither exists, return None.

    Parameters
    ----------
    string_annotation: annotation string that starts with <annotation>

    description: bool
        If True, also extract the description using <bqbiol:isDescribedBy>.

    Returns
    -------
    list-tuple (ontology type, ontology id)
         Return [] if none is provided.
    """

    # Define the bqbiol qualifiers to search
    bqbiol_qualifiers = ['is', 'isVersionOf']
    if description:
        bqbiol_qualifiers.append('isDescribedBy')

    result_identifiers = []
    for one_qualifier in bqbiol_qualifiers:
        res = divideExistingAnnotation(string_annotation, 'bqbiol:' + one_qualifier)
        if res is None:
            continue
        items = res['items']
        res = extract_ontology_from_items(items)
        # remove empty lists and duplicates
        res = [r for r in res if r is not None and r != [] and r not in result_identifiers]
        result_identifiers.extend(res)

    return result_identifiers

# def getOntologyFromString(string_annotation,
#                           description = False):
#   """
#   Parse string and return string annotation,
#   marked as <bqbiol:is> or <bqbiol:isVersionOf>;
#   (and extract the description using <bqbiol:isDescribedBy>)
#   If neither exists, return None.

#   Parameters
#   ----------
#   string_annotation: str
#   description: bool
#       If True, also extract the description using <bqbiol:isDescribedBy>


#   Returns
#   -------
#   list-tuple (ontology type, ontology id)
#        Return [] if none is provided
  
#   """
#   bqbiol_qualifiers=['is', 'isVersionOf']
#   if description:
#     bqbiol_qualifiers.append('isDescribedBy')
#   combined_str = ''
#   for one_qualifier in bqbiol_qualifiers:
#     one_match = '<bqbiol:' + one_qualifier + \
#                 '[^a-zA-Z].*?<\/bqbiol:' + \
#                 one_qualifier + '>'
#     one_matched = re.findall(one_match,
#                   string_annotation,
#                   flags=re.DOTALL)
#     if len(one_matched)>0:
#       matched_filt = [s.replace("      ", "") for s in one_matched]
#       one_str = '\n'.join(matched_filt) 
#     else:
#       one_str = ''
#     combined_str = combined_str + one_str

#   identifiers_list = re.findall('identifiers\.org/.*/', combined_str)
#   result_identifiers = [(r.split('/')[1],r.split('/')[2].replace('\"', '')) \
#                         for r in identifiers_list]
#   if description:
#     identifiers_list.extend(re.findall(r'rdf:resource="urn:miriam:([^"]+)"', combined_str))
#     result_identifiers.extend([(r.rsplit(':', 1)[0], r.rsplit(':', 1)[1]) for r in identifiers_list])

#   return result_identifiers


def getOneOntologyFromString(input_str, ontology_type,
                          description = False):
  """
  Parses string and returns an identifier. 
  If not, return None.
  Ontology type is allowed to be either a string or a list of string. 
  Usage; getOneOntologyFromString(annotationstring, [cn.NCBI_GENE], description)

  Parameters
  ----------
  input_str: str: string_annotation
  ontology_type: str/list-str: (list of) ontology type for resources
  description: bool
      If True, also extract the description using <bqbiol:isDescribedBy>

  Returns
  -------
  str/list-str (ontology Id)
      Returns an empty list if none is provided
  """
  ontologies = getOntologyFromString(input_str, description)
  # To make sure it works, make it lower
  if isinstance(ontology_type, str):
    ontology_type_list = [val for val in ontologies if val[0].lower()==ontology_type.lower()]
  elif isinstance(ontology_type, list):
    lower_ontology_types = [q.lower() for q in ontology_type]
    ontology_type_list = [val for val in ontologies \
                      if val[0].lower() in lower_ontology_types]
  if ontology_type_list:
    return [val[1] for val in ontology_type_list]
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
    if num_predicted == 0:
        # Avoid division by zero
        precision[one_k] = 0.0
    else:
        precision[one_k] = num_intersection / num_predicted
  # return value is rounded up to the three decimal places
  if mean:
      if precision:
          return np.round(np.mean([precision[val] for val in precision.keys()]), 3)
      else:
          # No species to test, return 0.0
          return 0.0
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

def getTopKAccuracy(ref, pred, k=1, mean=True):
  """
  Calculate top-k accuracy for annotations.
  Top-k accuracy measures if any of the k most confident
  predictions match with any reference annotation.

  Parameters
  ----------
  ref: dict
      {id: [str-annotation, e.g., formula/Rhea]}
      Annotations from reference. Considered 'correct'
  pred: dict
      {id: [str-annotation, e.g., formula/Rhea]}
      Annotations to be evaluated. Assumes the order
      represents confidence (first item is highest confidence)
  k: int
      Number of top predictions to consider
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
  # select elements that can be evaluated
  elements_to_test = ref_keys.intersection(pred_keys)
  top_k_accuracy = dict()
  
  # go through each element
  for one_k in elements_to_test:
    # Take only top k predictions
    top_k_preds = pred[one_k][:min(k, len(pred[one_k]))]
    # Check if any of the top k predictions match any reference
    matches = set(top_k_preds).intersection(set(ref[one_k]))
    # If there's at least one match, accuracy is 1, otherwise 0
    top_k_accuracy[one_k] = 1.0 if len(matches) > 0 else 0.0
  
  if mean:
    if top_k_accuracy:
      return np.round(np.mean([top_k_accuracy[val] for val in top_k_accuracy.keys()]), cn.ROUND_DIGITS)
    else:
      # No elements to test, return 0.0
      return 0.0
  else:
    return {val:np.round(top_k_accuracy[val], cn.ROUND_DIGITS) for val in top_k_accuracy.keys()}

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

def get_gene_char_count_df(tax_id):
  """
  Get the gene character count dataframe for a given tax_id
  
  Parameters
  ----------
  tax_id: str
      Taxonomy ID to load data for
      
  Returns
  -------
  DataFrame
      Gene character count dataframe for the taxonomy
      
  Raises
  ------
  FileNotFoundError
      If the data file for the given tax_id does not exist
  """
  filepath = os.path.join(cn.REF_GENE_CHARCOUNT_DIR, f"charcount_gene_df_scaled_{tax_id}.lzma")
  if not os.path.exists(filepath):
      raise FileNotFoundError(f"Reference data not found for tax_id {tax_id}")
  return compress_pickle.load(filepath)


def getCrossReference(gene_ids, fields = "uniprot"):
    """
    Query the mygene.info API to retrieve EntrezGene IDs and UniProt data 
    for given gene ids (ncbigene, or entrezgene), then process the results into a dictionary.

    Parameters:
    - gene_ids (list): List of gene ids (ncbigene, or entrezgene) to query.
    - fields (str): Fields to query (e.g., "uniprot"). Single field or comma-separated.
        - For available fields, see https://docs.mygene.info/en/latest/doc/data.html

    Returns:
    - dict: A dictionary structured as {field: {gene_id: value}}
            If 'uniprot' is included, 'Swiss-Prot' is prioritized compared to 'TrEMBL'.
    """
    if gene_ids is None:
        return {}
    url = "https://mygene.info/v3/query"
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }
    data = {
        "q": gene_ids,
        "fields": fields
    }
    response = requests.post(url, json=data, headers=headers)
    if response.status_code != 200:
        return {}
    results = response.json()

    field_list = fields.split(",") if "," in fields else [fields]
    # If the field is a nested field, e.g., ensembl.gene, refseq.rna, etc.,
    # only the first part of the field is used for the dictionary key.
    field_list = [field.split(".")[0] if "." in field else field for field in field_list]
    gene_dict = {field: {} for field in field_list}

    for entry in results:
        query_key = str(entry.get("query"))
        for field in field_list:
            if field == "uniprot":
                uniprot_info = entry.get("uniprot", {})
                swiss_prot = uniprot_info.get("Swiss-Prot")
                trembl = uniprot_info.get("TrEMBL", [])
                if swiss_prot:
                    gene_dict["uniprot"][query_key] = [swiss_prot] if isinstance(swiss_prot, str) else swiss_prot
                elif trembl:
                    gene_dict["uniprot"][query_key] = [trembl] if isinstance(trembl, str) else trembl
                else:
                    gene_dict["uniprot"][query_key] = None
            elif field.startswith("ensembl"):
                ensembl_info = entry.get("ensembl", {})
                if field == "ensembl":
                    ensembl_gene = ensembl_info.get("gene", []) if isinstance(ensembl_info.get("gene", []), list) else [ensembl_info.get("gene", [])]
                    ensembl_protein = ensembl_info.get("protein", []) if isinstance(ensembl_info.get("protein", []), list) else [ensembl_info.get("protein", [])]
                    gene_dict["ensembl"][query_key] = ensembl_gene + ensembl_protein
                else:
                    ensembl_field = field.split(".")[1]
                    gene_dict["ensembl"][query_key] = ensembl_info.get(ensembl_field, []) if isinstance(ensembl_info.get(ensembl_field, []), list) else [ensembl_info.get(ensembl_field, [])]
            elif field.startswith("refseq"):
                refseq_info = entry.get("refseq", {})
                if field == "refseq":
                    refseq_rna = refseq_info.get("rna", []) if isinstance(refseq_info.get("rna", []), list) else [refseq_info.get("rna", [])]
                    refseq_protein = refseq_info.get("protein", []) if isinstance(refseq_info.get("protein", []), list) else [refseq_info.get("protein", [])]
                    gene_dict["refseq"][query_key] = refseq_rna + refseq_protein
                else:
                    refseq_field = field.split(".")[1]
                    gene_dict["refseq"][query_key] = refseq_info.get(refseq_field, []) if isinstance(refseq_info.get(refseq_field, []), list) else [refseq_info.get(refseq_field, [])]
            else:
                cross_ref = entry.get(field, None)
                gene_dict[field][query_key] = [cross_ref] if isinstance(cross_ref, str) else cross_ref
    return gene_dict