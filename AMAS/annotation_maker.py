# annotation_maker.py
"""
Create string annotations for
AMAS recommendation.
"""

import itertools
import re
from AMAS import tools
from AMAS import constants as cn

MATCH_SCORE_BY = {'species': 'by_name',
                  'reaction': 'by_component',
                  'genes':'by_name',
                  'qual_species':'by_name'}

KNOWLEDGE_RESOURCE = {'species': 'chebi',
                      'reaction': 'rhea',
                      'genes': 'ncbigene',
                      'qual_species':'ncbigene'}


class AnnotationMaker(object):

  def __init__(self,
  	           element,
  	           prefix='bqbiol:is'):
    """
    Parameters
    ----------
    element: str
        Either 'species' or 'reaction' or 'genes' or 'qual_species'
        This will determine 
        the type of match score
        and the knowledge resource used. 
    """
    self.qualifiers = ['bqbiol:is', 'bqbiol:isDescribedBy']
    self.prefix = prefix
    self.knowledge_resource = KNOWLEDGE_RESOURCE[element]
    # Below is only used when annotation line is created; 
    self.version = 'v1'
    self.element = element
    self.score_by = MATCH_SCORE_BY[element]


  def createAnnotationContainer(self, items):
    """
    Create an empty annotation container
    that will hold the annotation blocks

    Parameters
    ----------
    items: str-list

    Returns
    -------
    list-str
    """
    container =[]
    for one_item in items:
      one_t = self.createTag(one_item)
      container = self.insertList(insert_from=one_t,
                                  insert_to=container)
    return container

  def createAnnotationItem(self,
                           knowledge_resource,
                           identifier):
    """
    Create a one-line annotation,
    e.g., <rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:15414"/>

    Parameters
    ----------
    knowledge_resource: str

    identifier: str

    Returns
    -------
    str
    """
    annotation_items = ['identifiers.org',
                        knowledge_resource,
                        identifier]
    res = '<rdf:li rdf:resource="http://' + \
          '/'.join(annotation_items)  +\
          '"/>'
    return res

  def createTag(self,
                tag_str):
    """
    Create a tag based on the given string
   
    Parameters
    ---------
    str: inp_str
  
    Returns
    -------
    list-str
    """
    head_str = tag_str
    tail_str = tag_str.split(' ')[0]
    res_tag = ['<'+head_str+'>', '</'+tail_str+'>']
    return res_tag

  def getAnnotationString(self,
                          candidates,
                          meta_id,
                          cross_reference = None):
    """
    Get a string of annotations,
    using a list of strings.
    (of candidates)
    Can replace a whole annotation. 

    Parameters
    ----------
    candidates: list-str
        e.g., ['CHEBI:12345', 'CHEBI:98765']

    meta_id: str
        Meta ID of the element to be included in the annotation. 
    
    cross_reference: str
        Cross reference to be used for adding additional annotations.
        Single field or comma-separated. E.g., 'uniprot' or 'uniprot,HGNC'

    Returns
    -------
    str
    """
    # get the cross reference first if provided
    if cross_reference:
        all_cross_dict = tools.getCrossReference(candidates, fields = cross_reference)
        field_list = cross_reference.split(",") if "," in cross_reference else [cross_reference]
        
    # First, construct an empty container
    container_items = ['annotation', 
                       cn.RDF_TAG,
                       'rdf:Description rdf:about="#' + str(meta_id) + '"',
                       self.prefix,
                       'rdf:Bag']
    empty_container = self.createAnnotationContainer(container_items)
    # Next, create annotation lines
    items_from = []
    for one_cand in candidates:
        items_from.append(self.createAnnotationItem(KNOWLEDGE_RESOURCE[self.element],
                                                  one_cand))
        if cross_reference:
            for field in field_list:
                one_cand_cross = all_cross_dict[field][one_cand]
                if isinstance(one_cand_cross, list):
                    for i in one_cand_cross:
                        items_from.append(self.createAnnotationItem(field, i))    
                else:
                    items_from.append(self.createAnnotationItem(field, one_cand_cross))     

    result = self.insertList(insert_to=empty_container,
                                insert_from=items_from)
    return ('\n').join(result)


  def getIndent(self, num_indents=0):
    """
    Parameters
    ----------
    num_indents: int
      Time of indentation
    
    Returns
    -------
    :str
    """
    return '  ' * (num_indents)

  def insertEntry(self, 
                  inp_str,
                  inp_list=[],
                  insert_loc=None):
    """
    Insert a string into a list
  
    Parameters
    ----------
    inp_str: str
  
    inp_list: list
      New entry will be inserted in the middle.
      If not specified, will create a new list

    insert_loc: int
       If None, choose based on the middle of inp_list

    insert: bool
        If None, just return the create tag

    Returns
    -------
    : list-str
    """
    if insert_loc:
      idx_insert = insert_loc
    else:
      idx_insert = int(len(inp_list)/2)
    val2insert = [self.getIndent(idx_insert) + inp_str]
    return inp_list[:idx_insert] + val2insert + inp_list[idx_insert:]

  def insertList(self,
                 insert_to,
                 insert_from,
                 start_loc=None):
    """
    Insert a list to another list.

    Parameters
    ----------
    insert_to:list
        List where new list is inserted to

    inser_from: list
        A list where items will be inserted from

    start_loc: int
        If not given, insert_from will be 
        added in the middle of insert_to
    """
    if start_loc is None:
      start_loc = int(len(insert_to)/2)
    indents = self.getIndent(start_loc)
    insert_from_indented = [indents+val for val in insert_from]
    return insert_to[:start_loc] + \
           insert_from_indented + \
           insert_to[start_loc:]


  def addAnnotation(self, 
                    terms,
                    annotation,
                    meta_id=None):
    """
    Add terms to existing annotations
    (meta id is supposed to be included
    in the existing annotation)

    Parameters
    ----------
    terms: str-list
        List of terms to be added

    annotation: str
        Existing element annotation

    meta_id: str
        Optional argument; 
        if not provided and is needed,
        it'll extract appropriate one from annotation.

    Returns
    -------
    :str
    """
    qualifier_res = [q for q in self.qualifiers if q != self.prefix][0]
    # Attempt to divide the existing annotation
    annotation_dict = tools.divideExistingAnnotation(annotation, qualifier = self.prefix)
    annotation_dict_res = tools.divideExistingAnnotation(annotation, qualifier = qualifier_res)

    # If there is no existing annotation, create a new one
    if not annotation_dict:
      if not meta_id: 
          meta_id = self.extractMetaID(annotation)
          
      # if there exists an annotation for the other qualifier, add the terms to that annotation
      if annotation_dict_res:
        empty_container = self.createEmptyContainerWithQualifiers(meta_id)
        items = self.cleanItems(annotation_dict_res['items'])
        container = tools.insertItemsBackToContainer(empty_container, items, qualifier = qualifier_res)
        items = [self.createAnnotationItem(KNOWLEDGE_RESOURCE[self.element], one_cand) for one_cand in terms]
    
    # if there is no annotation for the other qualifier, create a new one
      else:
        return self.getAnnotationString(terms, meta_id)
    
    else:
      # Process existing annotations
      container = annotation_dict['container']
      existing_items = annotation_dict['items']
      existing_identifiers = []
      for val in existing_items:
          url = re.findall('"(.*?)"', val)[0]
          existing_identifiers.append(url.split('/')[-1])
      
      # Add only new terms that are not already present
      additional_identifiers = [val for val in terms if val not in existing_identifiers]
      new_items = [self.createAnnotationItem(KNOWLEDGE_RESOURCE[self.element], one_cand) for one_cand in additional_identifiers]
      items = existing_items + new_items 

    # remove duplicates and sort alphabetically
    items = sorted(list(set(items)))
    res = tools.insertItemsBackToContainer(container, items, qualifier = self.prefix)
    return res


  def deleteAnnotation(self,
                       terms,
                       annotation):
    """
    Remove entire annotation by 
    returning a null string.

    Parameters
    ----------
    terms: str-list
        List of terms to be removed
      
    annotation: str
        Existing element annotation

    Returns
    -------
    :str
    """
    annotation_dict = tools.divideExistingAnnotation(annotation, qualifier = self.prefix)
    # if cannot parse annotation, return the original annotation
    if annotation_dict is None:
      return annotation
    container = annotation_dict['container']
    exist_items = annotation_dict['items']
    # finding remaining items
    rem_items = []
    for val in exist_items:
      if all([k not in val for k in terms]):
        rem_items.append(val)
    if rem_items:
      res = tools.insertItemsBackToContainer(container, rem_items, qualifier = self.prefix)
      return res
    # if all items were deleted, return an empty string
    else:
      return ''

  def extractMetaID(self,
                    inp_str): 
    """
    Extract meta id from
    the given annotation string, by searching for
    two strings: '#metaid_' and '">'.
    If none found, return an emtpy string

    Parameters
    ----------
    inp_str: str
        Annotation string

    Returns
    -------
    :str
        Extracted meta id
    """
    metaid_re = re.search('rdf:about="#(.*)">', inp_str)
    if metaid_re is None:
      return ''
    else:
      return metaid_re.group(1)

  def convertIsDescribedByToIs(self, inp_str, meta_id):
    """
    Convert items in <bqbiol:isDescribedBy> to <bqbiol:is>
    Will not move pubmed annotations

    Parameters
    ----------
    inp_str: str
        Existing annotation string

    meta_id: str
        Meta ID of the element to be included in the annotation. 

    Returns
    -------
    :str
        The updated annotation string with <bqbiol:is> items.
    """
    # find all <bqbiol:isDescribedBy> items
    des_dict = tools.divideExistingAnnotation(inp_str, qualifier = 'bqbiol:isDescribedBy')
    # find all <bqbiol:is> items
    is_dict = tools.divideExistingAnnotation(inp_str, qualifier = 'bqbiol:is')
    # construct a new annotation container
    new_container = self.createEmptyContainerWithQualifiers(meta_id)
    # convert each item to <bqbiol:is>
    if is_dict is not None:
      if des_dict is not None:
        items = des_dict['items']+is_dict['items'] 
      else:
        items = is_dict['items']
    else:
      if des_dict is not None:
        items = des_dict['items']
      else:
        items = []
    # formatting items
    items = self.cleanItems(items)
    pubmed_items = [item for item in items if 'pubmed' in item]
    non_pubmed_items = [item for item in items if 'pubmed' not in item]
    container = tools.insertItemsBackToContainer(new_container, pubmed_items, qualifier = 'bqbiol:isDescribedBy')
    res = tools.insertItemsBackToContainer(container, non_pubmed_items, qualifier = 'bqbiol:is')
    return res

  def createEmptyContainerWithQualifiers(self, meta_id):
      """
      Create a new empty annotation container that includes both bqbiol:is
      and bqbiol:isDescribedBy qualifiers.

      Parameters
      ----------
      meta_id : str
          The meta ID of the element to include in the annotation.

      Returns
      -------
      str
          The annotation string with empty blocks for bqbiol:is and bqbiol:isDescribedBy.
      """
      qualifiers = ['bqbiol:is', 'bqbiol:isDescribedBy']

      # Construct the container with RDF_TAG and meta_id
      container = f'<annotation>\n  <{cn.RDF_TAG}>\n    <rdf:Description rdf:about="#{meta_id}">\n'

      for qualifier in qualifiers:
          container += f'      <{qualifier}>\n        <rdf:Bag>\n        </rdf:Bag>\n      </{qualifier}>\n'

      container += '    </rdf:Description>\n  </rdf:RDF>\n</annotation>'
      return container

  def cleanItems(self, items):
    """
    Clean items by formatting them.
    Before; '<rdf:li rdf:resource="urn:miriam:ncbigene:1489671"/>'
    After; '<rdf:li rdf:resource="http://identifiers.org/ncbigene/1489671"/>'
    """
    ontologies = tools.extract_ontology_from_items(items)
    cleaned_items = []
    for ontology in ontologies:
        res = self.createAnnotationItem(ontology[0], ontology[1])
        cleaned_items.append(res)
    cleaned_items = sorted(list(set(cleaned_items))) # remove duplicates and sort alphabetically
    return cleaned_items

  def addCrossReference(self, inp_str, fields):
    """
    Add cross reference of NCBI Gene ID to existing annotations.

    Parameters
    ----------
    inp_str: str
        Existing annotation string

    fields: str
        Fields to add to annotations. Single field or comma-separated. E.g., 'uniprot' or 'uniprot,HGNC'.

    Returns
    -------
    :str
        The updated annotation string with cross reference items given in fields.
    """
    # find all existing items, only <bqbiol:is> items
    annotation_dict = tools.divideExistingAnnotation(inp_str, qualifier = 'bqbiol:is')
    if annotation_dict is None:
      return inp_str
    items = annotation_dict['items']
    # find all existing ncbi gene ids
    ontology_list = tools.extract_ontology_from_items(items)
    ontology_type_list = [val for val in ontology_list if val[0].lower()=='ncbigene']
    ncbi_ids = [val[1] for val in ontology_type_list]
    # get the cross reference items of the ncbi ids
    cross_dict = tools.getCrossReference(ncbi_ids, fields)
    # add the cross reference items to the annotation
    field_list = fields.split(",") if "," in fields else [fields]
    field_list = [field.split(".")[0] if "." in field else field for field in field_list]
    for field in field_list:
      for ncbi_id in ncbi_ids:
        if cross_dict[field][ncbi_id] is not None:
          for one_item in cross_dict[field][ncbi_id]:
              items.append(self.createAnnotationItem(field, one_item))
    items = sorted(list(set(items))) # remove duplicates and sort alphabetically
    container = annotation_dict['container']
    res = tools.insertItemsBackToContainer(container, items, qualifier = 'bqbiol:is')
    return res

  # def extractSBMLQualItems(self, inp_str):
  #     """
  #     Extract unique <rdf:li> items from SBML-qual string annotation.

  #     Parameters
  #     ----------
  #     inp_str : str
  #         Input annotation string.

  #     Returns
  #     -------
  #     list
  #         List of unique <rdf:li> items.
  #     """
  #     items = set()  # Use a set to ensure uniqueness

  #     # Split the input string into lines
  #     exist_anot_list = inp_str.split('\n')

  #     for line in exist_anot_list:
  #         stripped_line = line.strip()
  #         # Extract <rdf:li> lines
  #         if stripped_line.startswith('<rdf:li'):
  #             items.add(stripped_line)

  #     # Return items as a list
  #     return list(items)

  # def makeAnnotationString(self,
  #                         items,
  #                         meta_id):
  #   """
  #   Get a string of annotations,
  #   using a list of items.
  #   Can replace a whole annotation. 

  #   Parameters
  #   ----------
  #   items: list-str
  #       e.g., ['<rdf:li rdf:resource="http://identifiers.org/ncbigene/12345"/>']

  #   meta_id: str
  #       Meta ID of the element to be included in the annotation. 
  #   Returns
  #   -------
  #   str
  #   """
  #   # First, construct an empty container
  #   container_items = ['annotation', 
  #                      RDF_TAG,
  #                      'rdf:Description rdf:about="#' + str(meta_id) + '"',
  #                      self.prefix,
  #                      'rdf:Bag']
  #   empty_container = self.createAnnotationContainer(container_items)
  #   # Next, create annotation lines
  #   items_from = []
  #   for one_item in items:
  #     items_from.append(one_item)
                                                  
  #   result = self.insertList(insert_to=empty_container,
  #                            insert_from=items_from)
  #   return ('\n').join(result)

