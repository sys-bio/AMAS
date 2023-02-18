# annotation_maker.py
"""
Create string annotations for
AMAS recommendation.
"""


RDF_TAG_ITEM = ['rdf:RDF',
                'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"',
                'xmlns:dcterms="http://purl.org/dc/terms/"',
                'xmlns:vcard4="http://www.w3.org/2006/vcard/ns#"',
                'xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"',
                'xmlns:bqmodel="http://biomodels.net/model-qualifiers/"']
RDF_TAG = ' '.join(RDF_TAG_ITEM)

MATCH_SCORE_BY = {'species': 'by_name',
                  'reaction': 'by_component'}
KNOWLEDGE_RESOURCE = {'species': 'chebi',
                      'reaction': 'rhea'}


class AnnotationMaker(object):

  def __init__(self,
  	           element,
  	           prefix='bqbiol:is'):
    """
    Parameters
    ----------
    element: str
        Either 'species' or 'reaction'
        This will determine 
        the type of match score
        and the knowledge resource used. 
    """
    self.version = 'v1'
    self.prefix = prefix
    self.score_type = MATCH_SCORE_BY[element]
    self.knowledge_resource = KNOWLEDGE_RESOURCE[element]

  def createAnnotationBlock(self,
                            identifier,
                            match_score,
                            knowledge_resource=None,
                            score_type=None,
                            nested_prefix='bqbiol:hasProperty'):
    """
    An annotation block is a list of strings,
    with an annotation item at first
    and a list of string such as 
    <tag> <bag> <rdf:li 'amas_score'> </bag> </tag>.
    A block per each candidate is created,
    and they will be inserted into 
    the 'annotation container'.

    Parameters
    ----------
    identifier: str
        actual identifiers, e.g., 'CHEBI:59789'
    match_score: float/int/str
        match score
    knowledge_resource: str
        'chebi' or 'rhea'
        Optional unless want to use a different one
    score_type: str
        'by_name' for species or 
        'by_component' for reactions
        Optional unless user wants to specify it
    nested_prefix: str
        Default is 'bqbiol:hasProperty'

    Returns
    -------
    list-str
    """
    if knowledge_resource is None:
      knowledge_resource = self.knowledge_resource
    if score_type is None:
      score_type = self.score_type

    tags_nested = [nested_prefix, 'rdf:Bag']
    one_annotation = self.createAnnotationItem(knowledge_resource,
                                               identifier)   
    one_score = self.createScoreItem(match_score, self.score_type)
    # outer tags to be nested
    nest_container = []
    for idx_indent, one_str in enumerate(tags_nested):
      nest_container = self.insertEntry(one_str, nest_container, insert_loc=idx_indent)
    nested_block = self.insertEntry(inp_str=one_score,
                                    inp_list=nest_container,
                                    insert_loc=2,
                                    is_prefix=False)
    # Indentation for nested content unnecessary;
    # libsbml saved without it
    return [one_annotation] + nested_block

  def createAnnotationContainer(self, cont_items):
    """
    Create an empty annotation container
    that will hold the annotation blocks

    Parameters
    ----------
    cont_items: str-list

    Returns
    -------
    list-str
    """
    container =[]
    for one_str in cont_items:
      container = self.insertEntry(inp_str=one_str,
      	                           inp_list=container)
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


  def createScoreItem(self,
  	                  match_score,
  	                  score_type):
    """
    Create a one-line score annotation,
    e.g., <rdf:li rdf:resource="http://reproduciblebiomodels.org/amas/v1/by_name/1.00"/>
    Score type, e.g., name,
    is determined in __init__,
    when the type of model element was
    determined by either as species or reaction. 

    Parameters
    ----------
    match_score: float/str/int

    score_type: str

    Returns
    -------
    str
    """
    score_items = [self.version,
                   self.score_type,
                   format(match_score, '.2f')]
    res = '<rdf:li rdf:resource="http://reproduciblebiomodels.org/amas/' +\
          '/'.join(score_items)  +\
          '"/>'
    return res


  def createTag(self,
  	            tag_str,
  	            indent_val):
    """
    Create a tag based on the given string,
    with indent_val
   
    Parameters
    ---------
    str: inp_str
    indent_val: int
  
    Returns
    -------
    list-str
    """
    indent2pad = self.getIndent(indent_val)
    head_str = tag_str
    tail_str = tag_str.split(' ')[0]
    res_tag = [indent2pad+'<'+head_str+'>', indent2pad+'</'+tail_str+'>']
    return res_tag

  def getAnnotationString(self,
                          candidates,
                          meta_id):
    """
    Get a string of annotations,
    using a list of tuples (annotation, match_score).

    Parameters
    ----------
    candidates: list-tuple
        e.g., [(CHEBI:12345, 1.0), (CHEBI:98765, 0.8)]

    meta_id: str
        Meta ID of the element to be included in the annotation. 
    Returns
    -------
    str
    """
    items_from = []
    for one_cand in candidates:
      items_from = items_from + \
                   self.createAnnotationBlock(one_cand[0], one_cand[1])
    #
    container_items = ['annotation', 
                       RDF_TAG,
                       'rdf:Description rdf:about="#'+meta_id+'"',
                       self.prefix,
                       'rdf:Bag']
    empty_container = self.createAnnotationContainer(container_items)
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
                  insert_loc=None,
                  is_prefix=True):
    """
    Create an entry
  
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

    is_prefix: bool
        If False, insert as it is;
        If True, create <prefix> </prefix> to be inserted

    Returns
    -------
    : list-str
    """
    if insert_loc:
      idx_insert = insert_loc
    else:
      idx_insert = int(len(inp_list)/2)
    if is_prefix: 
      val2insert = self.createTag(tag_str=inp_str, indent_val=idx_insert)
    else:
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







