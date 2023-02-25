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
                          meta_id):
    """
    Get a string of annotations,
    using a list of strings
    (of candidates)

    Parameters
    ----------
    candidates: list-str
        e.g., ['CHEBI:12345', 'CHEBI:98765']

    meta_id: str
        Meta ID of the element to be included in the annotation. 
    Returns
    -------
    str
    """
    # First, construct an empty container
    container_items = ['annotation', 
                       RDF_TAG,
                       'rdf:Description rdf:about="#'+meta_id+'"',
                       self.prefix,
                       'rdf:Bag']
    empty_container = self.createAnnotationContainer(container_items)
    # Next, create annotation lines
    items_from = []
    for one_cand in candidates:
      items_from.append(self.createAnnotationItem(KNOWLEDGE_RESOURCE[self.element],
                                                  one_cand))
    #
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







