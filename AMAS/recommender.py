# recommender.py
# Recomender for running annotation predictions

import libsbml
import os
import pickle

from AMAS import constants as cn
from AMAS import tools
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra

with open(os.path.join(cn.REF_DIR, 'kegg2rhea_bi.pickle'), 'rb') as handle:
  ref_kegg2rhea_bi = pickle.load(handle)


class Recommender(object):

  def __init__(self, libsbml_fpath=None, exist_qualifier=cn.RHEA):
  	# creates both species & reaction annotation class instances
  	# 
    pass
    # First of all, collect model information from libsbml model
    # and send the informaton to create species/reaction annotations
    if libsbml_fpath:
      reader = libsbml.SBMLReader()
      document = reader.readSBML(libsbml_fpath)
      self.model = document.getModel()
      # Create species_annotation instance
      exist_spec_annotation_raw = {val.getId():tools.getQualifierFromString(val.getAnnotationString(), cn.CHEBI) \
                                   for val in self.model.getListOfSpecies()}
      exist_spec_annotation_filt = {val:exist_spec_annotation_raw[val] for val in exist_spec_annotation_raw.keys() \
                                    if exist_spec_annotation_raw[val] is not None}
      species_names = {val.getId():val.name for val in self.model.getListOfSpecies()}
      species_tuple = (species_names, exist_spec_annotation_filt)
      self.species = sa.SpeciesAnnotation(inp_tuple=species_tuple)
      # Create reaction_annotation instance
      # Annotation of Rhea
      reac_dict_raw_rhea = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.RHEA) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_rhea = {k:reac_dict_raw_rhea[k] \
                                 for k in reac_dict_raw_rhea.keys() \
                                 if reac_dict_raw_rhea[k] is not None}
      reac_dict_format_rhea = {k:['RHEA:'+val for val in reac_dict_raw_filt_rhea[k]] \
                                 for k in reac_dict_raw_filt_rhea.keys()}
      reac_dict_rhea = dict()
      for one_id in reac_dict_format_rhea.keys():
        one_itm = list(set([ref_rhea2bi[val] for val in reac_dict_format_rhea[one_id] \
                   if val in ref_rhea2bi.keys()]))
        if len(one_itm) > 0:
          reac_dict_rhea[one_id] = one_itm
      # Annotation of KEGG (mapped to corresponding Rhea BI term) 
      reac_dict_raw_kegg = {r.getId():tools.getQualifierFromString(r.getAnnotationString(), cn.KEGG_REACTION) \
                           for r in self.model.getListOfReactions()}
      reac_dict_raw_filt_kegg = {k:reac_dict_raw_kegg[k] \
                                 for k in reac_dict_raw_kegg.keys() \
                                 if reac_dict_raw_kegg[k] is not None}
      reac_dict_kegg = {k:[ref_kegg2rhea_bi[val] \
                           for val in reac_dict_raw_filt_kegg[k] if val in ref_kegg2rhea_bi.keys()] \
                        for k in reac_dict_raw_filt_kegg.keys()}
      reac_exist_annotation = reac_dict_rhea
      for one_id in reac_dict_kegg.keys():
        if one_id in reac_exist_annotation.keys():
          reac_exist_annotation[one_id] = list(set(reac_exist_annotation[one_id] + reac_dict_kegg[one_id]))
        else:
          reac_exist_annotation[one_id] = list(set(reac_dict_kegg[one_id]))
      # Next, reaction components for each reaction
      reac_components = {val.getId():list(set([k.species for k in val.getListOfReactants()]+[k.species for k in val.getListOfProducts()])) \
                         for val in self.model.getListOfReactions()}
      reaction_tuple = (reac_components, reac_exist_annotation)
      self.reactions = ra.ReactionAnnotataion(inp_tuple=reaction_tuple)




    # species_annotation = sa.SpeciesAnnotation(libsbml_fpath)
    # reaction_anotation = ra.ReactionAnnotataion(libsbml_fpath)
    