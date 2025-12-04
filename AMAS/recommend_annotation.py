#!/usr/bin/env python

# recommend_annotation.py
"""
Predicts annotations of species, reactions, and genes (if exist) using a local XML file
This is a combined version of recommend_species, recommend_reaction, and recommend_genes,
but is more convenient because user will just get the updated XML file or whole recommendations. 

Usage: python recommend_annotation.py files/e_coli_core.xml --tax 511145 --cutoff 0.6 --save csv --outfile rec.csv 
"""

import argparse
import itertools
import libsbml
import numpy as np
import os
from os.path import dirname, abspath
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import iterator as it
from AMAS import species_annotation as sa
from AMAS import reaction_annotation as ra
from AMAS import recommender
from AMAS import tools

def main():
  parser = argparse.ArgumentParser(description='Recommend annotations of an SBML model ' +\
                                               '(for species, reactions, and genes) and save results.') 
  parser.add_argument('model', type=str, help='SBML model file (.xml).')
  # taxonomy id for gene annotation, by default is human
  parser.add_argument('--tax', type=str, help='Taxonomy ID of the species', nargs='?', default='9606') 
  parser.add_argument('--cutoff', type=float, help='Match score cutoff.', nargs='?', default=0.0)
  parser.add_argument('--optimize', type=str, help='Whether to optimize or not. ' +\
                                                   'If y or yes is given, predictions will be ' +\
                                                   'optimized. N or no will not optimize predictions.',
                                              nargs='?',
                                              default='no')
  parser.add_argument('--mssc', type=str,
                                help='Match score selection criteria (MSSC). ' +\
                                     'Choose either "top" or "above". "top" recommends ' +\
                                     'the best candidates that are above the cutoff, ' +\
                                     'and "above" recommends all candidates that are above ' +\
                                     'the cutoff. Default is "top"',
                                nargs='?',
                                default='top')
  parser.add_argument('--save', type=str, 
                                help='Either "sbml" or "csv". ' +\
                                     'If "sbml" is chosen, model will be automatically ' +\
                                     'annotated with recommended candidates and saved. ' +\
                                     'If "csv" is chosen, recommendations will be saved ' +\
                                     'as a csv file. Default is "sbml".',
                                nargs='?',
                                default='sbml')
  parser.add_argument('--outfile', type=str, help='Path to save an output file.', nargs='?')
  args = parser.parse_args()
  one_fpath = args.model
  tax = str(args.tax)
  cutoff = args.cutoff
  optim_raw = args.optimize
  if optim_raw.lower() in ['y', 'yes']:
    optim = True
  else:
    optim = False
  mssc = args.mssc.lower()
  save = args.save
  outfile = args.outfile

  recom = recommender.Recommender(libsbml_fpath=one_fpath)
  if recom.model_type == 'SBML-qual':
    print("...\nAnnotating SBML-qual model\n")
    qual_species = recom.getQualSpeciesIDs()
    if not qual_species:
      print("\nNo qualitative species found in the model.")
    else:
      print("...\nAnalyzing %d qualitative species...\n" % len(qual_species))
    res_tab = recom.recommendQualitativeSpecies(tax_id=tax,
                                  ids=qual_species,
                                  mssc=mssc,
                                  cutoff=cutoff,
                                  outtype='table')
  
  else:
    specs = recom.getSpeciesIDs()
    if not specs:
      print("\nNo species found in the model.")
    else:
      print("...\nAnalyzing %d species...\n" % len(specs))
    reacts = recom.getReactionIDs()
    if not reacts:
      print("\nNo reactions found in the model.")
    else:
      print("...\nAnalyzing %d reaction(s)...\n" % len(reacts))
    genes = recom.getGeneIDs()
    if not genes:
      print("\nNo genes found in the model.")
    else:
      print("...\nAnalyzing %d gene(s)...\n" % len(genes))

    res_tab = recom.recommendAnnotation(tax_id = tax,
                                        mssc=mssc,
                                        cutoff=cutoff,
                                        optimize=optim,
                                        outtype='table')
  if save == 'csv':
    if outfile is None:
      outfile = os.path.join(os.getcwd(), 'recommendations.csv')
    recom.saveToCSV(res_tab, outfile)
  else:
    if outfile is None:
      outfile = os.path.join(os.getcwd(), 'updated_model.xml')
    res_sbml = recom.getSBMLDocument(sbml_document=recom.sbml_document,
                                     chosen=res_tab,
                                     auto_feedback=True)
    libsbml.writeSBMLToFile(res_sbml, outfile)   
  
  print("Recommendations saved as:\n%s\n" % os.path.abspath(outfile))


if __name__ == '__main__':
  main()