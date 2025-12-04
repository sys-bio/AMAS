#!/usr/bin/env python

# recommend_genes.py
"""
Predicts annotations of gene using a local SBML file
and the gene ID. 
Usage: python recommend_genes.py files/e_coli_core.xml --tax 511145 --min_len 2 --cutoff 0.6 --outfile genes_rec.csv
Usage (SBML-qual model): python recommend_genes.py files/111_APOPTOSIS_source.sbml
"""

import argparse
import os
from os.path import dirname, abspath
import pandas as pd
import sys
sys.path.insert(0, dirname(dirname(abspath(__file__))))

from AMAS import constants as cn
from AMAS import recommender


def main():
  parser = argparse.ArgumentParser(description='Recommend gene annotations of an SBML model and save results') 
  parser.add_argument('model', type=str, help='SBML model file (.xml or .sbml)')
  parser.add_argument('--tax', type=str, help='Taxonomy ID of the species', default='9606')
  # One or more gene IDs can be given
  parser.add_argument('--genes', type=str, help='ID(s) of genes to be recommended. ' +\
                                                  'If not provided, all genes will be used', nargs='*')
  parser.add_argument('--min_len', type=int, help='Minimum length of gene names to be used for prediction. ' +\
                                                  'gene with names that are at least as long as this value ' +\
                                                  'will be analyzed. Default is zero', nargs='?', default=0)
  parser.add_argument('--cutoff', type=float, help='Match score cutoff', nargs='?', default=0.0)
  parser.add_argument('--mssc', type=str,
                                help='Match score selection criteria (MSSC). ' +\
                                     'Choose either "top" or "above". "top" recommends ' +\
                                     'the best candidates that are above the cutoff, ' +\
                                     'and "above" recommends all candidates that are above ' +\
                                     'the cutoff. Default is "top"',
                                nargs='?',
                                default='top')
  parser.add_argument('--outfile', type=str, help='File path to save recommendation.', nargs='?',
                      default=os.path.join(os.getcwd(), 'genes_rec.csv'))
  args = parser.parse_args()

  try:
      recom = recommender.Recommender(libsbml_fpath=args.model)
  except ValueError as e:
      print(f"Error loading model: {e}")
      sys.exit(1)

  genes = args.genes
  min_len = args.min_len
  cutoff = args.cutoff
  mssc = args.mssc.lower()
  outfile = args.outfile
  if args.tax in cn.REF_TAXID2LABEL.keys():
    tax_id = args.tax
    tax_name = cn.REF_TAXID2LABEL[tax_id]
  else:
    print(f"Taxonomy ID {args.tax} not found in the reference data.")
    sys.exit(1)

  # if nothing is given, predict all IDs
  if genes is None:
    genes = recom.getGeneIDs()
    if not genes:
      print("No genes found in the model.")
      sys.exit(1)
  if recom.model_type == 'SBML-qual':
    print("...\nAnnotating SBML-qual model\n")
  print("...\nRunning for organism: %s (tax_id: %s)\n" % (tax_name, tax_id))
  print("...\nAnalyzing %d genes...\n" % len(genes))

  if recom.model_type == 'SBML':
    res_tab = recom.recommendGene(tax_id=tax_id,
                                  ids=genes,
                                  mssc=mssc,
                                  cutoff=cutoff,
                                  min_len=min_len,
                                  outtype='table')
  else:
    res_tab = recom.recommendQualitativeSpecies(tax_id=tax_id,
                                  ids=genes,
                                  mssc=mssc,
                                  cutoff=cutoff,
                                  min_len=min_len,
                                  outtype='table')
  recom.saveToCSV(res_tab, outfile)
  if isinstance(res_tab, pd.DataFrame):
      print("Recommendations saved as:\n%s\n" % os.path.abspath(outfile))


if __name__ == '__main__':
  main()