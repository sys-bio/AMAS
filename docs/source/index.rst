.. AMAS documentation master file, created by
   sphinx-quickstart on Sun Nov  6 20:01:20 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AMAS: Automatic Model Annotation System
=======================================

``AMAS`` (Automatic Model Annotation System) predicts and recommends annotations of systems biology models. The current version of ``AMAS`` recommends species and reaction annotations in SBML models. The algorithm was developed and tested using curated models from `BioModels <https://www.ebi.ac.uk/biomodels/>`_ and `BiGG <http://bigg.ucsd.edu/>`_.


Overview
========
``AMAS`` provides a collection of command-line methods to obtain annotations and update SBML model files. There are four commands: ``recommend_annotation``, ``recommend_speces``, ``recommend_reactions``, and ``update_annotation``. Annotations are recommended in `ChEBI <https://www.ebi.ac.uk/chebi/>`_ terms for species and in `Rhea <https://www.rhea-db.org/>`_ terms for reactions. 

Please check out the contents below for more information on the commands and how to use them. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   basics
   recom
   upd
   tbls

