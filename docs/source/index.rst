.. AMAS documentation master file, created by
   sphinx-quickstart on Sun Nov  6 20:01:20 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AMAS: Automatic Model Annotation System
=======================================

``AMAS`` (Automatic Model Annotation System) predicts and recommends annotations of systems biology models. Current version predicts annotations of species and reactions in SBML models, and the algorithm can be tested using BioModels or BiGG models.

Check out the :doc:`prediction` section for more information. 


Overview
========
AMAS is a collection of methods to predict and recommend annotations for SBML model elements. Current version focuses on predicting species and reaction annotations of metabolic models, such as that can be found in BiGG and BioModels repositories. Algorithm uses CHEBI for species annotations, and Rhea for reaction annotations.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   recommendation
