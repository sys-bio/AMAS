# AMAS
AMAS (Automatic Model Annotation System) predicts and recommends annotations of systems biology models. Current version focuses on species and reactions of models in SBML format, which can be tested using BioModels or BiGG models. 

## Overview
``AMAS`` is a collection of methods to predict and recommend annotations for SBML model elements. Current version focuses on predicting species and reaction annotations of metabolic models, such as that can be found in BiGG and BioModels repositories. Algorithm uses CHEBI for species annotations, and Rhea for reaction annotations.

## Example
First, class instance should be created using ``AMAS.recommender.Recommender``.

