# AMAS
AMAS (Automatic Model Annotation System) predicts and recommends annotations of systems biology models. Current version focuses on species and reactions of models in SBML format, which can be tested using BioModels or BiGG models. 

## Overview
``AMAS`` is a collection of methods to predict and recommend annotations for SBML model elements. Current version focuses on predicting species and reaction annotations of metabolic models, such as that can be found in BiGG and BioModels repositories. Algorithm uses CHEBI for species annotations, and Rhea for reaction annotations.

## Example
First, class instance should be created using ``AMAS.recommender.Recommender`` with an existing SBML file (optional). 
<img src="https://github.com/woosubs/AMAS/raw/main/png/create_instance.png" width="800"/>

Next, to get recommendation for species, user can use `.getSpeciesAnnotation` method. 
<img src="https://github.com/woosubs/AMAS/raw/main/png/getspecies_annotation.png" width="800"/>

Result is a namedtuple "Recommendation", with attributes including id, credibility, recommended CHEBI terms, and the urls of such terms. 
<img src="https://github.com/woosubs/AMAS/raw/main/png/spec_recommendation.png" width="800"/>

