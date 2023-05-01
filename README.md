# AMAS: Automatic Model Annotation System
``AMAS`` (Automatic Model Annotation System) predicts and recommends annotations of systems biology models. It uses information of species and reactions of SBML models and provides recommendations. ``AMAS`` can be tested using BioModels or BiGG models. 

## Overview
``AMAS`` is a collection of methods to predict and recommend annotations for SBML model elements. Current version focuses on predicting species and reaction annotations of metabolic models, such as that can be found in BiGG and BioModels repositories. The algorithm uses CHEBI for species annotations, and Rhea for reaction annotations. 

For detailed information, please check our
[Read the Docs](https://amas.readthedocs.io/en/latest/
) page.

## Example
First, class instance should be created using ``AMAS.recommender.Recommender`` with an existing SBML file (optional). 
<img src="https://github.com/woosubs/AMAS/raw/main/png/create_instance.png" width="800"/>

Next, to get recommendation for species, user can use the `.getSpeciesAnnotation` method. 

<img src="https://github.com/woosubs/AMAS/raw/main/png/getspecies_annotation_id.png" width="800"/>

When the `pred_id` argument is used, Recommender will search the model to find an available display name for prediction. Alternatively, user can use the `pred_str` argument to make a direct prediction, which does not need a pre-loaded model in the constructor. Below, 'S-adenosyl-L-methionine' is the display name of the species 'SAM' in the model file, so the result will be the same. 

<img src="https://github.com/woosubs/AMAS/raw/main/png/getspecies_annotation_str.png" width="800"/>


Result is a namedtuple 'Recommendation', with attributes including id, credibility, recommended CHEBI terms, and the urls of such terms. 
<img src="https://github.com/woosubs/AMAS/raw/main/png/spec_recommendation.png" width="800"/>

Similarly, recommendation of a reaction can be also obtained using the `.getReactionAnnotation` method. 

<img src="https://github.com/woosubs/AMAS/raw/main/png/getreaction_annotation.png" width="800"/>

Recommendation of a reaction uses Rhea database. 
<img src="https://github.com/woosubs/AMAS/raw/main/png/reac_recommendation.png" width="800"/>
