

Getting Recommendations
=======================


This section describes how ``AMAS`` can be used to get recommendations of annotations. 

To install the package, use:
      ``pip install AMAS-sb``

Currently, ``AMAS`` can run either as a command-line tool or as a regular Python package. Most of the time, you will need to provide an exsiting SBML model file (XML) to run it. 

Once installed, you can run the command ``recommend_annotation`` to get recommendations for all species and reactions existing in the model file. The example below uses a model file 


Recommendations will be made using ChEBI identifiers to annotate species and Rhea to annotate reactions. The algorithm predicts species annotations using their IDs or display names; the ChEBI therms with the smallest edit distance will be chosen.

-- placeholder of example--

Reactions are predicted based on the annotations of species. Annotations of the components of a reaction (that is, the chemical species consisting  reactants and products) are compared with that of the Rhea databases; the terms with the most matches will be chosen as candidates.

-- example --
