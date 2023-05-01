

Getting Recommendations
=======================


This section describes how ``AMAS`` can be used to get recommendations of annotations. 

To install the package, use:
      ``pip install AMAS-sb``

Currently, ``AMAS`` can run either as a command-line tool or as a regular Python package. Most of the time, you will need to provide an exsiting SBML model file (XML) to run it. 

Once installed, you can run the command ``recommend_annotation`` to get recommendations for all species and reactions existing in the model file. The example below assumes a model file ``BIOMD0000000190.xml`` is in the current working directory and uses it. 

.. code-block:: console
 
   $ recommend_annotation BIOMD0000000190.xml --outfile res.csv
   ...
   Analyzing 11 species...

   ...
   Analyzing 13 reaction(s)...

   Annotation recommended for 11 species:
   [A, AcCoA, CoA, D, Met, ORN, P, S, SAM, aD, aS]

   Annotation recommended for 13 reaction(s):
   [MAT, ODC, PAO_for_aD, PAO_for_aS, P_efflux, SAMdc, SSAT_for_D, SSAT_for_S, SpdS, SpmS, VCoA, VacCoA, aD_efflux]
  
   Recommendations saved as:
   /Users/amas/rec.csv


Recommendations will be made using ChEBI identifiers to annotate species and Rhea to annotate reactions. The algorithm predicts species annotations using their IDs or display names; the ChEBI therms with the smallest edit distance will be chosen.

-- placeholder of example--

Reactions are predicted based on the annotations of species. Annotations of the components of a reaction (that is, the chemical species consisting  reactants and products) are compared with that of the Rhea databases; the terms with the most matches will be chosen as candidates.

-- example --
