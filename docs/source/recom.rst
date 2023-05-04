

Getting Recommendation
======================


This section explains more detailed usage of ``AMAS``. The ``recommend_annotation`` command presented in the :doc:`basics` section has several optional arguments and you can use them as below:

.. code-block:: console
 
   $ recommend_annotation BIOMD0000000190.xml --cutoff 0.9 --save sbml --outfile result.sbml
   ... 
   Analyzing 11 species...

   ...
   Analyzing 13 reaction(s)...

   Annotation recommended for 11 species:
   [A, AcCoA, CoA, D, Met, ORN, P, S, SAM, aD, aS]

   Annotation recommended for 2 reaction(s):
   [ODC, P_efflux]
  
   Recommendations saved as:
   /Users/amas/result.sbml


In this example, ``AMAS`` automatically detected all existing species and reactions and made predictions, but only recommended annotations for the elements with match score of 0.9 or higher (--cutoff option). In addition, by choosing `sbml` for the --save option, a new SBML model file with updated annotations was created and saved. The default value of --save is `csv`. 

There are two additional commands to get recommendations for species and reactions, respectively. ``recommend_species`` and ``recommend_reactions`` take similar arguments as that of the above command, but you can explicitly choose the elements to be recommended; in addition, you can set the minimum length of names (species) or the minimum number of components (reactions) to improve overall accuracy of the predictions. 


--To be added--
