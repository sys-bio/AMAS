

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


In this example, ``AMAS`` automatically detected all existing species and reactions and made predictions, but only recommended annotations for the elements with match score of 0.9 or higher (``cutoff`` option). In addition, by choosing *sbml* for the ``save`` option, a new SBML model file with updated annotations was created and saved. The default value of ``save`` is *csv*. 

You can also choose to optimize predictions using the ``optimize`` option:

.. code-block:: console
 
   $ recommend_annotation BIOMD0000000015.xml --optimize y --outfile opt.csv
   ... 
   Analyzing 18 species...

   ...
   Analyzing 37 reaction(s)...

   Optimizing predictions...

   Annotation recommended for 18 species:
   [ATP, Ade, DNA, GTP, Gua, HX, IMP, PRPP, Pi, 
    R5P, RNA, SAM, SAMP, UA, XMP, Xa, dATP, dGTP]

   Annotation recommended for 37 reaction(s):
   [ada, ade, adna, adrnr, ampd, aprt, arna, asli,
    asuc, dada, den, dgnuc, dnaa, dnag, gdna, gdrnr,
    gmpr, gmps, gnuc, gprt, grna, gua, hprt, hx, hxd,
    impd, inuc, mat, polyam, prpps, pyr, rnaa, rnag,
    trans, ua, x, xd]
  
   Recommendations saved as:
   /Users/amas/opt.csv

In the above example, ``AMAS`` compared predictions of species and reactions and updated recommendations of them based on the comparison. Recommendations of species such as *IMP* and that of reactions such as *dada* have been updated. Note that it might take a significant amount of time if the numbers of species and reactions are large. 


There are two additional commands to get recommendations for species and reactions, respectively. ``recommend_species`` and ``recommend_reactions`` take similar arguments as that of the above command, but you can explicitly choose the elements to be recommended; in addition, you can set the minimum length of names (species) or the minimum number of components (reactions) to improve overall accuracy of the predictions. The example below shows how these arguments are used:


.. code-block:: console
 
   $ recommend_species BIOMD0000000190.xml --species SAM CoA --min_len 5 --cutoff 0.9
   ...
   Analyzing 2 species...

   Annotation recommended for 1 species:
   [SAM]

   Recommendations saved as:
   /Users/amas/species_rec.csv
  

In the example above, user asked recommendations only for two species, SAM and CoA, and chose to save recommendations if their query name was at least 5 and match score at least 0.9. The species CoA did not meet all of the criteria, so only annotations of SAM were recommended and saved. Since the path of the output file was not specified, recommendation was saved as ``species_rec.csv`` (default file name) in the current working directory. 

Getting recommendation for reactions using the ``recommend_reactions`` command is similar:


.. code-block:: console
 
   $ recommend_reactions BIOMD0000000190.xml --cutoff 0.5 --mssc above --outfile reactions.csv
   ...
   Analyzing 13 reaction(s)...

   Annotation recommended for 10 reaction(s):
   [MAT, ODC, P_efflux, SAMdc, SSAT_for_D, SSAT_for_S, SpdS, SpmS, VCoA, VacCoA]

   Recommendations saved as:
   /Users/amas/reactions.csv


This time, no reaction ID was listed; thus, ``AMAS`` will detect all existing reactions and make recommendations for those with match score of 0.5 or above. By choosing *above* for the ``mssc`` option, ``AMAS`` will recommend all of the predicted candidates with match score at or above the cutoff. If *top* (default value) was chosen instead, ``AMAS`` would report only those with the highest match score that is at or above the cutoff. 
