

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


In this example, ``AMAS`` automatically detected all existing species and reactions and made predictions, but only recommended annotations for the elements with match score of 0.9 or higher (``cutoff`` option). In addition, by choosing *sbml* for the ``save`` option, a new SBML model file with updated annotations was created and saved. The default value of ``save`` is *csv*, which will create a comma-separated value (csv) file. 

Here, we explain what the term *match score* means. In ``AMAS``, match score represents the measure to compute the similarity between the information from species/reactions and that from the databases of annotations, such as ChEBI and Rhea. For species, the match score represents cosine similarity between vector-based representations of the query and the reference species. For reactions, it is the number of overlapping species between the query and the reference reactions, normalized to the minimum number of species in those reference reactions that show the largest overlap with the query reaction. In short, ``AMAS`` tries to sort the possible candidates and tries to recommend most likely annotations for the user. 

You can also choose to optimize predictions using the ``optimize`` option. When this option is incurred, ``AMAS`` compares once-predicted annotations of species and reactions and iteratively updates them. To be more specific, ``AMAS`` tries to match components of predicted Rhea annotations with predicted species annotation, and if there is an unmatched species, it tries to replace its annotation with annotation from Rhea. The update will be accepted if the newly calculated match score improves. The example below illustrates how one can use this option:

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
  

In the example above, user asked recommendations only for two species, *SAM* and *CoA*, and chose to save recommendations if their query name was at least 5 and match score at least 0.9. The species *CoA* did not meet all of the criteria, so only annotations of SAM were recommended and saved. Since the path of the output file was not specified, recommendation was saved as ``species_rec.csv`` (default file name) in the current working directory. 

Getting recommendation for reactions using the ``recommend_reactions`` command is similar:


.. code-block:: console
 
   $ recommend_reactions BIOMD0000000190.xml --cutoff 0.5 --mssc above --outfile reactions.csv
   ...
   Analyzing 13 reaction(s)...

   Annotation recommended for 10 reaction(s):
   [MAT, ODC, P_efflux, SAMdc, SSAT_for_D, SSAT_for_S, SpdS, SpmS, VCoA, VacCoA]

   Recommendations saved as:
   /Users/amas/reactions.csv


This time, no reaction ID was listed; thus, ``AMAS`` will detect all existing reactions and make recommendations for those with match score of 0.5 or above. ``mssc`` means Match Score Selection Criteria, which helps the algorithm make automatic selection based on the match scores computed for all possible candidates. There are two options: *top* and *above*. By choosing *above* for the ``mssc`` option, ``AMAS`` will recommend all of the predicted candidates with match score at or above the cutoff. If *top* (default value) was chosen instead, ``AMAS`` would report only those with the highest match score that is at or above the cutoff. 
