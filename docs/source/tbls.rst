

Command-line Arguments
======================


This section summarizes arguments for the commands described in the other sections. Arguments with two dashes (``--``) indicate that they are optional arguments. All arguments of ``update_annotation`` are positional (i.e., required). 

.. list-table:: Arguments for ``recommend_annotation``
   :widths: 35 50 70 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \-\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \-\-optimize
     - string (*y, yes*)
     - optimizing predictions
     - *no*
   * - \-\-mssc
     - string (*top* or *above*)
     - match score selection criteria
     - *top*
   * - \-\-save
     - string (*sbml* or *csv*)
     - type of file to be saved
     - *csv*
   * - \-\-outfile
     - string 
     - path to save file
     - *upated model.xml* / *recommendations.csv*


.. list-table:: Arguments for ``recommend_species``
   :widths: 35 50 70 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \-\-species
     - string (one or more)
     - list of species IDs
     - all existing species
   * - \-\-min_len
     - integer
     - element filter (minimum length of name)
     - 0
   * - \-\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \-\-mssc
     - string (*top* or *above*)
     - match score selection criteria
     - *top*
   * - \-\-outfile
     - string 
     - path to save file
     - *species_rec.csv*


.. list-table:: Arguments for ``recommend_reactions``
   :widths: 35 50 70 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \-\-reactions
     - string (one or more)
     - list of reaction IDs
     - all existing reactions
   * - \-\-min_len
     - integer
     - element filter (minimum number of components)
     - 0
   * - \-\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \-\-mssc
     - string (*top* or *above*)
     - match score selection criteria
     - *top*
   * - \-\-outfile
     - string 
     - path to save file
     - *reaction_rec.csv*


.. list-table:: Arguments for ``update_annotation``
   :widths: 35 50 70 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - infile
     - string
     - path of the original model file
     - N/A
   * - feedback
     - string
     - file with feedback (*UPDATE ANNOTATION* column)
     - N/A
   * - outfile
     - string
     - path of the new file with updated annotations
     - N/A

