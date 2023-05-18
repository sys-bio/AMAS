

Command-line Arguments
======================


This section summarizes arguments for the commands described in the other sections. 

.. list-table:: Arguments for ``recommend_annotation``
   :widths: 35 50 50 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \\-\\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \\-\\-method
     - string (*top* or *above*)
     - mode of selection
     - *top*
   * - \\-\\-save
     - string (*sbml* or *csv*)
     - type of file to be saved
     - *csv*
   * - \\-\\-outfile
     - string 
     - path to save file
     - upated model.xml / recommendations.csv


.. list-table:: Arguments for ``recommend_species``
   :widths: 35 50 50 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \\-\\-species
     - string (one or more)
     - list of species IDs
     - all existing species
   * - \\-\\-reject
     - integer
     - element filter
     - 0
   * - \\-\\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \\-\\-method
     - string (*top* or *above*)
     - mode of selection
     - *top*
   * - \\-\\-outfile
     - string 
     - path to save file
     - species_rec.csv


.. list-table:: Arguments for ``recommend_reactions``
   :widths: 35 50 50 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - model
     - string
     - SBML model file
     - N/A
   * - \\-\\-reactions
     - string (one or more)
     - list of reaction IDs
     - all existing reactions
   * - \\-\\-reject
     - integer
     - element filter
     - 0
   * - \\-\\-cutoff
     - float (0.0 - 1.0)
     - match score cutoff
     - 0.0
   * - \\-\\-method
     - string (*top* or *above*)
     - mode of selection
     - *top*
   * - \\-\\-outfile
     - string 
     - path to save file
     - reaction_rec.csv


.. list-table:: Arguments for ``update_annotation``
   :widths: 35 50 50 50 
   :header-rows: 1

   * - Name
     - Type
     - Description
     - Default value
   * - infile
     - string
     - path of the original model file
     - N/A
   * - csv_select
     - string
     - file with feedback (*UPDATE ANNOTATION*)
     - N/A
   * - outfile
     - string
     - path of the new file with annotation
     - N/A

