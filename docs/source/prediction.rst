

Predicting Annotations
======================

Currently, ``AMAS`` uses ChEBI to annotate species and Rhea to annotate reactions. The algorithm predicts species annotations using their IDs or display names; the ChEBI therms with the smallest edit distance will be chosen.

-- placeholder of example--

Reactions are predicted based on the annotations of species. Annotations of the components of a reaction (that is, the chemical species consisting  reactants and products) are compared with that of the Rhea databases; the terms with the most matches will be chosen as candidates.

-- example --
