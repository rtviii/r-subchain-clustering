#### So it goes:

1. Get your ```.cif``` model ready.
2. Pipe your rnas from post-pfam *pdbid* profile and construct a neighbour tree for a given **GLOBAL PARAMETER** *radius* excluding those *rna*s.
3. Reduce neighbour tree to clusters of neighbours. **HAS TO BE MADE RECURSIVE**.
4. Save to file.
5. Iterate (1-4) for radii in range 1.5-8 angstrom to see full range of clustering.
6. 