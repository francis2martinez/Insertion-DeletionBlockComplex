# Finding Minimal 2 Cycles
The accompanying Jupyter notebook contains the computations necessary to prove that any set of words whose insertion chain complex has nontrivial 2-homology must consist of at least 8 words. The code systematically eliminates all smaller complexes by generating every possible 2-dimensional skeleton whose 2-cells are squares and 1-cells are edges, with the boundary operator defined as in our paper. Leveraging the results of Lemmas A.7 and A.8, each 2-skeleton is algorithmically analyzed and, in the final cases, verified by hand. This exhaustive process establishes Theorem 5.7.

Our algorithm,at a high level, proceeds as follows (this explanation is also provided in the manuscript):

1. We generate all possible 1-skeletons (graphs) with $n=5,6,7$ vertices using the function `generate\_1\_skeletons(n)`. 
2. We generate all possible 2-skeletons with the given list of possible 1-skeletons (`Possible\_Graphs`) that have nontrivial second homology, and that don't contain any of the forbidden patterns from Lemma A.8, using the function `generate\_2\_skeletons(Possible\_Graphs, forbidden\_patterns)`. This function only returns the corresponding sub-complexes with non-zero second homology, as pairs of the form `(Graph, squares)`, where `Graph` is the 1-skeleton, and `squares` is the list of 2-faces. 
3. For the 2-complexes that still remain after the previous two steps, we count how many different symbols must be inserted in the squares using the function `count\_inserted\_symbols(G,squares)` and discard those that require only one single symbol inserted, since the desired complex cannot be obtained in this way. 
4. Finally, 5 possible 2-skeletons survive, but applying Lemma A.7 to the pairs of squares sharing 2 edges, ultimately shows that these are either impossible or that they have vanishing 2-homology. 

