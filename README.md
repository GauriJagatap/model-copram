# copram

Code for phase retrieval of signals with underlying structured sparsity level s (Compressive Phase Retrieval with Alternating Minimization, i.e. CoPRAM). Please cite the following paper:

Gauri Jagatap and Chinmay Hegde, "Fast, Sample-Efficient Algorithms for Structured Phase Retrieval." In Advances in Neural Information Processing Systems 2017 (pp. 4924-4934).
Full paper at: http://papers.nips.cc/paper/7077-fast-sample-efficient-algorithms-for-structured-phase-retrieval

For further details on implementation refer the arXiv version at: https://arxiv.org/abs/1705.06412


Main code:
phase_retrieval_sparse.m 

Runs and analyzes performances of the following sparse phase retrieval algorithms:
1. CoPRAM
2. AltMinSparse 
(implemented based on the paper https://arxiv.org/abs/1306.0160)
3. Thresholded Wirtinger Flow for sparse phase retrieval (ThWF)
(implemented based on the paper https://arxiv.org/abs/1506.03382)
4. Sparse Phase Retrieval using Truncated Amplitude Flow (SparTA)
(implemented based on the paper https://arxiv.org/abs/1611.07641)

# block-copram

Phase retrieval of block-sparse signals of uniform block length b and overall sparsity s (Block Compressive Phase Retrieval with Alternating Minimization, i.e. Block CoPRAM). Please cite the following paper:

Gauri Jagatap and Chinmay Hegde, "Fast, Sample-Efficient Algorithms for Structured Phase Retrieval." In Advances in Neural Information Processing Systems 2017 (pp. 4924-4934).
Full paper at: http://papers.nips.cc/paper/7077-fast-sample-efficient-algorithms-for-structured-phase-retrieval

For further details on implementation refer the arXiv version at: https://arxiv.org/abs/1705.06412

Main code:
phase_retrieval_block.m

Runs and analyzes performance of Block CoPRAM algorithm w.r.t CoPRAM for block sparse signals.

# tree-copram

Phase retrieval of tree-sparse signals of sparsity s (Tree Compressive Phase Retrieval with Alternating Minimization, i.e. Tree CoPRAM). Please cite the following paper:

G. Jagatap, C. Hegde, “Towards Sample-Optimal Methods for Solving Random Quadratic Equations with Structure”, IEEE International Symposium on Information Theory (ISIT), 2018
Full paper at: https://gaurijagatap.github.io/assets/ISIT18.pdf

Main code:
phase_retrieval_trees.m

Runs and analyzes performance of Tree CoPRAM algorithm w.r.t CoPRAM for tree sparse signals.



