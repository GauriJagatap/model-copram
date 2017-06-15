#copram

Code for phase retrieval of signals with underlying structured sparsity level s (Compressive Phase Retrieval with Alternating Minimization, i.e. CoPRAM). Refer the paper:

"Phase Retrieval Using Structured Sparsity: A Sample Efficient Algorithmic Framework", by Gauri Jagatap and Chinmay Hegde,

for further details on implementation: https://arxiv.org/abs/1705.06412

phase_retrieval_sparse.m 

Runs and analyzes performances of the following sparse phase retrieval algorithms:
1. CoPRAM
2. AltMinSparse 
(implemented based on the paper https://arxiv.org/abs/1306.0160)
3. Thresholded Wirtinger Flow for sparse phase retrieval (ThWF)
(implemented based on the paper https://arxiv.org/abs/1506.03382)
4. Sparse Phase Retrieval using Truncated Amplitude Flow (SparTA)
(implemented based on the paper https://arxiv.org/abs/1611.07641)

# model-copram

Phase retrieval of block-sparse signals of uniform block length b and overall sparsity s (Block Compressive Phase Retrieval with Alternating Minimization, i.e. Block CoPRAM). Refer the paper:

"Phase Retrieval Using Structured Sparsity: A Sample Efficient Algorithmic Framework", by Gauri Jagatap and Chinmay Hegde,

for further details on implementation: https://arxiv.org/abs/1705.06412

phase_retrieval_block.m

Runs and analyzes performance of Block CoPRAM algorithm w.r.t CoPRAM for block sparse signals.

