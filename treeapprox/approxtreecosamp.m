function chat_tree = approxtreecosamp(recflag, b, k, d, A, At, tol, maxIter, VERBOSE, tol_early)

%%%%%%%%%%%
% 


temp = At(b);
N = length(temp);
err = 1; 
aa= zeros(N,1); 
s_cosamp = zeros(N,1);

verbose = false;
if VERBOSE
    verbose = true;
end

iterCnt = 0;
numAtoms = 2;

while ((err>tol)&&(iterCnt<=maxIter))
    iterCnt = iterCnt + 1;
    
    rr = b - A(s_cosamp);
    proxy = At(rr); %proxy = proxy(:);
    %------Estimation
    switch recflag
        case 'sparse'
            supp = thresh(proxy,numAtoms*k);
        case 'tree_approx' 
            supp = tree_approx(proxy,numAtoms*k,numAtoms*k*1.1,d,VERBOSE);
            supp = supp(:);
        case 'tree_exact' 
            supp = exact_tree_proj(proxy,d,numAtoms*k);
            supp = abs(supp) > 0;
        case 'tree_approx_mex'
            proxy = proxy(:); proxy = proxy';
            supp = tree_approx_cpp_mex(proxy,numAtoms*k,d);
            supp = supp(:);
        case 'tree_approx_mex_ludwig'
            opts.verbose = verbose;
            opts.tree_layout = 'wavelet_tree';
            supp = treeapprox_binsearch(proxy(:).^2,d,numAtoms*k,numAtoms*round(1.1*k),opts);
        case 'tree_exact_mex_ludwig'
            supp = treeexact_smalltable_wvtree(proxy(:).^2,d,numAtoms*k);
        case 'tree_exact_fulltable_mex_ludwig'
            supp = treeexact_fulltable_wvtree(proxy(:).^2,d,numAtoms*k);
    end
    tt = union(find(ne(s_cosamp,0)),find(supp));
    
    %------Least-squares
    PP_tt = @(z) A_I(A,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
    qq = PP_transpose_tt(b);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    w = cgsolve(PPtranspose_PP_tt,qq,1e-6, 100, 0);
    bb= 0*s_cosamp; bb(tt)= w;
    
    %------Prune
    switch recflag
        case 'sparse'
            supp = thresh(bb,k);
        case 'tree_approx'
            supp = tree_approx(bb,k,k*1.1,d,VERBOSE);
            supp = supp(:);
        case 'tree_exact' 
            supp = exact_tree_proj(bb,d,k);
            supp = abs(supp) > 0;
        case 'tree_approx_mex'
            bb = bb(:); bb = bb';
            supp = tree_approx_cpp_mex(bb,k,d);
            supp = supp(:);
        case 'tree_approx_mex_ludwig'
            opts.verbose = verbose;
            opts.tree_layout = 'wavelet_tree';
            supp = treeapprox_binsearch(bb(:).^2,d,k,round(1.1*k),opts);
        case 'tree_exact_mex_ludwig'
            supp = treeexact_smalltable_wvtree(bb(:).^2,d,k);
        case 'tree_exact_fulltable_mex_ludwig'
            supp = treeexact_fulltable_wvtree(bb(:).^2,d,k);
    end
    bb=bb(:);
    s_cosamp = supp.*bb;
    
    aa = s_cosamp;
    
    errTemp = norm(A(s_cosamp)-b)/norm(b);
    if (VERBOSE)
        fprintf('Iter: %d Err: %f\n',iterCnt,errTemp);
    end
    
    if ((err-errTemp)/err <= tol_early) %Early termination condition
        err = errTemp;
        if (VERBOSE)
            fprintf('Terminating early.\n');
        end
        break 
    end
    err = errTemp;
end

chat_tree= aa;

end

function supp = thresh(x,k)
[~,ww] = sort(abs(x(:)),'descend');
supp = 0*x;
supp(ww(1:k)) = 1;
end