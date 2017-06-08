% cosamp.m
% Performs CS reconstruction using the CoSaMP algorithm
% (D. Needell and J. Tropp, "CoSaMP: iterative signal recovery from 
%  incomplete and inaccurate measurements" , ACHA, 2009)
% 
% INPUTS
% yy : measurement (M x 1)
% Phi: measurement matrix (M x N)
% K  : signal sparsity
% Its: number of iterations
% 
% OUTPUTS
% xhat   : Signal estimate (N x 1) 
% xcosamp: Matrix with N rows and at most Its columns; 
%          columns represent intermediate signal estimates   
% 
%
% CITE: Richard Baraniuk, Volkan Cevher, Marco Duarte, Chinmay Hegde
%       "Model-based compressive    sensing", submitted to IEEE IT, 2008.
% Created: Aug 2009.
% email: chinmay@rice.edu

function [xhat,xcosamp] = cosamp(yy, Phi, K, Its, s_init);

yy = yy(:); % 
[M,N] = size(Phi);

xcosamp = zeros(N,Its);
kk=1; 
maxiter= 100;
verbose= 0;
tol= 1e-4;

if nargin < 5
    s_cosamp = zeros(N,1);
else 
    s_cosamp = s_init;
end

while le(kk,Its),
    
    %-----Backprojection---%
    rcosamp = yy - Phi*(s_cosamp); %%updated residue
    proxy_cosamp = Phi'*(rcosamp); %%calculates proxy of signal
    %%(with good probability, corresponds to highest active columns)
    [trash,ww]= sort(abs(proxy_cosamp),'descend'); %% sort proxy vector in descending order 
    tt_cosamp= union(find(ne(s_cosamp,0)),ww(1:(2*K)));%%pick out indices corresponding to highest '2s' coefficients
    %%append to pre-existing support (will make the new support contain <= 3s terms = s + 2s)
    %------Estimate------%
    [w_cosamp] =  Phi(:,tt_cosamp)\yy;
%     [w_cosamp, res, iter] = cgsolve(Phi(:,tt_cosamp)'*Phi(:,tt_cosamp), Phi(:,tt_cosamp)'*yy,...
%                                         tol,maxiter, verbose); %%Use phi corresponding to new support and calculate
                                    %%support of approximation vector using
                                    %%cgsolve for resulting tall matrix (m x 3s) where m > 3s
    
    bb2= zeros(N,1); 
    bb2(tt_cosamp)= w_cosamp; %%store resultant approximation with at most 3s non-zero terms
    
    %---Prune----%
    kk = kk+1;   
    [trash,ww2]= sort(abs(bb2),'descend'); s_cosamp=0*bb2; %%prune signal to s terms where s is sparsity
    s_cosamp(ww2(1:K))= bb2(ww2(1:K));
    
    tt = ww2(1:K);

    xcosamp(:,kk) = s_cosamp; % current signal estimate
    if (norm(xcosamp(:,kk)-xcosamp(:,kk-1)) < 1e-2*norm(xcosamp(:,kk))) %break if error less than tolerance
       break;
    end
    
end

xcosamp(:,kk+1:end)=[];
xhat = xcosamp(:,end);



