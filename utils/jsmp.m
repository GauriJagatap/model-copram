% jsmp.m
% Performs CS reconstruction of 1D block-sparse signals
%
% INPUTS
% yy      : measurements (M x 1)
% Phi     : measurement matrix (M x N)
% J       : Block size
% K       : number of active blocks
% Its     : number of iterations
% 
% OUTPUTS
% xhat    : Signal estimate (N x 1) 
% xjsmp   : Matrix with N rows and at most Its columns; 
%           columns represent intermediate signal estimates   
% 
%
% CITE    : Richard Baraniuk, Volkan Cevher, Marco Duarte, Chinmay Hegde
%          "Model-based compressive sensing", submitted to IEEE IT, 2008.
% Created : Aug 2009.
% email   : chinmay@rice.edu
%
% The model for block-sparsity used here is equivalent to the JSM-2 
% model introduced in:
% Baron et al, "Distributed Compressive Sensing", preprint 2009

function [xhat,xjsmp] = jsmp(yy, Phi, J, K, Its, s_init)

yy = yy(:); % 
[M,N] = size(Phi);
num_blocks = round(N/J);

xjsmp = zeros(N,Its); 
kk=1; 
maxiter= 100;
verbose= 0;
tol= 1e-3;

if nargin < 6
    s_jsmp = zeros(N,1);
else 
    s_jsmp = s_init;
end

while le(kk,Its),
    rjsmp = yy - Phi*(s_jsmp);
    proxy_jsmp = Phi'*(rjsmp);
   
    %--------------------------------------------%
    %-----Insert approximation algorithm here----%
    %--------------------------------------------%
    %-----JSM-MP
    % pick the K largest blocks. Easy.
    proxy_jsmp_block = reshape(proxy_jsmp, J, num_blocks);
    [trash,blockww] = sort(sum(proxy_jsmp_block.^2,1),'descend');
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:K)) = 1;
    newsupp = reshape(newsupp, N, 1);
    tt_jsmp=union(find(ne(s_jsmp,0)),find(newsupp==1));  
    %-----/Insert----%

    %------Estimate------%
    [w_jsmp, res, iter] = cgsolve(Phi(:,tt_jsmp)'*Phi(:,tt_jsmp), Phi(:,tt_jsmp)'*yy,...
                                        tol,maxiter, verbose);
    bb1= zeros(N,1);
    bb1(tt_jsmp)= w_jsmp;
    %---Prune----%
    kk = kk+1;   
    
    bb1_block = reshape(bb1, J, num_blocks);
    [trash,blockww] = sort(sum(bb1_block.^2,1),'descend');
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:K)) = 1;
    newsupp = reshape(newsupp, N, 1);
    s_jsmp=0*s_jsmp;
    s_jsmp = bb1.*newsupp;
    xjsmp(:,kk) = s_jsmp;
    
    if (norm(xjsmp(:,kk)-xjsmp(:,kk-1)) < 5e-3*norm(xjsmp(:,kk)))
       break;
    end
end
xjsmp(:,kk+1:end)=[];
xhat = xjsmp(:,end);