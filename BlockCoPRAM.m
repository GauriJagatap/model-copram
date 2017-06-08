%% BlockCoPRAM
%"Fast, sample-efficient algorithms for structured phase retrieval"
% Gauri Jagatap and Chinmay Hegde, Iowa State University
% email: gauri@iastate.edu

% recover x given y_abs, A.
% y_abs = |Ax|.
% A is sensing matrix (m x n); Normal(0,1) distributed; typically m < n.
% x is s-sparse (n x 1); uniform block length b, effective block sparsity k=s/b.
% y_abs are magnitude only measurements (m x 1).
function [x,err_hist,p,x_init] =  BlockCoPRAM(y_abs,A,b,K,max_iter,tol1,tol2,z)
%%updated 5/31/2017

%% initialize parameters
[m,n] = size(A);
Nb = n/b;
%If ground truth is unknown
if nargin < 7
    z = zeros(n,1);
end
p = zeros(m,1); %phase vector
error_hist(1,1) = 1;
error_hist(1,2) = 1;
Marg = zeros(1,n); %marginals
BMarg = zeros(1,Nb); %block marginals
MShat = zeros(K*b); %truncated correlation matrix 
AShat = zeros(m,K*b); %truncated sensing matrix 
y_abs2 = y_abs.^2; %quadratic measurements
phi_sq = sum(y_abs2)/m;
phi = sqrt(phi_sq); %signal power

%% s-Truncated sensing vectors

%signal marginals
Marg = ((y_abs2)'*(A.^2))'/m; % n x 1

%block marginals
for i=1:Nb
    BMarg(i,1) = norm(Marg((i-1)*b+1:i*b));%Nb x 1
end

[BMg BMgs] = sort(BMarg,'descend');
S0 = BMgs(1:K); %pick top K block-marginals
Sbhat = sort(S0); %sorted block indices
Shat = zeros(1,K*b); %sorted signal indices
for i=1:K
   Shat(1,(i-1)*b+1:i*b) = (Sbhat(i)-1)*b+1:b*Sbhat(i); %1 x b
end
%supp(Shat) = 1; figure; plot(supp); %support indicator
AShat = A(:,Shat); %m x K*b %sensing sub-matrix

%% Truncated measurements
TAF = 'on'; %consider only truncated measurements m' < m ; gives marginally better performance 
TAF = 'off'; %consider all measurements = m ; aligns with code presented in paper

switch TAF
    case 'on'
        card_Marg = ceil(m/6);
        %large measurements - amplitude flow
        for i=1:m
            M_eval(i) = y_abs(i)/norm(AShat(i,:));
        end 
        [Mm MmS] = sort(M_eval,'descend');
        Io = MmS(1:card_Marg); %indices between 1 to m
    case 'off'
        card_Marg = m;
        Io = 1:card_Marg;
end

%% Initialize x
%compute top singular vector according to thresholded (and possibly truncated) sensing vectors 
for i = 1:card_Marg
    ii = Io(i);
    MShat = MShat + (y_abs2(ii))*AShat(ii,:)'*AShat(ii,:); % (K*b x K*b)
end

svd_opt = 'svd'; %more accurate, but slower for larger dimensions
svd_opt = 'power'; %approximate, faster for lager dimensions

switch svd_opt
    case 'svd'
        [u,sigma,v] = svd(MShat);
        v1 = u(:,1); %top singular vector of MShat, normalized - s x 1
    case 'power'
        v1 = svd_power(MShat);
end

v = zeros(n,1);
v(Shat,1) = v1;
x_init = phi*v; %ensures that the energy/norm of the initial estimate is close to actual
x = x_init;

%% start decent

fprintf('\n#iter\t|y-Ax|\t\t|x-z|\n')
for t=1:max_iter

    p = sign(A*x);
    x = jsmp(p.*y_abs/sqrt(m), A/sqrt(m),b, K,10,x); %Its = 10 
    err_hist(t+1,1) = norm(y_abs-abs(A*x))/norm(y_abs);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\n',t,err_hist(t+1,1),err_hist(t+1,2))

    if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2)
        break;
    end
    
end


end