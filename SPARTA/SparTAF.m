% code by Gauri Jagatap (gauri@iastate.edu)
% based on algorithm described in 
% G. Wang, G. B. Giannakis, J. Chen and M. Akçakaya,
% "SPARTA: Sparse phase retrieval via Truncated Amplitude flow", 
% ICASSP 2017.
function [x,err_hist,p,x_init] = SparTAF(y_abs,A,s,max_iter,tol1,tol2,z)
%updated 5/31/2017

%% initialize parameters
[m, n] = size(A);
%If ground truth is unknown
if nargin < 7
    z = zeros(n,1);
end
error_hist(1,1) = 1;
error_hist(1,2) = 1;
Marg = zeros(1,n); %marginals
MShat = zeros(s); %truncated correlation matrix
AShat = zeros(m,s); %truncated sensing matrix
y_abs2 = y_abs.^2;
phi_sq = sum(y_abs2)/m;
phi = sqrt(phi_sq); %signal power
%SPARTA parameters
mu = 1;
gamma = 0.7;

%% s-Truncated sensing vectors

%signal marginals
Marg = ((y_abs2)'*(A.^2))'/m; % n x 1
[Mg MgS] = sort(Marg,'descend');
S0 = MgS(1:s); %pick top s-marginals
Shat = sort(S0); %store indices in sorted order
%supp(Shat) = 1; figure; plot(supp); %support indicator
AShat = A(:,Shat); % m x s %sensing sub-matrix

%% Truncated measurements
card_Marg = ceil(m/6);
%large measurements - amplitude flow
for i=1:m
    M_eval(i) = y_abs(i)/norm(AShat(i,:));
end 
[Mm MmS] = sort(M_eval,'descend');
Io = MmS(1:card_Marg); %indices between 1 to m

%% Initialize x
%compute top singular vector according to thresholded sensing vectors and large measurements
for i = 1:card_Marg
    ii = Io(i);
    MShat = MShat + (y_abs2(ii))*AShat(ii,:)'*AShat(ii,:); % (s x s)
end

svd_opt = 'svd'; %more accurate, but slower for larger dimensions
svd_opt = 'power'; %approximate, faster

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

%% start descent 
fprintf('\n#iter\t|y-Ax|\t\t|x-z|\n')
for t = 1:max_iter
      
    It_act = A*x;
    It_mag = abs(It_act);
    It = find(It_mag > y_abs/(1+gamma));
    
    sum_TAF = zeros(n,1);
    for i=1:length(It)
        ii = It(i); 
        sum_TAF = sum_TAF + (It_act(ii,1) - y_abs(ii)*(It_act(ii,1))/It_mag(ii,1))*A(ii,:)'; %(n x 1)
    end
    
    grad_x = (mu/m)*sum_TAF;
    arg_TAF = x - grad_x;
    x = truncated_AF(arg_TAF,s);
    
    %store error history
    err_hist(t+1,1) = norm(y_abs-abs(A*x))/norm(y_abs);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\n',t,err_hist(t+1,1),err_hist(t+1,2))
    if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2) | (err_hist(t+1,2)>5)
        break; %tends to diverge sometimes, hence third condition
    end
    
end
p = sign(A*x);
end