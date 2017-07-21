% code by Gauri Jagatap (gauri@iastate.edu)
% based on algorithm described in
% Cai, T. Tony, Xiaodong Li, and Zongming Ma. 
% "Optimal rates of convergence for noisy sparse phase retrieval 
% via thresholded Wirtinger flow." 
% The Annals of Statistics 44.5 (2016): 2221-2251.
function [x,err_hist,p,x_init] = Thresholded_WF(y_twf,A,K,max_iter,tol1,tol2,z)
%updated 5/31/2017

%% initialize parameters
[m,n] = size(A);
%If ground truth is unknown
if nargin < 7 
    z = zeros(n,1);
end
error_hist(1,1) = 1;
error_hist(1,2) = 1;
Marg = zeros(1,n); %marginals
phi_sq = sum(y_twf)/m;
phi = sqrt(phi_sq); %signal power
%Thresholded WF parameters
alpha = 1.5;
mu = 0.23; 
thres_param = (1+alpha*sqrt(log(m*n)/m))*phi_sq ;

%% Thresholded sensing vectors

%signal marginals
Marg = (y_twf'*(A.^2))'/m; % n x 1
S0 = find(Marg>thres_param);
ss = length(S0);
Shat = sort(S0);
MShat = zeros(ss);
AShat = zeros(m,ss);
%supp(Shat) = 1; figure; plot(supp); %support indicator
AShat = A(:,Shat); % m x ss

%compute top singular vector according to thresholded sensing vectors
for i = 1:m
    MShat = MShat + (y_twf(i))*AShat(i,:)'*AShat(i,:); % (ss x ss)
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

%% gradient descent
fprintf('\n#iter\t|y-Ax|\t\t|x-z|\n')

for t = 1:max_iter
    %threshold level for gradient descent  
    thres = threshold_lvl(x,y_twf,A);
    
    %gradient
    gradf_x = wirtinger_gradient(x,y_twf,A);
    
    %update
    xx = x-(mu/phi_sq)*gradf_x;   
    x = soft_threshold(xx,mu*thres/phi_sq);
   
    %store error history
    err_hist(t+1,1) = norm(y_twf-(A*x).^2)/norm(y_twf);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\n',t,err_hist(t+1,1),err_hist(t+1,2))
    if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2) | (err_hist(t+1,2)>5)
        break; %tends to diverge sometimes, hence third condition
    end
end
p = sign(A*x);
end