%% CoPRAM

function [x,err_hist,p,x_init] =  CoPRAM(y_abs,A,s,max_iter,tol1,tol2,z)
%%updated 5/31/2017

%% initialize parameters
[m,n] = size(A);
%If ground truth is unknown
if nargin < 7
    z = zeros(n,1);
end
p = zeros(m,1); %phase vector
error_hist(1,1) = 1; %stores error in measurement model
error_hist(1,2) = 1; %stores relative error b/w iterations
Marg = zeros(1,n); %marginals
MShat = zeros(s); %truncated correlation matrix
AShat = zeros(m,s); %truncated sensing matrix
supp = zeros(1,n); %indicator for initial support Shat
y_abs2 = y_abs.^2; %quadratic measurements
phi_sq = sum(y_abs2)/m;
phi = sqrt(phi_sq); %signal power

%% s-Truncated sensing vectors

%signal marginals
Marg = ((y_abs2)'*(A.^2))'/m; % n x 1
[Mg MgS] = sort(Marg,'descend');
S0 = MgS(1:s); %pick top s-marginals
Shat = sort(S0); %store indices in sorted order
%supp(Shat) = 1; figure; plot(supp); %support indicator
AShat = A(:,Shat); % m x s %sensing sub-matrix

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
        [~,MmS] = sort(M_eval,'descend');
        Io = MmS(1:card_Marg); %indices between 1 to m
    case 'off'
        card_Marg = m;
        Io = 1:card_Marg;
end

%% initialize x
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
for t=1:max_iter 
    
    p = sign(A*x);
    x = cosamp(p.*y_abs/sqrt(m), A/sqrt(m), s,10,x); %Its = 10
    
    err_hist(t+1,1) = norm(y_abs-abs(A*x))/norm(y_abs);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\n',t,err_hist(t+1,1),err_hist(t+1,2))
    if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2)
        break;
    end

end


end