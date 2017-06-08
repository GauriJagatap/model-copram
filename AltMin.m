function [x,err_hist,p,x_init] = AltMin(y_abs,A,iter,tol,tol2,z)
%%updated 5/31/2017

%% Initialize parameters
[m,n] = size(A);
%If ground truth is unknown
if nargin < 6
    z = zeros(n,1);
end
p = zeros(m,1); %phase vector
M = zeros(n);   %correlation matrix
error_hist(1,1) = 1; %error in measurement model
error_hist(1,2) = 1; %relative error b/w subsequent iterations
y_abs2 = y_abs.^2; %quadratic measurements

%% Initialize x
for i = 1:m
    M = M + y_abs2(i)*A(i,:)'*A(i,:); % (n x n)
end
M = M/m;

svd_opt = 'svd'; %accurate, but slow
svd_opt = 'power'; %approximate, faster

switch svd_opt
    case 'svd'
        [u,sigma,v] = svd(M);
        x_init = u(:,1); %normalized top singular vector of M (n x 1)
    case 'power'
        x_init = svd_power(M);
end

x = x_init;

fprintf('\n#iter\t|y-Ax|\t\t|x-z|\n')
%% start descent
for t=1:iter

    p = sign(A*x);         
    x = cgsolve(A'*A, A'*(p.*y_abs), 1e-3,20,0); %tol = 1e-3;max_iter = 20;verbose = 0
    
    err_hist(t+1,1) = norm(y_abs-abs(A*x))/norm(y_abs);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\n',t,err_hist(t+1,1),err_hist(t+1,2))

    if (err_hist(t+1,1) < tol) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2)
        break;
    end
end

end