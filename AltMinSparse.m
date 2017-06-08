function [x,err_history,p,x_init] = AltMinSparse(y_abs,A,s,iter,tol1,tol2,z)
%%edited 5/31/2017

%% Initialize parameters
[m,n] = size(A);
%If ground truth is unknown
if nargin < 7
    z = zeros(n,1);
end
p = zeros(m,1); %phase vector
Marg = zeros(n,1); %signal marginals 
x = zeros(n,1); 

%% Initialize x
Marg = abs(A')*abs(y_abs);
[S1 S2] = sort(Marg,'descend');
S = S2(1:s); %s
As = A(:,S); %m x n --> m x s
ys = y_abs; %m x 1
zs = z(S,1); %s x 1

if zs == 0 %if support picked up contains none of the true support
    fprintf('\nDid not pick up anything from true support in the initial stage...\n')
    err_history = 0;
    x_init = x;
else    
    [xs,err_history,p,x_init] = AltMin(ys,As,iter,tol1,tol2,zs);
    x(S,1) = xs; %s x 1
end

end