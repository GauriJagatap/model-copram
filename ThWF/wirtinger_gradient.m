function gradf_x = wirtinger_gradient(x,y,A)   
    [m,n] = size(A);
    yy = A*x; 
    coeff = (yy.^2-y).*(yy);
    gradf_x = A'*coeff/m;
end