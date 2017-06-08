function thres = threshold_lvl(x,y,A)
    [m,n] = size(A);
    beta = 0.3; %adjust 0.01
    yy = (A*x).^2; 
    summ = sum(((yy- y).^2).*(yy));
    tau_sq = (beta*log(m*n)/m^2)*(summ);
    thres = sqrt(tau_sq);
end