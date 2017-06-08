function [err, err_ind] = approx_err(x,z)
%edited 2/15/2017
    err1 = norm(x-z)/norm(z);
    err2 = norm(x+z)/norm(z); %if phase is flipped
    [err, err_ind] = min([err1,err2]);    
end