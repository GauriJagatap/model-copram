function x_out = truncated_AF(x_in,s)
    %updated 4/18/2017
    [~,ii] = sort(abs(x_in),'descend');
    x_out = zeros(size(x_in));
    x_out(ii(1:s),1) = x_in(ii(1:s),1);
end