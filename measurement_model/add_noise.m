function y_noisy = add_noise(z,y_clean,NSR)

  sigma2 = NSR*norm(z)^2; 
  eps = sqrt(sigma2)*randn(size(y_clean)); 
  y_noisy = y_clean+eps;

end
