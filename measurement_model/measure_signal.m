function [y_abs,y_ph,A] = measure_signal(m,z)
%edited 2/15/2017
n = length(z);
%% signal measurement
A = randn(m,n);
y = A*z; %measurements
y_abs = abs(y);
y_ph = sign(y); %actual phase
end