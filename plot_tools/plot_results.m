function [] = plot_results(z,x,y_ph,C)
%% plot recovered signal w.r.t original signal
%%updated 1/18/2017
figure;
hold on;
plot(z);
plot(x);
legend('original signal - z','recovered signal - x')

if y_ph == C
    fprintf('\nPhase recovered - Overall phase - correct\n');
elseif y_ph == -C
    fprintf('\nPhase recovered - Overall phase - flipped\n');
else
    fprintf('\nPhase not recovered\n');
end

end
