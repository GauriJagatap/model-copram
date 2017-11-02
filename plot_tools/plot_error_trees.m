function [] = plot_error_trees(stringload)

%%updated 5/30/2017

cd('../results')
load(stringload);
errors = zeros(length(fspan),length(mspan),2);

for f_iter = 1:length(fspan)
%average (mean) over trials 
    errors(f_iter,:,1) = mean(err_sig(f_iter,:,:,1),3); %recover error - TreeCoPRAM
    errors(f_iter,:,2) = mean(err_sig(f_iter,:,:,2),3); %recover error - CoPRAM

    figure;
    hold on;
    f = fspan(f_iter);
    t = ['mean error in retrieval w.r.t ground truth: s = ',num2str(round(f*n)),', n = ',num2str(n)];
    xlab = 'no. of measurements (m)';
    ylab = 'mean error';
    plot(mspan,errors(f_iter,:,1),'o-b')
    plot(mspan,errors(f_iter,:,2),'x-m')
    axisfortex(t,xlab,ylab);
    legend('TreeCoPRAM','CoPRAM')
    box on
    grid on

end

cd('../phase-retrieval')

end