%% plots algorithm running time
function [] = plot_runtime_trees(string)

%%updated 5/31/2017

    cd('../results')
    load(string)
    err_thres = 0.05;

    ttime = zeros(length(fspan),length(mspan),2);

    for f_iter = 1:length(fspan)
    
    ttime(f_iter,:,1) = mean(ttimer(f_iter,:,:,1),3); %recover error - TreeCoPRAM
    ttime(f_iter,:,2) = mean(ttimer(f_iter,:,:,2),3); %recover error - CoPRAM

    figure;
    hold on;
    f = fspan(f_iter);
    t = ['Mean running time : s = ',num2str(round(f*n)),', n = ',num2str(n)];
    xlab = 'no. of measurements (m)';
    ylab = 'running time in seconds';
    plot(mspan,ttime(f_iter,:,1),'o-b')
    plot(mspan,ttime(f_iter,:,2),'x-m')
    axisfortex(t,xlab,ylab);
    legend('TreeCoPRAM','CoPRAM')
    box on
    grid on
    
    end
    
    cd('../phase-retrieval')
    
end