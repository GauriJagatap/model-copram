%% plots phase transition plot
function [] = phase_trans_trees(string)

%%updated 5/30/2017

    cd('../results')
    load(string)
    err_thres = 0.05;    
    prob_err = zeros(length(fspan),length(mspan),2);
       
    for f_iter = 1:length(fspan)
        
    err_trees = err_sig(f_iter,:,:,1);
    prob_err(f_iter,:,1) = sum(err_trees>err_thres,3)/trials_M;
    err_copram = err_sig(f_iter,:,:,2);
    prob_err(f_iter,:,2) = sum(err_copram>err_thres,3)/trials_M;
    
    figure;
    hold on;
    plot(mspan,1-prob_err(f_iter,:,1),'o-b');
    plot(mspan,1-prob_err(f_iter,:,2),'x-m');
    f = fspan(f_iter);
    t = ['phase transition: s = ',num2str(round(f*n)),', n = ',num2str(n)];
    xlab='no. of measurements (m)';
    ylab='probability of success in recovery';
    axisfortex(t,xlab,ylab)
    legend('TreeCoPRAM','CoPRAM')
    grid on
    box on  
    
    end
    cd('../phase-retrieval')
    
end