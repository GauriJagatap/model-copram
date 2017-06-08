%% plots phase transition plot
function [] = phase_trans(string)

%%updated 5/31/2017

    cd('results')
    load(string)
    cd('..')
    err_thres = 0.05; %can vary 
    prob_err = zeros(length(kspan),length(mspan),4);
    

    for k_iter=1:length(kspan)    %displays phase transitions for a fixed 's'
        err_copram = err_sig(k_iter,:,:,1);
        prob_err(k_iter,:,1) = sum(err_copram>err_thres,3)/trials_M;
        err_AMSparse = err_sig(k_iter,:,:,2);
        prob_err(k_iter,:,2) = sum(err_AMSparse>err_thres,3)/trials_M;
        err_ThWF = err_sig(k_iter,:,:,3);
        prob_err(k_iter,:,3) = sum(err_ThWF>err_thres,3)/trials_M;
        err_SparTA = err_sig(k_iter,:,:,4);
        prob_err(k_iter,:,4) = sum(err_SparTA>err_thres,3)/trials_M;

        figure;
        hold on;
        plot(mspan,1-prob_err(k_iter,:,1),'o-b');
        plot(mspan,1-prob_err(k_iter,:,2),'s-g');
        plot(mspan,1-prob_err(k_iter,:,3),'d-r');
        plot(mspan,1-prob_err(k_iter,:,4),'x-m');
        t = ['phase transition: s = ',num2str(kspan(k_iter)),', n = ',num2str(n)];
        xlab='no. of measurements (m)';
        ylab='probability of success in recovery';
        axisfortex(t,xlab,ylab)
        legend('CoPRAM','AltMinSparse','ThWF','SPARTA')
        grid on
        box on  
    end
    
    %phase transition diagram for variable 's' and 'm'; comment if studying
    %variation of phase transition with block length.
    cd('plot_tools')
    prob_suc = 1-prob_err;
    t = [ 'CoPRAM : n = ',num2str(n)];
    plot_diagram(t,mspan,kspan,prob_suc(:,:,1));
    t = [ 'ThWF : n = ',num2str(n)];
    plot_diagram(t,mspan,kspan,prob_suc(:,:,3));
    t = [ 'SPARTA : n = ',num2str(n)];
    plot_diagram(t,mspan,kspan,prob_suc(:,:,4));
    cd('..')
    
end