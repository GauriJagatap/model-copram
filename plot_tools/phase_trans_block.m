%% plots phase transition plot
function [] = phase_trans_block(string)

%%updated 5/31/2017

    cd('results')
    load(string)
    cd('..')
    err_thres = 0.05;
    
    kb = size(K_b,2);
    prob_err = zeros(kb,length(mspan),2);
    
    for k=1:kb    
        err_block = err_sig(k,:,:,1);
        prob_err(k,:,1) = sum(err_block>err_thres,3)/trials_M;
        err_copram = err_sig(k,:,:,2);
        prob_err(k,:,2) = sum(err_copram>err_thres,3)/trials_M;

        figure;
        hold on;
        plot(mspan,1-prob_err(k,:,1),'o-b');
        plot(mspan,1-prob_err(k,:,2),'x-m');
        t = ['phase transition: K = ',num2str(K_b(1,k)),', J = ',num2str(K_b(2,k)),', N = ',num2str(n)];
        xlab='no. of measurements (m)';
        ylab='probability of success in recovery';
        axisfortex(t,xlab,ylab)
        legend('BlockCoPRAM','CoPRAM')
        grid on
        box on  
    end
    
    kspan = K_b(1,:).*K_b(2,:);
    b = mean(K_b(2,:));
    if prod((K_b(2,:)-b)==0)==1
        %phase transition diagram for variable 's' and 'm'; if block length is constant
        cd('plot_tools')
        prob_suc = 1-prob_err;
        t = [ 'Block CoPRAM : n = ',num2str(n),', b = ',num2str(b)];
        plot_diagram(t,mspan,kspan,prob_suc(:,:,1));
        t = [ 'CoPRAM : n = ',num2str(n)];
        plot_diagram(t,mspan,kspan,prob_suc(:,:,2));
    end
    cd('..')
    
end