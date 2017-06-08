%% plots algorithm running time
function [] = plot_runtime_block(string)

%%updated 5/31/2017

    cd('results')
    load(string)
    cd('..')
    err_thres = 0.05;
    
    kb = size(K_b,2);
    ttime = zeros(kb,length(mspan),2);
    
    for k=1:kb
        ttime(k,:,1) = mean(ttimer(k,:,:,1),3); %recover error - BlockCoPRAM
        ttime(k,:,2) = mean(ttimer(k,:,:,2),3); %recover error - CoPRAM
        
        figure;
        hold on;
        t = ['Mean running time : K = ',num2str(K_b(1,k)),', J = ',num2str(K_b(2,k)),', n = ',num2str(n)];
        xlab = 'no. of measurements (m)';
        ylab = 'running time in seconds';
        plot(mspan,ttime(k,:,1),'o-b')
        plot(mspan,ttime(k,:,2),'x-m')
        axisfortex(t,xlab,ylab);
        legend('BlockCoPRAM','CoPRAM')
        box on
        grid on
    end
    
    
end