%% plots algorithm running time
function [] = plot_runtime(string)

%%updated 5/31/2017

    cd('results')
    load(string)
    cd('..')
    err_thres = 0.05;
    ttime = zeros(length(kspan),length(mspan),4);
    
    for k=1:length(kspan)
        %average (mean) over trials for fixed M
        ttime(k,:,1) = mean(ttimer(k,:,:,1),3); %recover error - CoPRAM
        ttime(k,:,2) = mean(ttimer(k,:,:,2),3); %recover error - AltMinSparse
        ttime(k,:,3) = mean(ttimer(k,:,:,3),3); %recover error - ThWF
        ttime(k,:,4) = mean(ttimer(k,:,:,4),3); %recover error - SparTA

        figure;
        hold on;
        t = ['Mean running time : s = ',num2str(kspan(k)),', n = ',num2str(n)];
        xlab = 'no. of measurements (m)';
        ylab = 'running time in seconds';
        plot(mspan,ttime(k,:,1),'o-b')
        plot(mspan,ttime(k,:,2),'+-g')
        plot(mspan,ttime(k,:,3),'d-r')
        plot(mspan,ttime(k,:,4),'s-m')
        axisfortex(t,xlab,ylab);
        legend('CoPRAM','AltMinSparse','ThWF','SparTA')
        box on
        grid on
    end
    
end