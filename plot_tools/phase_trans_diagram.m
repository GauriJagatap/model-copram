%phase trans diagram

function [] = phase_trans_diagram(string)

%%updated 5/30/2017

    cd('results')
    load(string)
    err_thres = 0.05;    
    prob_err = zeros(length(fspan),length(mspan),2);
       
    for f_iter = 1:length(fspan)
        
    err_trees = err_sig(f_iter,:,:,1);
    prob_err(f_iter,:,1) = sum(err_trees>err_thres,3)/trials_M;
    err_copram = err_sig(f_iter,:,:,2);
    prob_err(f_iter,:,2) = sum(err_copram>err_thres,3)/trials_M;
    
    end

    prob_suc_treecopram = 1-prob_err(:,:,1);
    prob_suc_copram = 1-prob_err(:,:,2);
   
    figure; 
    sspan = round(fspan*n);
    mm = (0.05:0.05:0.5);%(floor(100*mspan/n)/100)/4;
    imagesc(mm,fspan/4,prob_suc_treecopram); 
    colormap(gray); 
    set(gca,'YDir','normal'); 
    set(gca,'FontSize',20); 
    ax = gca;
    ax.XTick = mm(2:2:10);
    ax.YTick = fspan(2:2:10)/4;
    axis square;
    xlabel('No. of samples/n = m/n')
    ylabel('Sparsity/n = s/n')
    title('Tree CoPRAM')
   
    figure; 
    imagesc(mm,fspan/4,prob_suc_copram); 
    colormap(gray); 
    set(gca,'YDir','normal'); 
    set(gca,'FontSize',20); 
    ax = gca;
    ax.XTick = mm(2:2:10);
    ax.YTick = fspan(2:2:10)/4;
    axis square;
    xlabel('No. of samples/n = m/n')
    ylabel('Sparsity/n =  s/n')
    title('CoPRAM')
    cd('..')
    
end