function [] = plot_diagram(t,mspan,kspan,prob_suc)

    figure; 
    imagesc(mspan,kspan,prob_suc); 
    colormap(gray); 
    set(gca,'YDir','normal'); 
    set(gca,'FontSize',20); 
    axis square;
    xlabel('No. of samples (m)')
    ylabel('Sparsity (s)')
    title(t);
    
end