function [z,z_ind] = generate_signal(n,K,b)
%edited 2/15/2017
N_bl = floor(n/b);
save_str = ['block_sparse_sig_',num2str(b),'_',num2str(K),'_',num2str(n),'.mat'];
cd('signal_model')
    if exist(save_str) == 0 %if signal with params b,K,n does not exist, create and save it
        bl_z_ind = randperm(N_bl,K); %creates sparse block vector
        z_ind = [];
        for i=1:length(bl_z_ind) %selects out all indices corresponding to selected blocks
            z_ind = [z_ind (b*(bl_z_ind(i)-1)+1):((b*(bl_z_ind(i)-1)+1)+b-1)];
        end
        z = zeros(n,1);
        z(z_ind) = randn(K*b,1); %generate K*b sparse signal (n x 1)
        z = z/norm(z);   %need not be normalized     
        save(save_str,'z','z_ind','b','K','n');
    end
    
    %load saved signal with params b,K,n
    load(save_str);
    cd('..')
end