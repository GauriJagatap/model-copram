%%updated 5/31/2017
%% phase retrieval for block sparse signals

if ~exist('results','dir')
    mkdir('results')
end
if ~exist('signal_model','dir')
    mkdir('signal_model')
end
addpath('utils','measurement_model','ThWF','SPARTA','plot_tools','results','signal_model')
close all;
clc;
clear all;

%% signal model
n = 3000; %signal length
%study variation of phase transition of Block CoPRAM with block length
%K_b = [10 5 4 2 1;2 4 5 10 20]; %K = block sparsity; b = block length;
%study phase transition of Block CoPRAM w.r.t CoPRAM
K_b = [1 2 3 4 5 6 7 8 9 10;5 5 5 5 5 5 5 5 5 5];
kb = size(K_b,2);

%% measurement params
mspan = 200:200:2000; %no. of measurements
ml = length(mspan);

%% recovery validity
trials_M = 1; %set to 50 or 100
err_sig = zeros(kb,ml,trials_M,2);
ttimer = zeros(kb,ml,trials_M,2);

for kj_iter = 1:kb
    for tr = 1:trials_M
        for m_iter = 1:ml
            K = K_b(1,kj_iter);
            b = K_b(2,kj_iter);
            m = mspan(m_iter);
            fprintf('\nTrial no. :%d\nNo. of measurements m :%d\nBlock sparsity K :%d\nBlock length b :%d\n',tr,m,K,b);

            %% generate signal and measurements  
            [z,z_ind] =  generate_signal(n,K,b);
            [y_abs,y_ph,A] = measure_signal(m,z);

            %% initialize paramters 
            iter = 30; %maximum iterations for CoPRAM based algorithms
            tol1 = 1e-3; %error tolerance for measurements
            tol2 = 1e-5; %relative error tolerance between subsequent iterations

            %% use CoPRAM/Block CoPRAM - recover x1/x2
            fprintf('\nRunning Block CoPRAM . . .\n');
            tic;
            [x1,err_hist1,C1,x1_init] = BlockCoPRAM(y_abs,A,b,K,iter,tol1,tol2,z);   
            ttimer(kj_iter,m_iter,tr,1) = toc;
            
            fprintf('\nRunning CoPRAM . . .\n');            
            tic;
            [x2,err_hist2,C2,x2_init] = CoPRAM(y_abs,A,b*K,iter,tol1,tol2,z);   
            ttimer(kj_iter,m_iter,tr,2) = toc;            
            
            [err_sig(kj_iter,m_iter,tr,1) err_ind1] = approx_err(x1,z);
            [err_sig(kj_iter,m_iter,tr,2) err_ind2] = approx_err(x2,z);

            %global phase compensation
            x1 = x1*(-1)^(err_ind1-1);
            x2 = x2*(-1)^(err_ind2-1);
            
            %% Block CoPRAM - results
            fprintf('\nError in approximation using Block CoPRAM is %2.4f\n',err_sig(kj_iter,m_iter,tr,1));
            fprintf('\nError in approximation using CoPRAM is %2.4f\n',err_sig(kj_iter,m_iter,tr,2));
            
            %plot if algorithm is run for 1 instance of s and m
            if (trials_M*ml*kb) == 1
                plot_results(z,x1,y_ph,C1)
                title('Signal recovery using Block CoPRAM')
                plot_results(z,x2,y_ph,C2)
                title('Signal recovery using CoPRAM')
            end
        end
    end
end

%% save results
cd('results')
    strr = ['BlockSparsePR_errors_','n',num2str(n),'_trials',num2str(trials_M),'.mat'];
    save(strr,'err_sig','n','K_b','mspan','trials_M','ttimer');
cd('..')

%% display results - phase transition plot
phase_trans_block(strr)

%% plot running time statistics 
% plot_runtime_block(strr)
