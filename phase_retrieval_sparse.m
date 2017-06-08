%% phase retrieval for sparse signals 
%stores signal models (.mat) and results (.mat/fig/jpg). not committed to git repo.
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
n = 1000; %signal length
kspan = 3:3:30; %sparsity vector
kl = length(kspan);
b = 1; %block length (non-essential for evaluating standard sparse models; trivially b=1)

%% measurement params
mspan = 100:100:1000; %no. of measurements
ml = length(mspan);

%% recovery validity
trials_M = 1; %set to 50 or 100
err_sig = zeros(kl,ml,trials_M,4);
ttimer = zeros(kl,ml,trials_M,4);

for k_iter = 1:kl
    for tr = 1:trials_M
        for m_iter = 1:ml
            s = kspan(k_iter)*b;
            m = mspan(m_iter);
            fprintf('\nTrial no. :%d\nNo. of measurements M :%d\nSparsity K :%d\n',tr,m,s);

            %% generate signal and measurements  
            [z,z_ind] =  generate_signal(n,s,b);
            [y_abs,y_ph,A] = measure_signal(m,z);
            
            %add noise to measurements in CoPRAM / SPARTA
            NSR = 0; %noise-signal ratio          
            y = add_noise(z,y_abs.*y_ph,NSR);            
            y_abs = abs(y); y_ph = sign(y);
            
            %measurements required for Thresholded Wirtinger Flow + Noise
            y_twf = y_abs.^2;  
            y_twf = add_noise(z,y_twf,NSR);

            %% initialize paramters 
            iter = 30; %maximum iterations for AltMin based algorithms
            iter_gd = 30; %maximum iterations for WF/AF based algorithms (typically take higher no. to convg)
            tol1 = 1e-3; %error tolerance for measurements
            tol2 = 1e-5; %error tolerance between subsequent iterations

            %% use CoPRAM/AltMinSparse/ThWF/SPARTA - recover x1/x2/x3/x4 
            fprintf('\nRunning CoPRAM . . .\n');
            tic;
            [x1,err_hist1,C1,x1_init] = CoPRAM(y_abs,A,s,iter,tol1,tol2,z);   
            ttimer(k_iter,m_iter,tr,1) = toc;
            
            fprintf('\nRunning AltMinSparse . . .\n');
            tic;
            [x2,err_hist2,C2,x2_init] = AltMinSparse(y_abs,A,s,iter,tol1,tol2,z);
            ttimer(k_iter,m_iter,tr,2) = toc;
            
            fprintf('\nRunning Thresholded Wirtinger Flow . . .\n');      
            tic;
            [x3,err_hist3,C3,x3_init] = Thresholded_WF(y_twf,A,s,iter_gd,tol1,tol2,z);
            ttimer(k_iter,m_iter,tr,3) = toc;
            
            fprintf('\nRunning Sparse Truncated Amplitude Flow . . .\n');
            tic;
            [x4,err_hist4,C4,x4_init] = SparTAF(y_abs,A,s,iter_gd,tol1,tol2,z);
            ttimer(k_iter,m_iter,tr,4) = toc;
            
            %error w.r.t ground truth
            [err_sig(k_iter,m_iter,tr,1) err_ind1] = approx_err(x1,z);
            [err_sig(k_iter,m_iter,tr,2) err_ind2] = approx_err(x2,z);
            [err_sig(k_iter,m_iter,tr,3) err_ind3] = approx_err(x3,z);
            [err_sig(k_iter,m_iter,tr,4) err_ind4] = approx_err(x4,z);

            %global phase compensation
            x1 = x1*(-1)^(err_ind1-1);
            x2 = x2*(-1)^(err_ind2-1);
            x2 = x2*(-1)^(err_ind3-1);
            x3 = x3*(-1)^(err_ind4-1);

            %% print recovery error w.r.t ground truth
            fprintf('\nError in approximation using CoPRAM is %2.4f\n',err_sig(k_iter,m_iter,tr,1));
            fprintf('\nError in approximation using AltMinSparse is %2.4f\n',err_sig(k_iter,m_iter,tr,2));
            fprintf('\nError in approximation using Thresholded Wirtinger Flow is %2.4f\n',err_sig(k_iter,m_iter,tr,3));
            fprintf('\nError in approximation using Sparse Truncated Amplitude Flow is %2.4f\n',err_sig(k_iter,m_iter,tr,4));

            %plot reconstructed signal if algorithm is run for 1 instance of fixed s and m
            if (trials_M*ml*kl) == 1
                plot_results(z,x1,y_ph,C1)
                title('Signal recovery using CoPRAM')
                plot_results(z,x2,y_ph,C2)
                title('Signal recovery using AltMinSparse')
                plot_results(z,x3,y_ph,C3)
                title('Signal recovery using Thresholded Wirtinger Flow')
                plot_results(z,x4,y_ph,C4)
                title('Signal recovery using Sparse Truncated Amplitude Flow')
            end
        end
    end
end

%% save results
cd('results')
    strr = ['SparsePR_errors_','n',num2str(n),'_trials',num2str(trials_M),'.mat'];
    save(strr,'err_sig','n','kspan','mspan','trials_M','ttimer');
cd('..')

%% display results - phase transition plot
phase_trans(strr)

%% plot running time statistics 
% plot_runtime(strr)
