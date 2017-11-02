%%updated 5/31/2017
%% phase retrieval for sparse signals using AltMin and CoSaMP (Compressive AltMinPhase)

path(path,'utils');
%path(path,'Signl');
path(path,'measurement_model');
path(path,'plot_tools');
path(path,'test_images');
path(path,'treeapprox');

close all;
clc;
clear all;

%% signal model
N = 512^2; %signal length
scal = 1/8;
n = N*(scal^2);

%load image
image = 'lovett';
if strcmp(image, 'lovett')
    im = imread('test_images/Lovett_Hall.jpg');
    im = im2double(im);
    im = im(:,(400-250):(400+250),:);
    nn = 512;
    im_sq = imresize(im,[nn nn],'bicubic');
    imgray = rgb2gray(im_sq);
elseif strcmp(image, 'stata')
    im = imread('test_images/stata.jpg');
    im = im2double(im);
    im = im(:,(1024 - 768 + 1):1024,:);
    nn = 512;
    im_sq = imresize(im,[nn nn],'bicubic');
    imgray = rgb2gray(im_sq);
else
    error('Unknown image')
end           
imgray = imresize(imgray,scal);
nn = scal*nn;
lev = floor(log(nn)/log(2));

%important: 'per' ensures you get a non-redundant transform.
wavmode = 'per';
dwtmode(wavmode,'nodisp')
[c,c_ind] = wavedec2(imgray,lev,'db8');
c_rearranged = rearrange_wavedec2(c, 'forward');

%% sparsity
fspan = 0.01:0.01:0.1;
fl = length(fspan);

%% measurement params
mspan = ceil((0.05:0.05:0.5)*n/4); %no. of measurements
ml = length(mspan);

%% recovery validity
trials_M = 1;
err_sig = zeros(fl,ml,trials_M,2);
ttimer = zeros(fl,ml,trials_M,2);

for tr = 1:trials_M
    for f_iter = 1:fl
    for m_iter = 1:ml
        m = mspan(m_iter);
        f = fspan(f_iter);

        %% load signal, set sparsity and measurements  

        s = round(f*n/4);
        fprintf('\nTrial no. :%d\nNo. of measurements M :%d Sparsity s :%d\n',tr,m,s);
        order = 4; % tree order; 4 for 2D images
        fprintf('Computing thresholded version of the image ...');
        supp = treeexact_smalltable_wvtree(c_rearranged.^2,order,s);
        chat = c_rearranged .* supp;
        imhat = waverec2(rearrange_wavedec2(chat, 'backward'),c_ind,'db8');
%         figure, imshow([imgray imhat]), axis image;
%         title('Actual v/s Thresholded');
        z = chat;
        z = z(:);
        z_ind = supp;

        [y_abs,y_ph,A] = measure_signal(m,z);

        %% initialize paramters of Tree CoPRAM and CoPRAM
        iter = 35; %maximum iterations for AltMin based algorithms
        tol1 = 1e-3; %error tolerance for measurements
        tol2 = 1e-5; %relative error tolerance between subsequent iterations

        %% use Tree CoPRAM/ CoPRAM - recover x1/x2
        fprintf('\nRunning Tree CoPRAM . . .\n');
        tic;
        [x1,err_hist1,C1,x1_init] = TreeCoPRAM(y_abs,A,s,iter,tol1,tol2,z);   
        ttimer(m_iter,tr,1) = toc;

        fprintf('\nRunning CoPRAM . . .\n');            
        tic;
        [x2,err_hist2,C2,x2_init] = CoPRAM(y_abs,A,s,iter,tol1,tol2,z);   
        ttimer(m_iter,tr,2) = toc;            

        [err_sig(f_iter,m_iter,tr,1) err_ind1] = approx_err(x1,z);
        [err_sig(f_iter,m_iter,tr,2) err_ind2] = approx_err(x2,z);

        %global phase compensation
        x1 = x1*(-1)^(err_ind1-1);
        x2 = x2*(-1)^(err_ind2-1);

        %% AltMinCoSaMP - results
        fprintf('\nError in approximation using Tree CoPRAM is %2.4f\n',err_sig(f_iter,m_iter,tr,1));
        fprintf('\nError in approximation using CoPRAM is %2.4f\n',err_sig(f_iter,m_iter,tr,2));

        %plot if algorithm is run for few instances of K and M
        if (trials_M*ml*fl) == 1
            im_tree_exact = waverec2(rearrange_wavedec2(x1, 'backward'),c_ind,'db8');
            im_approx = waverec2(rearrange_wavedec2(x2, 'backward'),c_ind,'db8');
            figure, clf
            imshow([im_tree_exact im_approx]);
            title('BlockCoPRAM v/s CoPRAM');
            axis image            
        end
    end
    end
end

%% save results
cd('../results')
    strr = ['TreeSparsePR_errors_','n',num2str(n),'_s',num2str(s),'_trials',num2str(trials_M),'.mat'];
    save(strr,'err_sig','n','fspan','mspan','trials_M','ttimer');
cd('../phase-retrieval')

%% display results - mean error in recovery
%plot_error_trees(strr);

%% display results - phase transition plot
phase_trans_trees(strr)
phase_trans_diagram(strr)

%% plot running time statistics 
%plot_runtime_trees(strr)
