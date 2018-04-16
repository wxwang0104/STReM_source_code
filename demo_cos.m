clear all
close all
clc

load('mask_record_ww.mat');%phase mask of different orientations

%% Generation of A matrix and the partial difference of A
[A, dx, dy, dz] = gene_PSF_matrix(mask_record);
%addpath('C:\Users\ww20\Desktop\06262015');
%% Recovery 
[im_raw, rec_pos] = run_recovery(A, dx, dy, dz, mask_record);
%% Plotting
run_plot_result;



