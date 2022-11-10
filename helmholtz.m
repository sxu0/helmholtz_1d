% running the 1D Helmholtz solver!
clear; close all; clc;

%% CONFIG
save_figs = false;

k_arr = linspace(0.1, 8, 40);
n_train = 20;
mag_inc_wave = 1;

p = 2;
nelem = 1000;

%% RB ACCURACY
% plot eigenvalues & associated errors
% examine + set ideal error tolerance below
[~, ~, ~] = snapshots(n_train, 1, p, nelem, 'log', mag_inc_wave);
if save_figs
    print(fullfile('figs', 'sorted_eigenvalues_and_errors_real.png'), '-dpng')
end

%% OFFLINE STEPS
err_tol = 2e-6;
[A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = rb_offline( ...
    n_train, err_tol);

%% ONLINE STEPS
[u_N_arr, s_N_arr, u_arr] = rb_online( ...
    A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mag_inc_wave, k_arr);

%% VISUALIZATION
if ~exist('N', 'var')
    N = size(u_N_arr, 1);
end
if save_figs
    output_gif = fullfile('figs', ...
        ['soln_ntrain=' num2str(n_train) '_N=' num2str(N) '.gif']);
else
    output_gif = '';
end
visualize( ...
    k_arr, u_arr, mesh, ref, mag_inc_wave, n_train, N, output_gif);
