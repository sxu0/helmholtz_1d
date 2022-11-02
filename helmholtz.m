%% running the 1D Helmholtz solver!
clear; close all; clc;

k_arr = linspace(0.1, 8, 40);
n_train = 20;
mag_inc_wave = 1;

p = 2;
nelem = 1000;

% plot eigenvalues & associated errors
% examine + set ideal error tolerance below
[~, ~, ~] = snapshots(n_train, 1, p, nelem, 'log', mag_inc_wave);

%% OFFLINE
err_tol = 2e-6;
[A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = rb_offline( ...
    n_train, err_tol);

%% ONLINE
[u_N_arr, s_N_arr, u_arr] = rb_online( ...
    A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mag_inc_wave, k_arr);

%% VISUALIZATION
if ~exist('N', 'var')
    N = size(u_N_arr, 1);
end
anim_soln_real = visualize(k_arr, u_arr, mesh, ref, mag_inc_wave, n_train, N);
anim_soln_imag = visualize(k_arr, u_arr, mesh, ref, mag_inc_wave, n_train, N, 'imag');
