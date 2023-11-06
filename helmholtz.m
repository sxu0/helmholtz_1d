% running the 1D Helmholtz solver!
clear; close all; clc;
addpath('lib')

%% CONFIG
save_figs = false;

kk_train = linspace(0, 2, 41) .* pi;
n_train = numel(kk_train);
wave_ampl = 1;

p = 2;
nelem = 1000;

%% RB ACCURACY
% plot eigenvalues & associated errors
% examine + set ideal error tolerance below
[~, ~, ~] = snapshots(kk_train, 1, p, nelem, 'log', wave_ampl);
if save_figs
    print(fullfile('figs', 'sorted_eigenvalues_and_errors_real.png'), '-dpng')
end

%% OFFLINE STEPS
err_tol = 5.5e-6;
[A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = rb_offline( ...
    kk_train, err_tol);

%% ONLINE STEPS
kk_test = sort(random('Uniform', 0, 2*pi, [1, 81]));
[uu_N, qq_N, uu] = rb_online( ...
    A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, wave_ampl, kk_test);

%% VISUALIZATION
if ~exist('N', 'var')
    N = size(uu_N, 1);
end
if save_figs
    output_gif = fullfile('figs', ...
        ['soln_ntrain=' num2str(n_train) '_N=' num2str(N) '.gif']);
else
    output_gif = '';
end
visualize( ...
    kk_test, uu, mesh, ref, wave_ampl, n_train, N, output_gif);

%% QUANTITY OF INTEREST
qq = quantity_of_interest(uu, wave_ampl);
figure;
plot(kk_test, qq, 'o')
xlim('padded')
ylim('padded')
xlabel('$k$')
ylabel('$q_N(k)$')
title(['RB output $q_N$ by $k$, $n_\mathrm{train}=' ...
    num2str(n_train) '$, $n_\mathrm{solve}=' num2str(numel(kk_test)) '$'])

%% FE ERROR CONVERGENCE
nnelem = 2 .^ (1:12);
[u_err_fe, q_err_fe] = error_fe(nnelem, 2*pi, wave_ampl, p, true);
% see error_fe.m for a couple of things to check

%% RB ERROR CONVERGENCE
NN = 2:size(uu_N, 1);
[u_err_NN, q_err_NN] = error_rb( ...
    NN, kk_test, kk_train, err_tol, wave_ampl, p, nelem, true);
