%% testing helmholtz_FE_solver() function
clear; close all;

mag_inc_wave = 1;
p = 2;
nelem = 1000;
plt = true;

figure
for k = (11:60)/10
[U, ~, ~] = helmholtz_FE_solver(k, mag_inc_wave, p, nelem, plt);
end
hold off

% confirmed functional + verified output against analytical soln


%% testing helmholtz_assemble_inprod_matrix() function
clear; close all;

p = 2;
nelem = 1000;

X = helmholtz_assemble_inprod_matrix(p, nelem);

x_diag = diag(X);
x_supdiag = diag(X, 1);
x_subdiag = diag(X, -1);
x_blank_diag = diag(X, 2);

% confirmed functional, yet to verify


%% testing helmholtz_off_on_setup() function
clear; close all;

mag_inc_wave = 1;
p = 2;
nelem = 1000;

[A1, A2, A3, F] = helmholtz_off_on_setup(mag_inc_wave, p, nelem);

A1_diag = diag(A1);
A1_supdiag = diag(A1, 1);
A1_subdiag = diag(A1, -1);
A1_blank_diag = diag(A1, 2);
A2_diag = diag(A2);
A2_supdiag = diag(A2, 1);
A2_subdiag = diag(A2, -1);
A2_blank_diag = diag(A2, 2);
A3_diag = diag(A3);
A3_blank_diag = diag(A3, 1);

% confirmed functional, yet to verify


%% testing the heart of the solver!
% (next 3 code sections)
% helmholtz_snapshots(), helmholtz_RB_offline(), & helmholtz_RB_online()
clear; close all

k_arr = linspace(0.1, 8, 40);
n_train = 20;
mag_inc_wave = 1;

p = 2;
nelem = 1000;

% plot eigenvalues & associated errors to find ideal error tolerance
[~, ~, ~] = helmholtz_snapshots(n_train, 1, p, nelem, true, mag_inc_wave);

%% OFFLINE
err_tol = 2.2e-7;
[A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = helmholtz_RB_offline( ...
    n_train, err_tol);

%% ONLINE & VISUALIZATION
[u_N_arr, s_N_arr, u_arr] = helmholtz_RB_online( ...
    A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mag_inc_wave, k_arr);
fig_soln = figure;
labels = cell(size(k_arr));
for j = 1:numel(k_arr)
    U = u_arr(:, j);
    % next 7 lines are sourced from `femmat/util/plot/plotfield1d()`
    nlp = size(ref.shp, 2);
    xxplot = linspace(0, 1, 10*mesh.p)';
    shp = shape1d(ref.p, ref.xint, xxplot);
    xtri = reshape(mesh.coord(mesh.tri), [nelem, nlp]);
    utri = reshape(U(mesh.tri), [nelem, nlp]);
    xx = shp * xtri';
    uu = shp * utri';
    % flatten to 1D array without repeated nodes
    xxx = unique(xx);
    uuu = zeros(size(xxx));
    for i = 1:numel(xxx)
        [i1, i2] = find(xx == xxx(i));
        uuu(i) = uu(i1(1), i2(1));
    end
    % draw waves
    hold on
    plot(xxx, real(uuu))
    labels{j} = ['$k=' num2str(k_arr(j), '%.4f') '$'];
end
hold off
legend(labels, location='eastoutside')
fig_width = 19 * numel(labels);
set(fig_soln, 'position', [250 50 fig_width fig_width])
if ~exist('N', 'var')
    N = size(A1_N_max, 1);
end
title(['1D Helmholtz Solutions, $A=' num2str(mag_inc_wave) ...
    '$; $n_\mathrm{train}=' num2str(n_train) '$, $N=' num2str(N) '$'])
xlabel('$x \in \Omega$')
ylabel('$u_h(x; k)$')
