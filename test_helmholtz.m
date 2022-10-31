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


%%

