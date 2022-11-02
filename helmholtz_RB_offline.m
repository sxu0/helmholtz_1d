function [A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = helmholtz_RB_offline( ...
    n_train, err_tol, mag_inc_wave, p, nelem ...
    )

% incoming wave magnitude (boundary condition)
if ~exist('mag_inc_wave', 'var')
    mag_inc_wave = 1;
end

% FE specifications
if ~exist('p', 'var')
    p = 2;
end
if ~exist('nelem', 'var')
    nelem = 1000;
end

% compute reduced basis
[Z_N_max, mesh, ref] = helmholtz_snapshots( ...
    n_train, err_tol, p, nelem, false, mag_inc_wave ...
    );

% assemble parameter-independent stiffness matrices & load vector
[A1, A2, A3, F1] = helmholtz_off_on_setup(mag_inc_wave, p, nelem);
% note: Theta_a1 = 1, Theta_a2 = -k^2, Theta_a3 = k, Theta_f1 = k
% thus A = A1 - k^2 .* A2 + k .* A3, F = k .* F1

% project A1, A2, A3, and F1 onto RB
A1_N_max = Z_N_max' * A1 * Z_N_max;
A2_N_max = Z_N_max' * A2 * Z_N_max;
A3_N_max = Z_N_max' * A3 * Z_N_max;
F1_N_max = Z_N_max' * F1;

end
