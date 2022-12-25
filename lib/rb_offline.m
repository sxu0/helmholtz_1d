function [A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mesh, ref] = rb_offline( ...
    n_train, err_tol, mag_inc_wave, p, nelem ...
    )
% Offline procedure for reduced-basis solver, including
% producing snapshots & setting up parameter-independent
% stiffness matrices/load vectors.
%
% inputs
% ------
% n_train (int): number of snapshots to produce
% err_tol (float): error threshold for truncating tailing eigenmodes
% mag_inc_wave (float): nondimensionalized magnitude of incident wave
% p (int): polynomial degree of FE Lagrange basis
% nelem (int): number of elements in FE discretization
%
% outputs
% -------
% A1_N_max, A2_N_max, A3_N_max (sparse matrices of floats): {N_max by
%   N_max} (hierarchical) parameter-independent RB stiffness matrices
% F1_N_max (vector of floats): {N_max by 1} (hierarchical) parameter
%   -independent RB load vector
% Z_N_max (matrix of floats): {ndof by N_max} (hierarchical in N) reduced
%   basis matrix
% mesh: FE mesh object from femmat
% ref: FE reference element object from femmat

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
[Z_N_max, mesh, ref] = snapshots( ...
    n_train, err_tol, p, nelem, false, mag_inc_wave ...
    );

% assemble parameter-independent stiffness matrices & load vector
[A1, A2, A3, F1] = param_indep_setup(mag_inc_wave, p, nelem);
% note: Theta_a1 = 1, Theta_a2 = -k^2, Theta_a3 = k, Theta_f1 = k
% thus A = A1 - k^2 .* A2 + k .* A3, F = k .* F1

% project A1, A2, A3, and F1 onto RB
A1_N_max = Z_N_max' * A1 * Z_N_max;
A2_N_max = Z_N_max' * A2 * Z_N_max;
A3_N_max = Z_N_max' * A3 * Z_N_max;
F1_N_max = Z_N_max' * F1;

end
