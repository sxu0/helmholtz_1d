function [u_N_arr, s_N_arr, u_arr] = rb_online( ...
    A1_N_max, A2_N_max, A3_N_max, F1_N_max, Z_N_max, mag_inc_wave, ...
    mu_arr, N ...
    )

N_max = size(A1_N_max, 1);
if ~exist('N', 'var')
    N = N_max;  % default to N_max accuracy
elseif N > N_max
    error("N cannot be larger than N_max")
end
% all these are hierarchical!
A1_N = A1_N_max(1:N, 1:N);
A2_N = A2_N_max(1:N, 1:N);
A3_N = A3_N_max(1:N, 1:N);
F1_N = F1_N_max(1:N);
Z_N = Z_N_max(:, 1:N);

% storage for soln RB coefficients
u_N_arr = zeros(N, size(mu_arr, 2));

for i = 1:size(mu_arr, 2)
    k = mu_arr(i);

    % assemble RB stiffness matrix
    A_N = A1_N - k^2 .* A2_N + k .* A3_N;

    % assemble RB load vector
    F_N = k .* F1_N;

    % solve for soln RB coefficients
    u_N = A_N \ F_N;

    % store u_N
    u_N_arr(:, i) = u_N;
end

% construct soln FE coefficients
u_arr = Z_N * u_N_arr;

% generate output
s_N_arr = 1/2 .* abs(mag_inc_wave - u_arr(1,:)) .^2;

end
