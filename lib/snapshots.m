function [Z_N_max, mesh, ref] = snapshots(kk_train, epsilon_tol, p, nelem, plt, mag_inc_wave)
% Produce snapshots sampled across parameter space, solved by FEM.
% Outputs reduced basis matrix, containing reduced basis functions.
%
% inputs
% ------
% kk_train (array of floats): training points sampled in parameter space;
%   a snapshot is produced for each k
% epsilon_tol (float): error threshold for truncating tailing eigenmodes
% p (int, optional): polynomial degree of FE Lagrange basis, defaults to 2
% nelem (int, optional): number of FE mesh elements in FE discretization,
%   defaults to 1000
% plt (bool, optional): whether to generate plots, defaults to false
% mag_inc_wave (float, optional): nondimensionalized incoming wave
%   magnitude, defaults to 1
%
% outputs
% -------
% Z_N_max: {ndof by N_max} reduced basis matrix
% mesh: FE mesh object
% ref: FE reference element object

if ~exist('plt', 'var')
    plt = false;
end

if ~exist('mag_inc_wave', 'var')
    mag_inc_wave = 1;
end

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

n_train = numel(kk_train);  % number of snapshots to produce

% FE specifications
if ~exist('p', 'var')
    p = 2;
end
if ~exist('nelem', 'var')
    nelem = 1000;
end
N_t = p * nelem + 1;  % dim of "truth" FE approx

%% solve problem for each mu
U_store = zeros(N_t, n_train);
for n = 1:n_train
    k = kk_train(n);
    % U are {p*nelem+1 by 1} vectors
    if n == 1
        [U_n, mesh, ref] = fe_solver(k, mag_inc_wave, p, nelem);
    else
        [U_n, ~, ~] = fe_solver(k, mag_inc_wave, p, nelem);
    end
    U_store(:, n) = U_n;
end

% assemble H1 inner product matrix (over spatial domain)
X_H1 = h1_inprod_matrix(p, nelem);

% assemble correlation matrix
C_POD = zeros(n_train, n_train);
for i = 1:n_train
    U_N_i = U_store(:, i);
    for j = 1:i
        U_N_j = U_store(:, j);
        C_POD(i, j) = 1/n_train * (U_N_i' * X_H1 * U_N_j);
    end
end
% build symmetric matrix from lower left triangular part
C_POD = C_POD + C_POD' - eye(size(C_POD)) .* C_POD;

%% solve for normalized eigenvectors
[psis, lambdas] = eig(C_POD);
% sort eigenvalues largest -> smallest
% (note eig returns eigenvalues in an arbitrary order!)
[~, ind] = sort(diag(lambdas), 'descend');
lambdas = lambdas(ind, ind);  % dimensions {n_train by n_train}
lambdas_ordered = diag(lambdas);
if any(imag(lambdas_ordered) > 1e-4)  % hardcoded threshold at the moment
    warning('C_POD not Hermitian!')
else
    % drop tiny imaginary parts from floating point error
    lambdas_ordered = real(lambdas_ordered);
end
psis = psis(:, ind);

% normalize eigenvectors such that their inner product is n_train * lambda_N
psis = psis ./ sqrt(repmat(lambdas_ordered.', n_train, 1) .* n_train);

%% determine basis functions to truncate based on minimum error tolerance

% calculate error up to each basis function
epsilons_POD = sqrt(cumsum(abs(lambdas_ordered), 'reverse'));  % L2 norm in param
% format longg
% epsilons_POD

Ns = 1:numel(lambdas_ordered);

if strcmp(plt, 'lin')
    figure
    yyaxis left
    plot(Ns, lambdas_ordered, 'o')
    ylabel("$\lambda^{\mathrm{POD},\,N}$")
    yyaxis right
    plot(Ns, epsilons_POD, '.')
    ylabel("$\bar{\bar{\epsilon}}_N^\mathrm{POD}$")
    xlabel("$N$")
    legend('Eigenvalue', 'Truncation Error up to $N$th Eigenmode')
    title(['Sorted Eigenvalues \& Associated Truncation Errors, ' ...
        '$n_\mathrm{train}=' num2str(n_train) '$'])
    xlim('padded')
    ylim('padded')
    set(gcf, 'position', [300 300 640 480])
elseif strcmp(plt, 'log')
    figure  % there is a kink!
    yyaxis left
    semilogy(Ns, lambdas_ordered, 'o')
    ylabel("$\lambda^{\mathrm{POD},\,N}$")
    yyaxis right
    semilogy(Ns, epsilons_POD, '.')
    ylabel("$\bar{\bar{\epsilon}}_N^\mathrm{POD}$")
    xlabel("$N$")
    legend('Eigenvalue', 'Truncation Error up to $N$th Eigenmode')
    title(['Sorted Eigenvalues \& Associated Truncation Errors, ' ...
        '$n_\mathrm{train}=' num2str(n_train) '$'])
    xlim('padded')
    ylim('padded')
    set(gcf, 'position', [300 300 640 480])
end

% find N_max within specified tolerance
N_pass = Ns(epsilons_POD < epsilon_tol);
if isempty(N_pass)
    error("error tolerance tighter than achievable!")
end
N_max = min(N_pass);

%% construct basis functions, store in N_max matrix (hierarchical)
Z_N_max = U_store(:, :) * psis(:, 1:N_max);

end
