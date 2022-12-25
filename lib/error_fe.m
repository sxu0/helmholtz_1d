function [u_errr, q_errr] = error_fe(nnelem, k, mag_inc_wave, p, plt)
% Studies error convergence by dof for FE snapshots.
% Computes H1 norm of solution error, and absolute value of output error.
%
% inputs
% ------
% nnelem (vector of ints): array of nelem (number of finite elements)
% k (float): nondimensionalized wavenumber, restricted to (0, inf]
% mag_inc_wave (float): nondimensionalized magnitude of incident wave
% p (int): polynomial degree of FE Lagrange basis
% plt (bool, optional): whether to plot computed errors by dof;
%   defaults to false
%
% outputs
% -------
% u_errr (array of floats): H1 norm of u_h error for each nelem
% q_errr (array of floats): absolute value of q_h error for each nelem

% TO CHECK:
% (1) does error from real & complex parts need to be treated separately?
% (2) {p=2, nelem=256, k=2*pi, mag_inc_wave=1} produces q_err of 0 - can this be the case?

if ~exist('plt', 'var')
    plt = false;
end

u_errr = nan(size(nnelem));
q_errr = nan(size(nnelem));
for i = 1:numel(nnelem)
    u_errr(i) = u_err_h1_norm(nnelem(i), k, mag_inc_wave, p);
    q_errr(i) = q_err_abs(nnelem(i), k, mag_inc_wave, p);
end

if plt
    plot_errors(u_errr, q_errr, nnelem, p)
end

function plot_errors(u_err, q_err, nnelem, p)
    dof = p .* nnelem + 1;
    figure;
    yyaxis left
    loglog(dof, u_err, 'o')
    ylabel('$\|u_h - u\|_{H^1(\Omega)}$')
    ylim('padded')
    yyaxis right
    loglog(dof, q_err, 'o')
    ylabel('$|q_h - q|$')
    ylim('padded')
    xlabel('dof')
    xlim('padded')
    title(['FE error convergence, $p=' num2str(p) '$'])
end

function q_err = q_err_abs(nelem, k, mag_inc_wave, p)
    u_exact_inlet = exact_soln(0, k, mag_inc_wave);
    q_exact = 1/2 * abs(u_exact_inlet - mag_inc_wave) ^ 2;

    uu_h = fe_solver(k, mag_inc_wave, p, nelem);
    u_h_inlet = uu_h(1);
    q_h = 1/2 * abs(u_h_inlet - mag_inc_wave) ^ 2;

    q_err = abs(q_h - q_exact);
end

function u_err = u_err_h1_norm(nelem, k, mag_inc_wave, p)
    bilinear_exact = 4 .* mag_inc_wave .^ 2 .* k ...
        .* (1i .* sin(k) + cos(k)) .^ 2 ...
        .* cos(k) .* (1i .* cos(k) - sin(k));

    [uu_h, ~, ~, A] = fe_solver(k, mag_inc_wave, p, nelem);
    u_h_mat = repmat(uu_h, 1, size(uu_h, 1));
    bilinear_fe = sum(u_h_mat .* u_h_mat.' .* A, 'all');

    u_err = sqrt(bilinear_exact - bilinear_fe);
end

function uu = exact_soln(xx, k, mag_inc_wave)
    uu = -2 .* mag_inc_wave .* (1i .* sin(k) + cos(k)) ...
        .* (cos(k) .* cos(k .* xx) + sin(k) .* sin(k .* xx));
end

end
