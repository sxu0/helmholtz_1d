function [u_err_NN, q_err_NN] = error_rb(NN, kk_test, n_train, err_tol, mag_inc_wave, p, nelem, plt)
% Studies error convergence by N (RB dimension).
% Error is defined as the maximum of those evaluated at each k.
% H1 norm is taken for solution error, and absolute value for output error.
%
% inputs
% ------
% NN (array of ints): array of N (RB dimension) at which to observe error
% kk_test (array of floats): array of k (nondimensionalized wavenumber)
%   at which solution is of interest
% n_train (int): number of snapshots to produce
% err_tol (float): error threshold for truncating tailing eigenmodes
% p (int): polynomial degree of FE Lagrange basis
% nelem (int): number of FE mesh elements in FE discretization
% plt (bool, optional): whether to plot computed errors by dof;
%   defaults to false
%
% outputs
% -------
% u_err_NN (array of floats): max u_N error for each N
% q_err_NN (array of floats): max q_N error for each N

if ~exist('plt', 'var')
    plt = false;
end

% offline steps
[offline.A1_N_max, offline.A2_N_max, offline.A3_N_max, ...
    offline.F1_N_max, offline.Z_N_max] = rb_offline(n_train, err_tol);

% online steps for each N
u_err_NN = zeros(size(NN));
q_err_NN = zeros(size(NN));
for i = 1:numel(NN)
    [u_err_NN(i), q_err_NN(i)] = max_error( ...
        NN(i), offline, kk_test, mag_inc_wave, p, nelem);
end

if plt
    plot_errors(NN, u_err_NN, q_err_NN, n_train)
end

function plot_errors(NN, u_err_N, q_err_N, n_train)
    figure;
    yyaxis left
    semilogy(NN, u_err_N, 'o')
    ylabel(['$\mathrm{max}_{k \in \Xi_\mathrm{test}} ' ...
        '\|u_h(k)-u_N(k)\|_{H^1(\Omega)}$'])
    ylim('padded')
    yyaxis right
    semilogy(NN, q_err_N, 'o')
    ylabel(['$\mathrm{max}_{k \in \Xi_\mathrm{test}} ' ...
        '|q_h(k)-q_N(k)|$'])
    ylim('padded')
    xlabel('N')
    xlim('padded')
    title(['RB error convergence, ' ...
        '$n_{\mathrm{train}}=' num2str(n_train) '$'])
end

function [u_err_max, q_err_max] = max_error(N, offline, kk_test, mag_inc_wave, p, nelem)
    % RB soln & output
    [~, q_RB_kk, u_RB_kk] = rb_online( ...
        offline.A1_N_max, offline.A2_N_max, offline.A3_N_max, ...
        offline.F1_N_max, offline.Z_N_max, ...
        mag_inc_wave, kk_test, N);

    % compute error for each k
    u_errr = zeros(size(kk_test));
    q_errr = zeros(size(kk_test));
    for j = 1:numel(kk_test)
        % FE soln
        [u_h, mesh, ref] = fe_solver(kk_test(j), mag_inc_wave, p, nelem);

        % FE output
        q_h = quantity_of_interest(u_h, mag_inc_wave);

        % plot FE against RB
        if 0
            [xx_h, uu_h] = higher_res(u_h, mesh, ref, 10);
            [~, uu_RB] = higher_res(u_RB_kk(:, j), mesh, ref, 10);
            figure
            plot(xx_h, uu_h, '-', xx_h, uu_RB, ':')
            legend("$u_h$", "$u_N$")
            xlabel("$x$")
            ylabel("$u$")
        end

        % compute errors
        u_errr(j) = hnorm_1d(u_RB_kk(:,j) - u_h, mesh, ref);
        q_errr(j) = abs(q_RB_kk(:,j) - q_h);
    end

    % deterine max errors
    u_err_max = max(u_errr);
    q_err_max = max(q_errr);
end

function hnorm = hnorm_1d(w, mesh, ref)
    % note gradient() & trapz() are not exact; hence 10x resolution is used
    [xx, ww] = higher_res(w, mesh, ref, 10);
    dx = xx(2) - xx(1);
    wx = gradient(ww, dx);
    hnorm = sqrt(trapz(xx, wx .* wx + ww .* ww));
end

function [xx_h, uu_h] = higher_res(u_h, mesh, ref, res)
    % reorders the FE coeff vector & refines soln resolution as desired
    % code adapted from femmat/util/plot/plotfield1d.m
    n_elem = size(mesh.tri, 1);
    nlp = size(ref.shp, 2);
    xxplot = linspace(0, 1, res*mesh.p)';
    shp = shape1d(ref.p, ref.xint, xxplot);
    xtri = reshape(mesh.coord(mesh.tri), [n_elem, nlp]);
    utri = reshape(u_h(mesh.tri), [n_elem, nlp]);
    xx_h = reshape(shp*xtri', 1, []);  % flatten to row vector
    % FE soln @ {res}x resolution
    uu_h = reshape(shp*utri', 1, []);  % flatten to row vector
end

end
