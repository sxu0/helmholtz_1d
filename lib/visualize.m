function anim_soln = visualize(k_arr, u_arr, mesh, ref, mag_inc_wave, n_train, N, gif_name)
% Visualizes real/imag parts of solution, animated by wavenumber.
%
% inputs
% ------
% k_arr (array of floats): wavenumbers (param values) targeted by solution
% u_arr (array of floats): soln FE coefficients at each wavenumber;
%   {ndof by numel(k_arr)} array
% mesh: FE mesh object from femmat
% ref: FE reference element object from femmat
% mag_inc_wave (float, optional): nondimensionalized incoming wave
%   magnitude
% n_train (int, optional): number of snapshots produced
% N (int, optional): dimension of reduced basis
% filename (str, optional): file path for saving output animated figure,
%   does not save if '' or unspecified
%
% outputs
% -------
% anim_soln (figure): animated solution figure

for j = 1:numel(k_arr)
    U = u_arr(:, j);
    % next 8 lines are sourced from `femmat/util/plot/plotfield1d()`
    nelem = size(mesh.tri, 1);
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
    if j == 1
        anim_soln = figure(2);
    end
    figure(2)
    plot(xxx, real(uuu));
    hold on
    plot(xxx, imag(uuu));
    hold off
    ymin = min(min(min(real(u_arr))), min(min(imag(u_arr))));
    ymax = max(max(max(real(u_arr))), max(max(imag(u_arr))));
    yrange = ymax - ymin;
    ylim([ymin - 0.1 * yrange, ymax + 0.1 * yrange])
    legend('real part', 'imag part')
    text(0.055, 0.920, ['$k=' num2str(k_arr(j), '%.4f') '$'], units='normalized')
    if exist('mag_inc_wave', 'var')
        BC_config = [': $A=' num2str(mag_inc_wave) '$'];
    else
        BC_config = '';
    end
    if exist('n_train', 'var') && exist('N', 'var')
        RB_params = ['; $n_\mathrm{train}=' num2str(n_train) '$, $N=' num2str(N) '$'];
    elseif exist('n_train', 'var')
        RB_params = ['; $n_\mathrm{train}=' num2str(n_train) '$'];
    elseif exist('N', 'var')
        RB_params = ['; $N=' num2str(N) '$'];
    else
        RB_params = '';
    end
    title(['1D Helmholtz Solution' BC_config RB_params])
    xlabel('$x \in \Omega$')
    ylabel('$u_\mathrm{RB}(x; k)$')
    % save figure to gif
    if exist('gif_name', 'var') && ~isempty(gif_name)
        if j == 1 && isfile(gif_name)
            delete(gif_name)
        end
        exportgraphics(anim_soln, gif_name, 'Append', true);
    end
    pause(0.025)
end

end
