function anim_soln = visualize(k_arr, u_arr, mesh, ref, mag_inc_wave, n_train, N, complex_part)
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
% complex_part (str, optional): choose 'real' or 'imag' part to plot,
%   defaults to real
%
% outputs
% -------
% anim_soln (figure): animated solution figure

figure;
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
    if ~exist('complex_part', 'var') || strcmp(complex_part, 'real')
        anim_soln = plot(xxx, real(uuu));
        ymin = min(min(real(u_arr)));
        ymax = max(max(real(u_arr)));
    elseif strcmp(complex_part, 'imag')
        anim_soln = plot(xxx, imag(uuu));
        ymin = min(min(imag(u_arr)));
        ymax = max(max(imag(u_arr)));
    end
    yrange = ymax - ymin;
    ylim([ymin - 0.1 * yrange, ymax + 0.1 * yrange])
    legend(['$k=' num2str(k_arr(j), '%.4f') '$'])
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
    if ~exist('complex_part', 'var') || strcmp(complex_part, 'real')
        title(['1D Helmholtz Solution Real Part' BC_config RB_params])
    elseif strcmp(complex_part, 'imag')
        title(['1D Helmholtz Solution Imaginary Part' BC_config RB_params])
    end
    xlabel('$x \in \Omega$')
    ylabel('$u_\mathrm{RB}(x; k)$')
    pause(0.1)
end

end
