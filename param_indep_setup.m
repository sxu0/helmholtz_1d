% helmholtz_fe_setup.m
% Assembles the parameter-independent stiffness matrices & load vectors.

function [A1, A2, A3, F1] = param_indep_setup(mag_inc_wave, p, nelem)

% arguments
% =========
% mag_inc_wave (float): nondimensionalized magnitude of incident wave
% p (int): polynomial degree of FE approximation
% nelem (int): number of elements in FE discretization
% 
% outputs
% =======
% A1 (matrix): parameter-independent stiffness matrix 1
% A2 (matrix): parameter-independent stiffness matrix 2
% A3 (matrix): parameter-independent stiffness matrix 3 (Q_a = 3 here)
% F (vector): parameter-independent load vector (Q_f = 1 here)

% number of gaussian quadrature points
pquad = 2 * p;

% setup reference element
ref = mkref1d(p, pquad);

% make mesh
mesh = mkmesh1d(nelem);
mesh = meshlin2p1d(mesh, ref);
mesh.bgrp = mkbgrp1d(mesh);

% number of shape functions per element
nshp = size(ref.shp, 2);

% compute and store local matrices
amat1 = zeros(nshp, nshp, nelem);
amat2 = zeros(nshp, nshp, nelem);
imat = zeros(nshp, nshp, nelem);
jmat = zeros(nshp, nshp, nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';

    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    jacq = ref.shpx(:,:,1) * xl;
    detJq = jacq;
    ijacq = 1 ./ jacq;

    % compute quadrature weight
    wqJ = ref.wq .* detJq;

    % compute basis
    phiq = ref.shp;
    phixq = bsxfun(@times, ref.shpx(:,:,1), ijacq);

    % compute stiffness matrices
    aaloc1 = phixq(:,:,1)' * diag(wqJ) * phixq(:,:,1);  % Theta_a1 = 1
    aaloc2 = phiq(:,:,1)' * diag(wqJ) * phiq(:,:,1);  % Theta_a2 = -k^2

    % prepare element-to-node connectivity map
    amat1(:,:,elem) = aaloc1;
    amat2(:,:,elem) = aaloc2;
    imat(:,:,elem) = repmat(tril, [1, nshp]);
    jmat(:,:,elem) = repmat(tril', [nshp, 1]);
end

% assemble param-indep stiffness matrices
ndof = size(mesh.coord, 1);  % p*nelem+1
A1 = sparse(imat(:), jmat(:), amat1(:), ndof, ndof);  % Theta_a1 = 1
A2 = sparse(imat(:), jmat(:), amat2(:), ndof, ndof);  % Theta_a2 = -k^2
A3 = sparse(1, 1, 1i, ndof, ndof);  % Theta_a3 = k

% set up param-indep load vector
F1 = zeros(ndof, 1);
F1(1, 1) = -2i * mag_inc_wave;  % Theta_f1 = k
end
