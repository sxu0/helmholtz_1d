% helmholtz_FE_solver.m
% 1D Helmholtz equation.
% Using femmat library in UTIAS ACEL gitlab.
% Shiqi Xu

function [U, mesh, ref] = helmholtz_FE_solver(k, mag_inc_wave, p, nelem, plt)
% arguments
% =========
% k (float): nondimensionalized wavenumber, on (0, inf]
% mag_inc_wave (float): nondimensionalized magnitude of incident wave
% p (int): polynomial degree of FE basis
% nelem (int): number of elements over physical domain
% plt (bool, optional): whether to generate plots, defaults to false
%
% outputs
% =======
% U (matrix(float)): nodal basis coefficients
% mesh: FE mesh object
% ref: FE reference element object

if ~exist('plt', 'var')
    plt = false;
end

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
amat_H1 = zeros(nshp, nshp, nelem);
imat = zeros(nshp, nshp, nelem);
jmat = zeros(nshp, nshp, nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';  % each row is indices of nodes on local element

    % compute mesh jacobians
    xl = mesh.coord(tril,:);  % local mesh coords
    jacq = ref.shpx(:,:,1) * xl;  % jacobian (eval. at each quadrature pt)
    detJq = jacq;  % det(J) (eval. at each quadrature pt)
    ijacq = 1 ./ jacq;  % inverse jacobian (eval. at each quadrature pt)

    % compute quadrature weight
    wqJ = ref.wq .* detJq;  % incl. detJq from area transformation (dx)

    % compute basis
    phiq = ref.shp;  % shape functions at quadrature points of physical element
    % shape functions same in ref & physical domains
    phixq = bsxfun(@times, ref.shpx(:,:,1), ijacq);  % derivs of shape fns at quadrature pts of physical elem
    % ref.shpx is derivatives of shape functions on reference interval
    % bsxfun() broadcasts & performs element-wise operation

    % compute local stiffness matrix
    aaloc = phixq(:,:,1)' * diag(wqJ) * phixq(:,:,1) ...
        - k^2 .* phiq(:,:,1)' * diag(wqJ) * phiq(:,:,1);

    % prepare element-to-node connectivity map
    amat_H1(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril, [1, nshp]);
    jmat(:,:,elem) = repmat(tril', [nshp, 1]);
end

% assemble H1 part of global stiffness matrix
ndof = size(mesh.coord, 1);  % p*nelem+1
A_H1 = sparse(imat(:), jmat(:), amat_H1(:), ndof, ndof);
% assemble BC part of global stiffness matrix
A_BC = sparse(1, 1, 1i * k, ndof, ndof);
% put together global stiffness matrix
A = A_H1 + A_BC;

% set up global load vector
F = zeros(ndof, 1);
F(1, 1) = -2i * k * mag_inc_wave;

% solve linear system
U = A \ F;

% plot solution
if plt == true
    figure(1), clf,
    plotfield1d(mesh, ref, U);
end

end
