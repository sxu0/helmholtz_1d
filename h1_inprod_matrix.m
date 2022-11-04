function X = h1_inprod_matrix(p, nelem)
% Assembles sparse H1 inner product matrix.
%
% inputs
% ------
% p (int): polynomial degree of FE Lagrange basis
% nelem (int): number of elements in FE discretization
%
% outputs
% -------
% X (sparse matrix of floats): {ndof by ndof} H1 inner product matrix

% number of quadrature points
pquad = 2 * p;

% setup reference element
ref = mkref1d(p, pquad);

% make mesh
mesh = mkmesh1d(nelem);
mesh = meshlin2p1d(mesh, ref);
mesh.bgrp = mkbgrp1d(mesh);

% number of shape functions per element
nshp = size(ref.shp, 2);

% compute & store local matrices
Xmat = zeros(nshp, nshp, nelem);
imat = zeros(nshp, nshp, nelem);
jmat = zeros(nshp, nshp, nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';

    % compute mesh jacobians
    xl = mesh.coord(tril,:);  % local mesh coords
    jacq = ref.shpx(:,:,1) * xl;  % 1d jacobian
    detJq = jacq;  % 1d determinant of jacobian
    ijacq = 1 ./ jacq;  % 1d inverse jacobian

    % compute quadrature weight
    wqJ = ref.wq .* detJq;  % detJq from area transformation (dx) included

    % compute basis
    phiq = ref.shp;  % shape functions at quadrature points of physical element
    phixq = bsxfun(@times, ref.shpx(:,:,1), ijacq);  % derivs of shape fns at quadrature pts of physical elem

    % compute stiffness matrix
    Xaloc = phixq(:,:,1)' * diag(wqJ) * phixq(:,:,1) ...
        + phiq(:,:,1)' * diag(wqJ) * phiq(:,:,1);

    % prep element-to-node connectivity map
    Xmat(:,:,elem) = Xaloc;
    imat(:,:,elem) = repmat(tril, [1, nshp]);
    jmat(:,:,elem) = repmat(tril', [nshp, 1]);
end

% assemble matrix
ndof = size(mesh.coord, 1);  % p*nelem+1
X = sparse(imat(:), jmat(:), Xmat(:), ndof, ndof);

end
