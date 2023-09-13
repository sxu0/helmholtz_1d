function [uq, uxq] = uq_uxq(ref, nelem, elem, k, mag_inc_wave)

% assumes 1D domain
% (Gaussian quadrature order used at the moment is not very high)

h = (ref.xint(end) - ref.xint(1)) / nelem;

xq = ref.xint(1) + (elem - 1) * h + ref.xq * h;

uq = u(xq, k, mag_inc_wave);
uxq = ux(xq, k, mag_inc_wave);

function u = u(x, k, A)
    u = -2 * A .* exp(1i*k) .* (cos(k) .* cos(k.*x) + sin(k) .* sin(k.*x));
end

function ux = ux(x, k, A)
    ux = 2 * A .* k .* exp(1i*k) .* (cos(k) .* sin(k.*x) - sin(k) .* cos(k.*x));
end

end
