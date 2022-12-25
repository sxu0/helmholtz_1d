%% double checking my derivations :-)
A = 1;
nn = linspace(-4, 4, 4001);
kk = pi .* nn;

real_part = -4 * A^2 .* kk .* sin(3 .* kk) .* cos(kk);
imag_part = 4 * A^2 .* kk .* cos(3 .* kk) .* cos(kk);

figure
plot(nn, real_part, nn, imag_part)
legend('real', 'imag')
title('$a(u, u; k)$')

a = 4i * A^2 .* kk .* exp(3i .* kk) .* cos(kk);
real_part = real(a);
imag_part = imag(a);

figure
plot(nn, real_part, nn, imag_part)
legend('real', 'imag')
title('$a(u, u; k)$')
