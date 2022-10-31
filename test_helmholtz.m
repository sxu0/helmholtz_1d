clear; close all;

%% testing helmholtz_FE_solver() function

mag_inc_wave = 1;
p = 2;
nelem = 1000;
plt = true;

figure
for k = (11:60)/10

[U, ~, ~] = helmholtz_FE_solver(k, mag_inc_wave, p, nelem, plt);

end
hold off


%%
