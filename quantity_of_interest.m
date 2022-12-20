function q_arr = quantity_of_interest(u_arr, mag_inc_wave)
% Computes the quantity of interest given the solution
% FE coefficients and the magnitude of the incoming wave.
%
% inputs
% ------
% u_arr (array of floats): soln FE coefficients at each wavenumber;
%   {ndof by numel(k_arr)} array
% mag_inc_wave (float): nondimensionalized incoming wave magnitude
%
% outputs
% -------
% q_arr (array of floats): quantity of interest (output) at each k;
%   {1 by numel(k_arr)} array

q_arr = 1/2 .* abs(u_arr(1,:) - mag_inc_wave) .^ 2;

end
