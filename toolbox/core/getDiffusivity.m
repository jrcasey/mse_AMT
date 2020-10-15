function D = getDiffusivity(MW,T)
% Calculates molecular diffusivity in seawater.

T_k = T+273.18; % temperature in kelvin
exoMetVol = MW .* (1/0.56); % angstroms cubed
% convert to cm3 mol-1
exoMetVol2 = exoMetVol * (1/1e24) * 6.022e23; % cm^3 mol^-1
A = -3.7188;
B = 579.919;
C = -137.546;
eta = exp(A + (B./(C+T_k))); % mPa*s OR centipoise
D = 1e-4 * 13.26e-5 ./ ((eta.^1.4).*(exoMetVol2.^0.589))'; %m2 s-1
end
