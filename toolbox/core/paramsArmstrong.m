function [ks_d, ks_p] = paramsArmstrong(Vmax, kcat, n, r, A, D, u)

% Compute parameters of substrate uptake following Armstrong's approximation to
% Pasciak and Gavis' quadratic equation.
%%
% Inputs
%
% Vmax     -       Double. Porter-limited limit of uptake rate for which 
%                  the number of transporters is known. [moles cell-1 s-1] 
% kcat     -       Double. Catalytic constant. [moles Tp-1 s-1]
% n        -       Double. Number of transporters corresponding to the porter-limited
%                  limit. [cell-1]
% r        -       Double. Cell radius. [m]
% A        -       Double. Patch area for the transporter determined by molecular
% modeling of protein sequence. [m2]
% D        -       Double. Diffusive constant for the substrate to be 
%                  transported. Determined at a known temperature and 
%                  salinity using getDiffusion(). [m2 s-1]
% u        -       Double. Shear velocity. [m s-1]

% Outputs
%
% v        -       Double. Uptake rate. [moles cell-1 s-1]
% ks_dl    -       Double. Half saturation concentration at the diffusive
%                  limit. [moles m-3]
% ks_pl    -       Double. Half saturation concentration at the porter
%                  limit. [moles m-3]





%% Solve for ks at the diffusive limit
ks_d = (n.*kcat)./ (r .*D); % moles m-3

%% Solve for ks at the porter limit
% mass transfer coefficient
mtc = D./r + u/2; % m s-1
% capture probability
alpha = ( (mtc).*sqrt(pi()*A) ) ./ (4*D); % dimensionless

ks_p = ( (pi() .* kcat) ./ (4*alpha .* A .* D) ) .* sqrt(A./pi()); % moles m-3



end
