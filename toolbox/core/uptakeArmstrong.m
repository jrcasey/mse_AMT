function [v] = uptakeArmstrong(S, Vmax, ks_d, ks_p)

%% Compute uptake rate of a substrate following Armstrong's approximation to
% Pasciak and Gavis' quadratic equation.

% Inputs
%
% S        -       Double. Substrate concentration. [moles m-3]
% Vmax     -       Double. Porter-limited limit of uptake rate for which 
%                  the number of transporters is known. [moles cell-1 s-1] 
% ks_dl    -       Double. Half saturation concentration at the diffusive
%                  limit. [moles m-3]
% ks_pl    -       Double. Half saturation concentration at the porter
%                  limit. [moles m-3]

% Outputs
%
% v        -       Double. Uptake rate. [moles cell-1 s-1]



%% Solve for uptake rate

v = (Vmax .* S)./(ks_d + ks_p + S);

end
