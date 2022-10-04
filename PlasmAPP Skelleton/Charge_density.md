```MATLAB
%% Function that calculates the charge density
function rho = Charge_density(charge,interp,dx)
    rho = full((charge/dx)*sum(interp))';
end
