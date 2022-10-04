```MATLAB
% Function for field equations
function [Phi,Eg] = Field(charge_density)
    switch Field_method
        case 'Finite Difference Method'
            Phi = Ax\(-charge_density(1:Ng-1)*dx^2);
            Phi = [Phi;0];
            Eg = ([Phi(Ng); Phi(1:Ng-1)]-[Phi(2:Ng);Phi(1)])/(2*dx);
        case 'Fast Fourier Transform (FFT)'
            rho_k = fft(charge_density(1:Ng));
            phi_k = rho_k./kap.^2;
            Phi = real(ifft(phi_k,'symmetric'));
            Phi = [Phi;0];
            Eg_k = -1i.*kap.*phi_k;
            Eg =  real(ifft(Eg_k,'symmetric'));
        case 'Direct Integration'
            charge_density(1) = charge_density(1) + charge_density(Ng);
            charge_density(Ng) = charge_density(1);
            G = 0;
            for k = 1:(Ng-1)
                G = G + (charge_density(k+1)+charge_density(k))*dx/2;
                G1(k) = G;
            end
            G1 = [0,G1];

            F = 0;
            for k = 1:(Ng-1)
                F = F + (G1(k+1)+G1(k))*dx/2;
                F1(k) = F;
            end
            F1 = [0,F1];

            Phi_prime = (1/L)*F1(Ng);
            Eg = -Phi_prime +  G1;

            Eg = Eg';

            Phi1 = 0;
            for k = 1:(Ng-1)
                Phi1 = Phi1 + (Eg(k+1)+Eg(k))*dx/2;
                Phi(k) = -Phi1;
            end
            Phi = [Phi';0];
            Phi(Ng) = Phi(1);
            Eg(Ng) = Eg(1);
    end
end
