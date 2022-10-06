```MATLAB
% Parameters
L = 64; dt = 0.1; Nt = 16000; Ng = 256; Num_beams = 2; 
% Beam 1
N1 = 10000; V01 = 5; Vth1 = 1; QM1 = -1; XP1 = 0; Mode1 = 0; WP = 1;
% Beam 2
N2 = 10000; V02 = -5; Vth1 = 1; QM1 = -1; XP1 = 0; Mode1 = 0; 
% Ions in the background
IB = 20000;

% Methods
Motion_method = 'Euler method';  % 'Leapfrog' 'Runge Kutta (RK4)'
Field_method = 'Direct Integration';  % 'Finite Difference Method' 'Fast Fourier Transform (FFT)'
Interpolation_method = 'Cloud in Cell (CIC)'; % 'Nearest Grid Point (NGP)'

% Size of the cell
dx = L/(Ng-1);                                      

% Charge
[Q1,Q2,rho_back] = Charge(QM1,QM2,IB,N1,N2,L);

% Initial loading
[xp1,vp1] = InitialLoading(N1,V01,Vth1,XP1,Mode1); 
[xp2,vp2] = InitialLoading(N2,V02,Vth2,XP2,Mode2);

% Auxiliarity vectors
p1 = AuxVector(N1); p2 = AuxVector(N2); 

% Poisson equation preparation
un = ones(Ng-1, 1); % Ng-1 * 1
Ax = spdiags([un -2*un un], [-1 0 1], Ng-1, Ng-1); % Matrix that gives the indices of the Poisson equation Ng-1 * Ng-1.
kap = (2*pi/L)*(-Ng/2:Ng/2-1);
kap = fftshift(kap');
kap(1) = 1;

% Computational cycle
for it = 1:Nt
    switch Motion_method
        case 'Leapfrog'
            vp1 = MotionV(p1,QM1,mat1,Eg,it,N1);
            vp2 = MotionV(vp2,QM2,mat2,Eg,it,N2);
    end

    %Updating positions
    xp1 = MotionX(xp1,vp1);
    xp2 = MotionX(xp2,vp2);

    %ly periodic boundary conditions:
    xp1 = PeriodicBC(xp1,0,L);
    xp2 = PeriodicBC(xp2,0,L);

    %Interpolation functions
    mat1 = InterpolationF(xp1,N1,p1);
    mat2 = InterpolationF(xp2,N2,p2);

    % Charge density:
    rho1 = Charge_density(Q1,mat1,dx);
    rho2 = Charge_density(Q2,mat2,dx);
    rhot = rho1 + rho2 + rho_back;

    % Field equations:
    [Phi,Eg] = Field(,rhot);

    %Updating velocity
    vp1 = MotionV(vp1,QM1,mat1,Eg,it,N1);
    vp2 = MotionV(vp2,QM2,mat2,Eg,it,N2);
end
```
