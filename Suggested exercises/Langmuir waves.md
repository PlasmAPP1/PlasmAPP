# LANGMUIR WAVES

Simulate Langmuir waves using the following parameters:
L = 1024, Nt = 8000, dt = 0.05, Ng = 8192, N1 = 50000, Ni = 50000, Vth1 = 1, QM1 = -1, WP = 1, and set the other parameters to zero:

<img src="https://user-images.githubusercontent.com/114958650/193858000-48e58e54-d577-4a0d-9b60-7f41540ed2fd.png" data-canonical-src="https://user-images.githubusercontent.com/114958650/193858000-48e58e54-d577-4a0d-9b60-7f41540ed2fd.png" width="300" height="200" />

<img src="https://user-images.githubusercontent.com/114958650/193858314-46c6c8fd-11d1-4a5e-b046-fbfcf96cfbb0.png" data-canonical-src="https://user-images.githubusercontent.com/114958650/193858314-46c6c8fd-11d1-4a5e-b046-fbfcf96cfbb0.png" width="300" height="200" />

<img src="https://user-images.githubusercontent.com/114958650/193858578-7059b408-33ee-485b-b306-6327991296e9.png" data-canonical-src="https://user-images.githubusercontent.com/114958650/193858578-7059b408-33ee-485b-b306-6327991296e9.png" width="300" height="200" />

<img src="https://user-images.githubusercontent.com/114958650/193858910-8d829c3c-faff-43d5-92c4-6b57663943c8.png" data-canonical-src="https://user-images.githubusercontent.com/114958650/193858910-8d829c3c-faff-43d5-92c4-6b57663943c8.png" width="300" height="200" />

Run the program and save the results:
![image](https://user-images.githubusercontent.com/114958650/193861727-97d02f99-f273-4f7c-9012-1b0e9b8d4141.png)

Theoretically, for Langmuir waves the dispersion relation is:
$$\omega^2 =\pm\omega_{pe} \left(1+\frac{3}{2}k^2\lambda_D^2\right)$$

where $\omega$ corresponds to the frequency associated to the propagating fields and $k$ the corresponding wave-number.

Load the results into MATLAB and create a program that plots the theoretical dispersion relation over the simulated one:
```MATLAB
%% DISPERSION RELATION EXPERIMENTAL 
n=nextpow2(N_g);
n1=nextpow2(NT);
kmin= -pi/Dx; kmax=pi/Dx;
wmin=2*pi/(Dt*NT)*(-NT/2); wmax=2*pi/(Dt*NT)*(NT/2);
kaxis=linspace(kmin,kmax,65536); waxis=linspace(wmin,wmax,65536);
fftEg=fftshift(fft2(Eg_hist'));

%% THEORICAL DISPERSION RELATION 
k = linspace(-kmax, kmax , 65536);
w = sqrt(1 + 3 * Vth_1^2 * (k).^2);
m = sqrt(3)*Vth_1*k;

%% PLOT THEORICAL AND SIMULATED
figure;
imagesc(kaxis,waxis,abs(fftEg))
hold on
p = plot(k,w,Color='Red',LineWidth=0.5);
p3 = plot(k,ones(1,65536),'--',LineWidth=1, Color='m');
legend('\omega^2=\omega_{pe}^2+3k^2v_{th,e}^2','\omega=\omega_{pe}',FontName='Times')
axis xy
yticks(0.5:0.25:1.5)
ylabel('Angular frequency');
xlabel('Wavenumber');
ylim([0.5 1.5])
xlim([0 0.3])
c = colorbar;
ytickformat('%.2f'); xtickformat('%.1f'); ax = gca; ax.FontName ='Times'; ax.FontSize = 10;
hold off
```

You should obtain:
![image](https://user-images.githubusercontent.com/114958650/193866884-19154599-9afa-4e0b-859d-781d7841f53c.png)


