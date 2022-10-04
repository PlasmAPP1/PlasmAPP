# ION-ACOUSTIC WAVES

Use the following parameters:
$\texttt{L}=2048$, $\texttt{Nt} = 8000$, $\texttt{dt}= 0.05$, $\texttt{Ng} =8192$, $\texttt{N1}=8000$, $\texttt{N2}=8000$, $\texttt{Ni}=0$,
$\texttt{V01}=0$, $\texttt{V02}=0$, $\texttt{Vth1}=1$, $\texttt{Vth2}=0$, $\texttt{QM1}=-1$, $\texttt{QM2}=0.01$, $\texttt{XP1}=0.2$, $\texttt{XP2}=0$, $\texttt{M1}=1$ y $\texttt{M2}=0$.
Then run the simulation and save the parameters.

The dispersion relation for Langmuir waves is:
$$\omega^2 =\pm\omega_{pe} \left(1+\frac{3}{2}k^2\lambda_D^2\right)$$
and for ion-acoustic waves is:
$$\frac{\omega^2}{k^2} = \frac{k_BT_e}{m_i}\left(\frac{1}{1+k^2\lambda_D^2}\right)+\frac{3k_BT_i}{m_i}$$

Load the results in MATLAB and plot the theoretical disperion relations over the simulated one. You must obtain:

![DRteoexp (2)](https://user-images.githubusercontent.com/114958650/193869357-64065141-127f-41a7-9821-e2d18fb89394.png)
