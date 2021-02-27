%%%%%%%%This code solves voriticity-streamfunction formulation %%%%%%%%%%%%%%%%%
%of the 2D Navier Stokes equation in a box by implementing a pseudo spectral algorithm%%%%%%
%%%%%%%%with RK3-IMEX scheme in time integration : implicit for linear%%%%%%
%%%%%%%%%%%%%and explicit for non-linear terms%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
close all;
%%%%%%%%%%%CONSTANTS FOR RKW3 IMEX (Rogers et. al (JCP), 1991) scheme %%%%%%% 
A(1) = 29/96;A(2) = -3/40;A(3) = 1/6;
B(1) = 37/160;B(2) = 5/24;B(3) = 1/6;
G(1) = 8/15;G(2) = 5/12;G(3) = 3/4;
Z(1) = -17/60;Z(2) = -5/12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Re = 1e6;%Reynolds number
prompt = 'Enter the simulation time (50 in this case)';
t_end = input(prompt);
prompt = 'Enter the spatial grid points (1024 in this case)';
Nx = input(prompt);
Ny = Nx;
Lx = 1;%%%Domain in the x-direction%%%%%%%%%
Ly = 1;%%%Domain in the y-direction%%%%%%%%%
prompt = 'Enter the time step';
del_t = input(prompt);
x = (0:Nx-1)*Lx/Nx;
y = (0:Ny-1)*Ly/Ny;
[X,Y] = meshgrid(x,y);
k = [0:Nx/2-1, -Nx/2:-1]*(2*pi)/Lx;  %wavenumbers corresponding to x
l = [0:Ny/2-1, -Ny/2:-1]*(2*pi)/Ly;  %wavenumbers corresponding to y
k(1) = 1e-4;l(1) = 1e-4;
[K,L] = meshgrid(k,l);
k_til = K.*K + L.*L;% total wave number%%%%%%%%%%%%%%%
mu = 0.0;    %mean
sigma = 1.0; %standard deviation
j = 1;
m = 0;
%E_in = k_til.^(3/4);
%u_hat = sqrt(E_in)./mean(sqrt(E_in));
%v_hat = sqrt(E_in)./mean(sqrt(E_in));
%mu = 1e-5;
Re = 250*5e4%Reynolds number
%u_amp = Re*mu;
%u_avg = mean(u_hat);
%u_rms = rms(u_hat - u_avg);
%v_hat = u_amp*v_hat./u_rms;%v_hat(1,1) = 0.0+1i*0.0;
%u_hat = u_amp*u_hat./u_rms;%u_hat(1,1) = 0.0+1i*0.0;
%rms(u_hat - u_avg)
%omg = ifft2(omeg_hat);
%[hC, hC] = contourf(X, Y, real(omg), 50); axis equal; 
omega = wgn(Ny, Nx,5*sigma);
omeg_hat = fft2(omega);
u_hat = 1i*L.*omeg_hat./k_til;u_hat(1,1) = 0.0+1i*0.0;
v_hat = -1i*K.*omeg_hat./k_til;v_hat(1,1) = 0.0+1i*0.0;
max(rms(u_hat - mean(u_hat)))
%omeg_hat = k_til.*u_hat./(1i*L);
%omeg_hat = omeg_hat./(mean(omeg_hat));
E_in = 0.5*sum(sum(u_hat.*conj(u_hat)+v_hat.*conj(v_hat)))/(Nx*Ny);
E = E_in
fprintf('Total energy at t=0 is %f\n',E_in);
omg = ifft2(omeg_hat);
m = m + 1;
[hC, hC] = contourf(X, Y, real(omg), 50); axis equal
set(hC,'LineStyle','none'); colorbar; colormap redblue(100); shading interp;drawnow
saveas(hC,sprintf('FIG_ext%d.png',m));
fp = fopen('transient.dat', 'w');  %Saving Energy and Enstrophy at each time
time = 0.0;
while(time < t_end)
    time = time+del_t;
    fprintf('time: %f, delt: %12.8f, energy = %f\n', time, del_t, E/E_in);
    %**********************************************************************
    %************************FIRST SUBSTEP********************************
    %**********************************************************************
    [NL] = Convol(omeg_hat,Nx,Ny,K,L);
    omeg_hat = (omeg_hat + del_t*(-A(1)*(k_til.*omeg_hat)/Re + G(1)*NL))./(1.0 + B(1)*del_t*k_til/Re); 
    %**********************************************************************
    %************************SECOND SUBSTEP********************************
    %**********************************************************************
    [NLp] = Convol(omeg_hat,Nx,Ny,K,L);
    omeg_hat = (omeg_hat + del_t*(-A(2)*(k_til.*omeg_hat)/Re + G(2)*NLp + Z(1)*NL))./(1.0 + B(2)*del_t*k_til/Re);
     %**********************************************************************
    %************************THIRD SUBSTEP********************************
    %**********************************************************************
    [NLdp] = Convol(omeg_hat,Nx,Ny,K,L);
    omeg_hat = (omeg_hat + del_t*(-A(3)*(k_til.*omeg_hat)/Re + G(3)*NLdp + Z(2)*NLp))./(1.0 + B(3)*del_t*k_til/Re);
    j = j + 1;
    u_hat = 1i*L.*omeg_hat./k_til; u_hat(1, 1) = 0.0+1i*0.0;
    v_hat = -1i*K.*omeg_hat./k_til; v_hat(1, 1) = 0.0+1i*0.0;
    E = 0.5*sum(sum(u_hat.*conj(u_hat)+v_hat.*conj(v_hat)))/(Nx*Ny);
    %calculating total enstrophy
    Es = 0.5*(sum(sum(omeg_hat.*conj(omeg_hat))))/(Nx*Ny);
    fprintf(fp, '%d    %14.8e    %14.8e    %14.8e\n', j, time, E, Es);
    if(mod(j, 100)==0)
        omg = ifft2(omeg_hat);
        [hC, hC] = contourf(X, Y, real(omg), 50); axis equal;
        m = m + 1; 
	%surf(X, Y, real(omg), 50)
        set(hC,'LineStyle','none'); colorbar; colormap redblue(100); shading interp;drawnow
        saveas(hC,sprintf('FIG_ext%d.png',m));
    end
end
%%%%%%%%%%%%%%%END OF THE MAIN CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Convolution function for non linear term%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[udu] =  Convol(omeg_hat,Nx,Ny,K,L)
    Px = int32(Nx/2);Py = int32(Nx/2); 
    M1 = zeros(Ny/2,Px) + 1i*zeros(Ny/2,Px); %number of wavenumbers padded with zeros
    M2 = zeros(Py,Nx+Px) + 1i*zeros(Py,Nx+Px); %number of wavenumbers padded with zeros
    sx1 = 1; ex1 = Nx/2; sx2 = Nx/2+Px+1; ex2 = Nx+Px;
    sy1 = 1; ey1 = Ny/2; sy2 = Ny/2+Py+1; ey2 = Nx+Py;
    k_til = K.*K + L.*L;
    temp = (1i*K).*omeg_hat; 
    temp_extd = [temp(1:Ny/2, 1:Nx/2), M1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          M2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), M1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    temp = (1i*L).*omeg_hat./k_til; temp(1, 1) = 0.0+1i*0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), M1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          M2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), M1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = -omg_conv_extd.*psi_conv_extd;
    temp = (1i*L).*omeg_hat;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), M1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          M2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), M1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    omg_conv_extd = ifft2(temp_extd);
    
    temp = (1i*K).*omeg_hat./k_til; temp(1, 1) = 0.0+1i*0.0;
    temp_extd = [temp(1:Ny/2, 1:Nx/2), M1, temp(1:Ny/2, Nx/2+1:Nx);... 
                                          M2;...
                 temp(Ny/2+1:Ny, 1:Nx/2), M1, temp(Ny/2+1:Ny, Nx/2+1:Nx)];
    psi_conv_extd = ifft2(temp_extd);
    conv_extd = conv_extd + omg_conv_extd.*psi_conv_extd;
    conv_hat_extd = fft2(conv_extd);
    udu = [conv_hat_extd(sy1:ey1, sx1:ex1), conv_hat_extd(sy1:ey1, sx2:ex2); ...
           conv_hat_extd(sy2:ey2, sx1:ex1), conv_hat_extd(sy2:ey2, sx2:ex2)];
end
%%%%%%%%%%%%%%%%End of the Non linear term calculation%%%%%%%%%%%%%%%%%%
