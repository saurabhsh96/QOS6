%Finding gamma
% Antenna array calculations
%% Defining inputs

%Distances in X and Y direction
dx = 15e-3;
dy = 15e-3;

%Dimensions of single dipole
w = 1e-3;
l = 14e-3;

%Defining the mesh
drad = pi/180;
dth = drad;
dph = drad;
%[th, ph] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);
th = 0;
ph = 0;

%Defining mx and my indexes
%Bounds
upper = 10;
lower = -10;
mx = lower:1:upper;
my = mx;

%Voltage
Vtm = 1;

%Defining frequency range
freq = 6e9:100e6:14e9;

c = 3e8; %Speed of causality :P

indZ = 1;
%Defining Z
z = zeros(1, size(freq, 2));
v = zeros(1, size(freq, 2));
escatt = zeros(1, size(freq, 2));
iBF = zeros(1, size(freq, 2));
zeta = 377;

%% Calculating einc
einc = Vtm/sqrt(dx*dy);

%% Looping over Theta cases and Frequency
for ind = 1:size(th, 2)
    for f = freq
        lam = c/f;
        k0 = 2*pi/lam;
        z(ind, indZ) = ZActive(k0,mx,my,th(ind),ph,l,w,dx,dy);
        kx = k0.*sin(th).*cos(ph);
        ky = k0.*sin(th).*sin(ph);
        
        Jt = sinc(-ky.*w/2/pi);
        It = (2.*k0).*(cos(-kx.*l./2) - cos(k0*l/2))...
            ./((k0.^2 + kx.^2).*sin(k0*l/2));

        v(ind, indZ) = einc.*Jt.*It;
        iBF(ind, indZ) = v(ind, indZ)./z(ind, indZ);

        Jt = sinc(ky.*w/2/pi);
        It = (2.*k0).*(cos(kx.*l./2) - cos(k0*l/2))...
            ./((k0.^2 - kx.^2).*sin(k0*l/2));

        bSGF = createSGF(k0, kx, ky, zeta, th);
        
        escatt(ind, indZ) = (iBF(ind, indZ)./(dx*dy)).*bSGF(1,1).*Jt.*It;
        indZ = indZ + 1;
    end
    indZ = 1;
end

gamma = escatt./einc;
T = 1+gamma;

%Plotting
figure(1);
plot(freq, abs(gamma), 'LineWidth', 1.5); hold on;
plot(freq, abs(T), 'LineWidth', 1.5);
title('Gamma Plot');
xlabel('Frequency in Hz');
ylabel('Gamma');
legend('gamma', 'T');

% figure(2);
% plot(freq, real(z(2,:)), 'LineWidth', 1.5); hold on;
% plot(freq, imag(z(2,:)), 'LineWidth', 1.5);
% title('Real and Imaginary Parts of Active Impedence \Theta = 30');
% xlabel('Frequency in Hz');
% ylabel('Real(Z) and Imag(Z) [in Ohm]');
% legend('Real', 'Imag');

%% Case 2
freq = 10e9;
[th, ph] = meshgrid(eps:drad:pi/2, eps:pi/2:pi/2);

z = zeros(size(th));
v = zeros(size(th));
escatt = zeros(size(th));
iBF = zeros(size(th));

lam = c/f;
k0 = 2*pi/lam;
        
kx = k0.*sin(th).*cos(ph);
ky = k0.*sin(th).*sin(ph);
        
for ind = 1:size(ph, 1)
    for indt = 1:size(th, 2)
        z(ind, indt) = ZActive(k0,mx,my,th(ind, indt),ph(ind, indt),l,w,dx,dy);
        
        Jt = sinc(-ky(ind, indt).*w/2/pi);
        It = (2.*k0).*(cos(-kx(ind, indt).*l./2) - cos(k0*l/2))...
            ./((k0.^2 + kx(ind, indt).^2).*sin(k0*l/2));

        v(ind, indt) = einc.*Jt.*It;
        iBF(ind, indt) = v(ind, indt)./z(ind, indt);

        Jt = sinc(ky(ind, indt).*w/2/pi);
        It = (2.*k0).*(cos(kx(ind, indt).*l./2) - cos(k0*l/2))...
            ./((k0.^2 - kx(ind, indt).^2).*sin(k0*l/2));

        bSGF = createSGF(k0, kx, ky, zeta, th(ind, indt));
        
        escatt(ind, indt) = (iBF(ind, indt)./(dx*dy)).*bSGF(1,1).*Jt.*It;
    end
    %indZ = 1;
end

gamma = escatt./einc;
T = 1+gamma;

%Plotting
figure(2);
plot(th(1,:)./drad, abs(gamma(1,:)), 'LineWidth', 1.5); hold on;
plot(th(1,:)./drad, abs(gamma(2,:)), 'LineWidth', 1.5);
title('Gamma Plot');
xlabel('Frequency in Hz');
ylabel('Gamma');
legend('Phi = 0', 'Phi = 90');

