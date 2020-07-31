%AEP Calculation

%% Defining inputs

%Geometric parameters
dx = 20e-3;
dy = 20e-3;
w = 1e-3;
l = 15e-3;

%EM parameters
freq = 10e9;
c = 3e8;
lam = c/freq;
k0 = 2*pi/lam;
zeta = 377;

%TM Incidance
drad = pi/180;
dth = drad;
th = eps:dth:pi/2-dth;
phi = 0*drad; %% E-plane
%th0 = 0;
%ph0 = 0;

%Array parameters
range = 5;
mx = -range:1:range-1;
my = mx;

%Voltage
Vtm = 1;

%% Calculating incident field
einc = Vtm/sqrt(dx*dy);

%% Calculating IBF

kx = k0.*sin(th).*cos(phi);
ky = k0.*sin(th).*sin(phi);
kz = k0.*cos(th);

z = zeros(size(th));
iBF = zeros(size(th));
v = einc.*B(k0, -kx, -ky, w, l);
% v = einc.*B(k0, -kx(1), -ky(1), w, l);
% z = ZActive(k0,mx,my,th0,ph0,l,w,dx,dy);
% iBF = v./z;
for ind = 1:size(th, 2)
    z(ind) = ZActive(k0,mx,my,th(ind),phi,l,w,dx,dy);
    iBF(ind) = v(ind)./(50+z(ind));
end

%% Calculating AEP

bSGF = createSGF(k0, kx, ky, zeta, th);
r_temp = 1000.*lam;
exp_term = exp(-1j.*k0.*r_temp)./(2*pi*r_temp);
AEPx1 = 1j.*kz.*(squeeze(bSGF(1,1,1,:)))'.*iBF.*B(k0, kx, ky, w, l).*exp_term;
AEPx2 = 1j.*kz.*(squeeze(bSGF(2,1,1,:)))'.*iBF.*B(k0, kx, ky, w, l).*exp_term;
AEPx3 = 1j.*kz.*(squeeze(bSGF(3,1,1,:)))'.*iBF.*B(k0, kx, ky, w, l).*exp_term;
AEPm = sqrt((abs(AEPx1)).^2 + (abs(AEPx2)).^2 + (abs(AEPx3)).^2);
%AEPm = sqrt((abs(AEPx1)).^2);
AEP_max = max(AEPm);

%Question: Why pow2db?
plot([-th(1,size(th, 2):-1:1)./drad th(1,:)./drad], [mag2db(abs(AEPm(1,size(AEPm,2):-1:1)./AEP_max))...
    mag2db(abs(AEPm(1,:)./AEP_max))], 'LineWidth', 1.5);
ylim([-25, 0]);
xlabel('Scan Angle \theta_0 (in deg)');
ylabel('AEP (in dB) \phi = 0 [rad]');
title('Active Element Pattern, TE-Incidence');
grid on;

%% Total pattern

%Obs points0
dtho = drad/4;
dpho = drad/4;
thov = eps:dtho:pi/2-dtho;
phov = eps:dpho:2*pi-dpho;
[tho, pho] = meshgrid(thov, phov);

thin = [0 31]*drad;
phin = 0*drad;

[AFx, AFy, eFieldx1, eFieldx2, eFieldx3] = AFfield(tho, pho, ...
    thin(1), phin, iBF, k0, w, l, th, zeta, mx, my, dx, dy, exp_term);

[AFx1, AFy1, eFieldx11, eFieldx21, eFieldx31] = AFfield(tho, pho, ...
    thin(2), phin, iBF, k0, w, l, th, zeta, mx, my, dx, dy, exp_term);

% EFfth = (EFfx.*cos(tho).*cos(pho)) + (EFfy.*cos(tho).*sin(pho)) - (EFfz.*sin(tho));
% EFfphi = (-EFfx.*sin(pho)) + (EFfy.*cos(pho));
    
figure(2);
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(AFx(size(tho,1)/2+1,size(tho, 2):-1:1)./max(max(AFx(:,:))))) ...
    mag2db(abs(AFx(1,:)./max(max(AFx(:,:)))))], 'LineWidth', 1.5); hold on;
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(AFy(round(2*size(tho,1)/3+1),size(tho, 2):-1:1)./max(max(AFy)))) ...
    mag2db(abs(AFy(91,:)./max(max(AFy))))], '--', 'LineWidth', 1.5);
ylim([-40, 0]);
title('Array Factor vs. theta observation');
xlabel('\theta_{obs} (in deg)');
ylabel('AFx and AFy (in dB)');
legend('AFx \phi = 0', 'AFy \phi = 90');
grid on;

eMag = sqrt((abs(eFieldx1)).^2 + (abs(eFieldx2)).^2 + (abs(eFieldx3)).^2);
eMagMax = max(max(eMag(:,:)));

eMag1 = sqrt((abs(eFieldx11)).^2 + (abs(eFieldx21)).^2 + (abs(eFieldx31)).^2);
eMagMax1 = max(max(eMag1(:,:)));

figure(3);
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(eMag(size(tho,1)/2+1,size(tho, 2):-1:1)./eMagMax)) ...
     mag2db(abs(eMag(1,:)./eMagMax))], 'LineWidth', 1.5); hold on;
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(eMag1(size(tho,1)/2+1,size(tho, 2):-1:1)./eMagMax)) ...
     mag2db(abs(eMag1(1,:)./eMagMax))], 'LineWidth', 1.5); hold on;
% plot(tho(1,:)./drad, mag2db(abs(eMag(1,:)./eMagMax)), 'LineWidth', 1.5); hold on;
ylim([-40, 0]);
title('Normalized Electric Far-Field Pattern');
xlabel('\theta_{obs}(in deg)');
ylabel('Normalized Electric Field (in dB)');
legend('Scan Angle \theta_i = 0 deg', 'Scan Angle \theta_i = 30 deg');
grid on;

figure(4);
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(eMag(size(tho,1)/2+1,size(tho, 2):-1:1)./eMagMax)) ...
     mag2db(abs(eMag(1,:)./eMagMax))], 'LineWidth', 1.5); hold on;
plot([-tho(1, size(tho, 2):-1:1)./drad tho(1,:)./drad], [mag2db(abs(eMag1(size(tho,1)/2+1,size(tho, 2):-1:1)./eMagMax)) ...
     mag2db(abs(eMag1(1,:)./eMagMax))], 'LineWidth', 1.5); hold on;
plot([-th(1,size(th, 2):-1:1)./drad th(1,:)./drad], [mag2db(abs(AEPm(1,size(AEPm,2):-1:1)./AEP_max))...
    mag2db(abs(AEPm(1,:)./AEP_max))], 'LineWidth', 1.5);
% plot(tho(1,:)./drad, mag2db(abs(eMag(1,:)./eMagMax)), 'LineWidth', 1.5); hold on;
ylim([-40, 0]);
title('Normalized Far-Field and AEP');
xlabel('\theta_{obs}(in deg)');
ylabel('Normalized Electric Field (in dB)');
legend('Scan Angle \theta_i = 0 deg', 'Scan Angle \theta_i = 30 deg');
grid on;