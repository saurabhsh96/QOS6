function [AFx, AFy, eFieldx1, eFieldx2, eFieldx3] = AFfield(tho, pho, ...
    thin, phin, iBF, k0, w, l, th, zeta, mx, my, dx, dy, exp_term)
    kxo = k0.*sin(tho).*cos(pho);
    kyo = k0.*sin(tho).*sin(pho);
    kzo = k0.*cos(tho);

    kxin = k0.*sin(thin).*cos(phin);
    kyin = k0.*sin(thin).*sin(phin);

    AFx = zeros(size(tho));
    AFy = zeros(size(tho));

    %AF calculations, Nx = Ny
    %for indj = 1:size(thin, 2)
    for ind = 0:size(mx, 2)-1
        AFx = AFx + exp(1j.*(kxo - kxin).*ind.*dx);
        AFy = AFy + exp(1j.*(kyo - kyin).*ind.*dy);
    end
    %end
    %Jw calculation
    %iBFReq = zeros(size(thin));
    B_term = B(k0, kxo, kyo, w, l);
    iBFReq = iBF(find(ismembertol(th,thin)));
    Jw = iBFReq.*B_term.*AFx.*AFy;
    %Jw = ;
    %use exp_term from last section
    bSGF = createSGF(k0, kxo, kyo, zeta, tho(1));
    eFieldx1 = 1j.*kzo.*squeeze(bSGF(1,1,:,:)).*Jw.*exp_term;
    eFieldx2 = 1j.*kzo.*squeeze(bSGF(2,1,:,:)).*Jw.*exp_term;
    eFieldx3 = 1j.*kzo.*squeeze(bSGF(3,1,:,:)).*Jw.*exp_term;
end