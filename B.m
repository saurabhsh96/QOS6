%Calculating B
function val = B(k0, kx, ky, w, l)
    Jt = sinc(-ky.*w/2/pi);
    It = (2.*k0).*(cos(-kx.*l./2) - cos(k0*l/2))...
            ./((k0.^2 - kx.^2).*sin(k0*l/2));
    val = Jt.*It; 
end