function Pe = PeLowINR(z, S, I)
    M = length(z)/2;
    x = z(1:M) + 1j*z(M+1:end);

    d = abs(x.' - x);   %MxM
    d(1:M+1:end) = inf;

    dmin = min(d(:));
    Nmin = sum((abs(d(:)-dmin) < 1e-9)) /M;

    Pe = Nmin*(1/(2*pi))* integral(@(theta) qfunc(sqrt(S*(dmin^2)/2) - sqrt(2*I)*cos(theta)), 0, 2*pi, 'AbsTol', 1e-10, 'RelTol', 1e-6);
end
