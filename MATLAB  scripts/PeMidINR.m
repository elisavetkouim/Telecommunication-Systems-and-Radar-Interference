function Pe = PeMidINR(z, S, I, N)
    M = length(z)/2;
    x = z(1:M) + 1j*z(M+1:end);
    X_idx = randi([1 M] ,1, N);
    X = x(X_idx);
    Theta = 2*pi*rand(1,N);
    Z = (wgn(N,1,0) + 1j*wgn(N,1,0))' /sqrt(2);
    Y = sqrt(S)*X + sqrt(I)*exp(1j*Theta) + Z;
    lmin = zeros(M,N);
    for m = 1:M
        d = abs(Y-sqrt(S)*x(m));
        logI0 = log(besseli(0, 2*sqrt(I)*d, 1)) + 2*sqrt(I)*d;
        lmin(m,:) = d.^2 - logI0;
    end
    [~,X_r_idx] = min(lmin, [], 1);
    Pe = mean(X_r_idx ~= X_idx);
end