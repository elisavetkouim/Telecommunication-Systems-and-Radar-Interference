epsilon = 10^(-0.82);
S_db = 30;
S = 10^(S_db/10);
I_db = 2 * S_db;
I = 10^(I_db/10);

M = 2;
N = 5e4;


while true
    nvars = 2 * M;
    x0 = -(M-1):2:(M-1);
    x0 = x0 / sqrt(mean(abs(x0).^2));
    z0 = [real(x0), imag(x0)];
    objfun = @(z) PeHighINR(z, S, I, N) + 0.1 * abs(mean(abs(z(1:M) + 1j*z(M+1:end)).^2) - 1);
    
    lb = -5 * ones(1, nvars);
    ub =  5 * ones(1, nvars);
        
    opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e3, 'MaxIterations', 1e3);
    

    problem = createOptimProblem('fmincon','objective', objfun,'x0', z0, 'lb', lb, 'ub', ub, 'options', opts);
        
    gs = GlobalSearch('NumStageOnePoints', 50, 'NumTrialPoints', 200, 'Display', 'off');
    [z_opt, Pe_val] = run(gs, problem);
    

    if Pe_val <= epsilon
        M_best = M;
        z_best = z_opt;
        M = M + 1;
    else
        break;
    end
end


x_opt = z_best(1:M_best) + 1j * z_best(M_best+1:end);
figure;
plot(real(x_opt), imag(x_opt), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
xlim([-2 2]);
ylim([-2 2]);
grid on; axis equal;
title('Sdb = 30, IdB = 2SdB & ε = 10^-0.82');
xlabel('Re'); ylabel('Im');