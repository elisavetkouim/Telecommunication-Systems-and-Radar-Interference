epsilon = 10^-5;
S_db = 20;
S = 10^(S_db/10);
I_db = 0.25*S_db;
I = 10^(I_db./10);

M = 2;

while true
    nvars = 2*M;
    x0 = randn(1, M) + 1j* randn(1, M);
    x0 = x0 / sqrt(mean(abs(x0).^2));
    z0 = [real(x0), imag(x0)];
    objfun = @(z) PeLowINR(z, S, I) + 0.005 * abs(mean(z(1:M).^2 + z(M+1:end).^2) - 1);

    lb = -5 * ones(1, nvars);
    ub =  5 * ones(1, nvars);
    
    opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e4);
    
    problem = createOptimProblem('fmincon', 'objective', objfun,'x0', z0, 'lb', lb, 'ub', ub, 'options', opts);

    gs = GlobalSearch("NumStageOnePoints", 50, 'NumTrialPoints', 100, 'Display', 'off');
    [z_opt, Pe_val] = run(gs, problem);


    if Pe_val <= epsilon
         M_best = M;
         z_best = z_opt;
         M = M + 1;
    else
        break;
    end
end


x_opt = z_best(1:M_best) + 1j*z_best(M_best+1:end);
figure;
plot(real(x_opt), imag(x_opt), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
grid on;
xlim([-2 2]);
ylim([-2 2]);
axis equal;
title('SdB = 20, IdB = 0.25SdB & ε = 10^-5' );