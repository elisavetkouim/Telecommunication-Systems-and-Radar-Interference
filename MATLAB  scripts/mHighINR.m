M = 16;
Sdb = 20;
S = 10^(Sdb/10);
Idb = 2*Sdb;
I = 10^(Idb/10);

nvars = 2*M;
N = 1e5;

x0 = randn(1, M) + 1j* randn(1, M);
x0 = x0 / sqrt(mean(abs(x0).^2));
z0 = [real(x0), imag(x0)];

objfun = @(z) PeHighINR(z, S, I, N) + 1000 * abs(mean(z(1:M).^2 + z(M+1:end).^2) - 1);
lb = -5 * ones(1, nvars);
ub =  5 * ones(1, nvars);

opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e8, 'MaxIterations', 1e8);

problem = createOptimProblem('fmincon', 'objective', objfun, 'x0', z0, 'lb', lb, 'ub', ub, 'options', opts);

gs = GlobalSearch("NumStageOnePoints", 500, 'NumTrialPoints', 10000, 'Display', 'off');
[z_opt, Pe_val] = run(gs, problem);

x_opt = z_opt(1:M) +1j* z_opt(M+1:end);

figure;
plot(real(x_opt), imag(x_opt), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
xlim([-2 2]);
ylim([-2 2]);
grid on; axis equal;
title('SdB = 20, IdB = 2*SdB M=16');
xlabel('Re'); ylabel('Im');