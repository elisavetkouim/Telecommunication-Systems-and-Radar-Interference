M = 16;
Sdb = 20;
S = 10^(Sdb/10);
Idb = 0.25*Sdb;
I = 10^(Idb/10);

nvars = 2*M;

i = (1:M);
[re, im] = meshgrid([-3 -1 1 3], [-3 -1 1 3]);
x0 = (re(:) + 1j*im(:))';
x0 = x0 / sqrt(mean(abs(x0).^2));
z0 = [real(x0), imag(x0)];

objfun = @(z) PeLowINR(z, S, I) + 10 * abs(mean(z(1:M).^2 + z(M+1:end).^2) - 1);
lb = -5 * ones(1, nvars);
ub =  5 * ones(1, nvars);

opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6);

problem = createOptimProblem('fmincon', 'objective', objfun, 'x0', z0, 'lb', lb, 'ub', ub, 'options', opts);

gs = GlobalSearch("NumStageOnePoints", 100, 'NumTrialPoints', 1000, 'Display', 'off');
[z_opt, Pe_val] = run(gs, problem);

x_opt = z_opt(1:M) +1j* z_opt(M+1:end);

figure;
plot(real(x_opt), imag(x_opt), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
xlim([-2 2]);
ylim([-2 2]);
grid on; axis equal;
title('SdB = 20, IdB = 0.25SdB M=16');
xlabel('Re'); ylabel('Im');