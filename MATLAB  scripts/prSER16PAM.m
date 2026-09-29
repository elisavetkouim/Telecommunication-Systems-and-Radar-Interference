S_db = 10;
S = 10^(S_db/10);
I_db = -20:2:60;
I = 10.^(I_db./10);

M = 16;

SER_plot_lowerBound = 2*(1-1/M)*(1/(2*pi))*qfunc(sqrt(6*S/(M^2 - 1))) *ones(size(I));

for idx = 1:length(I)
    SER = 2*(1-1/M)*(1/(2*pi))*integral(@(theta) qfunc(sqrt(6*S/(M^2 - 1)) - sqrt(2*0)*cos(theta)), 0, 2*pi);
    SER_plot_lowerBound(idx) = SER;
end

SER_plot_ML = zeros(size(I));
SER_plot_IC = zeros(size(I));
N = 5e5;
x = -(M-1):2:(M-1);
x = x / sqrt(mean(x.^2));

X_idx = randi([1 M] ,1, N);
X = x(X_idx);
Theta = 2*pi*rand(1,N);
Z = wgn(N,1,0) + 1j*wgn(N,1,0);
for idx = 1:length(I)
    Y = sqrt(S)*X + sqrt(I(idx))*exp(1j*Theta) + Z'/sqrt(2);
    lminML = zeros(M,N);
    lminIC = zeros(M,N);
    for m = 1:M
        d = abs(Y-sqrt(S)*x(m));
        logI0 = log(besseli(0, 2*sqrt(I(idx))*d, 1)) + 2*sqrt(I(idx))*d;
        lminML(m,:) = d.^2 - logI0;
        lminIC(m,:) = (abs(Y-sqrt(S)*x(m))-sqrt(I(idx))).^2;
    end
    [~,X_r_idxML] = min(lminML, [], 1);
    SER_plot_ML(idx) = mean(X_r_idxML ~= X_idx);
    [~,X_r_idxIC] = min(lminIC, [], 1);
    SER_plot_IC(idx) = mean(X_r_idxIC ~= X_idx);
end


SER_plot_TIN = zeros(size(I));
for idx = 1:length(I)
    SER = 2*(1-1/M)*(1/(2*pi))*integral(@(theta) qfunc(sqrt(6*S/(M^2 - 1)) - sqrt(2*I(idx))*cos(theta)), 0, 2*pi);
    SER_plot_TIN(idx) = SER;
end


SER_plot_asymptote = zeros(size(I));
for idx = 1:length(I)
    SER = 2*(1-(1/M))*(1/(2*pi))*integral(@(theta) qfunc(sqrt(6*S*(cos(theta)).^2/(M^2-1))), 0, 2*pi);
    SER_plot_asymptote(idx) = SER;
end




figure;
plot(I_db/S_db, log10(SER_plot_lowerBound), '-.');   hold on;
plot(I_db/S_db, log10(SER_plot_ML));   hold on;
plot(I_db/S_db, log10(SER_plot_IC), '--');   hold on;
plot(I_db/S_db, log10(SER_plot_TIN),'-.');   hold on;
plot(I_db/S_db, log10(SER_plot_asymptote), '-.');
xlabel('SNR (dB)');
ylabel('Symbol Error Rate');
legend('lower bound INR=0', 'ML decoder (optimal)' ,'IC decoder (suboptimal)', 'TIN decoder (suboptimal)', 'asymptote (IC) INR>>SNR>>1');
title('16-PAM: SER vs I in dB, for fixed S.');
grid on;
