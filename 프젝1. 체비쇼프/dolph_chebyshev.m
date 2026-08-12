N = input('배열 소자 수(N)를 입력하세요: ');
SLL_dB = input('원하는 부엽 수준(SLL, 양수 dB 단위)을 입력하세요: ');

freq = 10e9;
c = 3e8;
lambda = c / freq;
d = lambda / 2;

R = 10^(SLL_dB / 20);
z0 = cosh(acosh(R) / (N - 1));

W_k = zeros(1, N);
for k = 0:N-1
    u_k = (2 * pi * k) / N;
    x_k = z0 * cos(u_k / 2);
    m = N - 1;
    if abs(x_k) <= 1
        W_k(k+1) = cos(m * acosh(x_k));
    else
        W_k(k+1) = cosh(m * acosh(x_k));
    end
end

w = zeros(1, N);
for n = 1:N
    sum_val = 0;
    for k = 0:N-1
        sum_val = sum_val + W_k(k+1) * exp(1j * 2 * pi * k * (n - (N+1)/2) / N);
    end
    w(n) = real(sum_val);
end

w = w / max(w);

theta = -90:0.1:90;
AF = zeros(1, length(theta));

for n = 1:N
    element_pos = n - (N + 1) / 2;
    psi = (2 * pi * d / lambda) * sind(theta);
    AF = AF + w(n) * exp(1j * element_pos * psi);
end

AF_mag = 20 * log10(abs(AF) / max(abs(AF)));
AF_mag(AF_mag < -SLL_dB - 20) = -SLL_dB - 20;
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
plot(theta, AF_mag, 'b', 'LineWidth', 2);
grid on; hold on;
yline(-SLL_dB, 'r--', 'LineWidth', 2);
title(['[Dolph-Chebyshev] N=', num2str(N), ', SLL=-', num2str(SLL_dB), 'dB'], 'FontSize', 13);
xlabel('Angle (Degrees)', 'FontSize', 11);
ylabel('Normalized Magnitude (dB)', 'FontSize', 11);
xlim([-90 90]);
ylim([-SLL_dB-20 5]);
legend('Beam Pattern (Equal-Ripple)', 'Target SLL Line');
