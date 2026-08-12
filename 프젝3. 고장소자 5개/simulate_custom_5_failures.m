% =========================================================================
% 사용자 맞춤형 5개 고장 소자 시뮬레이터 (LS 보상 전/후 비교)
% =========================================================================
clear; clc; close all;

%% 1. 사용자 설정 (고장 소자 무작위 추출)
% 1~21 사이의 숫자 5개를 무작위로 추출합니다.
fault_idx = sort(randperm(21, 5)); % 예: [2 5 13 18 20]

%% 2. 기본 안테나 배열 및 체비쇼프 설계 파라미터
N = 21;             % 총 소자 수
freq = 3e9;         % 주파수: 3 GHz
c0 = 3e8;           % 빛의 속도
lambda = c0 / freq; % 파장
d = lambda / 2;     % 소자 간격 (반파장)
k = 2 * pi / lambda;% 파수(Wave number)
SLL = 25;           % 목표 부엽 수준 (dB)

% Chebyshev 타겟 패턴 유도
R = 10^(SLL / 20);
beta = cosh(acosh(R) / (N - 1));
sample_num = 1801;  % 연산 속도를 위해 1801 샘플 사용 (-90~90도 0.1도 간격)
psi = linspace(0, pi, 10000).';

target = zeros(10000, 1);
for i = 1:10000
    x = beta * cos(psi(i) / 2);
    if abs(x) <= 1
        target(i) = cos((N - 1) * acos(x));
    else
        target(i) = cosh((N - 1) * acosh(x));
    end
end
target = target / cosh((N - 1) * acosh(beta));

% 가중치 계산 (행렬 방정식)
M = (N - 1) / 2;
A_cheb = zeros(10000, M + 1);
A_cheb(:, 1) = 1;
for m = 1:M
    A_cheb(:, m + 1) = 2 * cos(m * psi);
end
coeff = A_cheb \ target;
w = [fliplr(coeff(2:end).'), coeff(1), coeff(2:end).'];
w = real(w);
w = w / sum(w); % 정규화

%% 3. 조향 행렬 (Steering Matrix) 생성
theta_deg = -90:0.1:90;
theta_rad = theta_deg * pi / 180;
element_position = ((0:N-1) - (N-1)/2) * d;

A_full = exp(1j * k * element_position.' * sin(theta_rad)); % 21 x 1801
b_desired = w * A_full; % 1 x 1801
BP_D_mag = abs(b_desired) / max(abs(b_desired));
BP_D_dB = 20 * log10(BP_D_mag + 1e-12);
side_lobe_idx = abs(theta_deg) > 10; % 주엽(Main lobe)을 제외한 부엽 구간 인덱스

%% 4. 고장 처리 및 LS 보상 연산
remain_idx = setdiff(1:N, fault_idx); % 생존 소자 인덱스
A_fault = A_full(remain_idx, :);      % 생존 소자만의 조향 행렬 (16 x 1801)

% [보상 전: Fault (No LS)] 기존 가중치를 그대로 사용
w_fault = w(remain_idx); 
b_fault = w_fault * A_fault;
BP_fault_mag = abs(b_fault) / max(abs(b_fault));
BP_fault_dB = 20 * log10(BP_fault_mag + 1e-12);

% [보상 후: Hermitian LS Compensated] 가중치 재계산
A_fault_T = A_fault.'; % 1801 x 16 (전치)
b_desired_T = b_desired.';
w_LS = (A_fault_T' * A_fault_T) \ (A_fault_T' * b_desired_T); % 최소제곱법 최적해 도출
b_LS = (A_fault_T * w_LS).';
BP_LS_mag = abs(b_LS) / max(abs(b_LS));
BP_LS_dB = 20 * log10(BP_LS_mag + 1e-12);

%% 5. 성능 지표(SLL, MSE) 산출 및 출력
sll_normal = max(BP_D_dB(side_lobe_idx));
sll_fault = max(BP_fault_dB(side_lobe_idx));
sll_LS = max(BP_LS_dB(side_lobe_idx));

mse_fault = mean((BP_D_mag - BP_fault_mag).^2);
mse_LS = mean((BP_D_mag - BP_LS_mag).^2);

fprintf('================================================\n');
fprintf('고장 소자 위치: %s\n', mat2str(fault_idx));
fprintf('================================================\n');
fprintf('[목표 성능] 정상 SLL: %.2f dB\n\n', sll_normal);
fprintf('[보상 전(No LS)] SLL: %.2f dB (증가: +%.2f dB) | MSE: %.4f\n', sll_fault, sll_fault - sll_normal, mse_fault);
fprintf('[보상 후(LS 적용)] SLL: %.2f dB (남은 오차: +%.2f dB) | MSE: %.4f\n', sll_LS, sll_LS - sll_normal, mse_LS);
fprintf('👉 SLL 개선 폭: %.2f dB 복구\n', sll_fault - sll_LS);
fprintf('================================================\n');

%% 6. 결과 시각화 플롯
figure('Color', 'w', 'Position', [100, 100, 900, 500]);
plot(theta_deg, BP_D_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, BP_fault_dB, 'r--', 'LineWidth', 1.5);
plot(theta_deg, BP_LS_dB, 'g-', 'LineWidth', 1.5);

yline(-SLL, 'k--', 'Target SLL (-25dB)', 'LineWidth', 1);
ylim([-80 0]); xlim([-90 90]);
grid on;
xlabel('Angle \theta (deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Array Factor [dB]', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('고장 복원 시뮬레이션 (제거 소자: %s)', mat2str(fault_idx)), 'FontSize', 14, 'FontWeight', 'bold');
legend('Desired Pattern', 'Fault (No LS)', 'LS Compensated', 'Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11);
