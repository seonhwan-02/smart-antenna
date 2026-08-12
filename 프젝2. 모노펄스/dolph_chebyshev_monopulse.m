clear; clc; close all;

%% 1. 파라미터 설정 (System Parameters)
N = 21;                        % 소자 수
SLL_dB = 30;                   % 원하는 부엽 수준 (30 dB)
freq = 1e9;                    % 주파수 (1GHz)
c = 3e8;                       % 빛의 속도 (m/s)
lambda = c / freq;             % 파장
d = lambda / 2;                % 소자 간격 (반파장)
k_wave = 2 * pi / lambda;      % 파수 (Wave number)

% 10도 간격 고정 빔 스캔 각도 레이아웃
scan_angles = -60:10:60;       % -60, -50, -40, ..., 50, 60 (총 13개 빔)
theta = -90:0.1:90;            % 관측 각도
theta_rad = deg2rad(theta);    % 라디안 변환

%% 2. 체비쇼프 가중치 계산 (이전 방식 유지)
R = 10^(SLL_dB / 20);
z0 = cosh(acosh(R) / (N - 1));
W_k = zeros(1, N);

for k = 0:N-1
    u_k = (2 * pi * k) / N;
    x_k = z0 * cos(u_k / 2);
    m = N - 1;
    if abs(x_k) <= 1
        W_k(k+1) = cos(m * acos(x_k));
    else
        W_k(k+1) = cosh(m * acosh(x_k));
    end
end

w_cheb = zeros(1, N);
for n = 1:N
    sum_val = 0;
    for k = 0:N-1
        sum_val = sum_val + W_k(k+1) * exp(1j * 2 * pi * k * (n - (N+1)/2) / N);
    end
    w_cheb(n) = real(sum_val);
end
w_cheb = w_cheb / max(w_cheb);

element_pos = (1:N) - (N+1)/2;

%% 3. 각 스캔 각도별 독립적인 체비쇼프 빔 패턴(Array Factor) 선행 계산
num_beams = length(scan_angles);
AF_matrix = zeros(num_beams, length(theta));

for idx = 1:num_beams
    theta0 = scan_angles(idx);
    theta0_rad = deg2rad(theta0);

    % 위상 조향 가중치 계산
    steering_phase = exp(-1j * k_wave * d * element_pos * sin(theta0_rad));
    w_sum = w_cheb .* steering_phase;

    % 공각 공간 선형 빔 수신 계산
    for t = 1:length(theta)
        psi = k_wave * d * sin(theta_rad(t));
        spatial_phase = exp(1j * element_pos * psi);
        AF_matrix(idx, t) = abs(sum(w_sum .* spatial_phase)); % 전압 크기 성분 추출
    end
    % 정규화
    AF_matrix(idx, :) = AF_matrix(idx, :) / max(AF_matrix(idx, :));
end

%% 4. 메인 디스플레이 창 생성 (좌측: 빔패턴 그물, 우측: 인접구간 교차 모노펄스)
figure('Color', 'w', 'Position', [100, 100, 1400, 600]);

% 왼쪽 : Beam Pattern 그물망 플롯]
subplot(1, 2, 1);
for idx = 1:num_beams
    AF_dB = 20 * log10(AF_matrix(idx, :));
    AF_dB(AF_dB < -50) = -50;
    plot(theta, AF_dB, 'b', 'LineWidth', 1.2); hold on;
    xline(scan_angles(idx), 'k:', 'LineWidth', 0.5);
end
grid on;
yline(-SLL_dB, 'r--', 'LineWidth', 1.5);
xlabel('Angle (Degrees)', 'FontSize', 12);
ylabel('Normalized Magnitude (dB)', 'FontSize', 12);
title('Dolph-Chebyshev Beam Lattice (10\circ Steps)', 'FontSize', 14);
xlim([-90 90]); ylim([-50 5]);

% 오른쪽 : 인접 빔 조합을 이용한 구간별 모노펄스 특성
subplot(1, 2, 2);

for idx = 1:(num_beams - 1)
    beam1_center = scan_angles(idx);        % ex: -60도
    beam2_center = scan_angles(idx + 1);     % ex: -50도
    mid_angle = (beam1_center + beam2_center) / 2; % 구간의 중간값 (ex: -55도)

    % 인접 두 빔의 데이터 추출
    AF_beam1 = AF_matrix(idx, :);
    AF_beam2 = AF_matrix(idx + 1, :);

    % 수신단 인접 빔 조합 모노펄스 연산 (차 패턴 및 합 패턴)
    Sum_pattern = AF_beam2 + AF_beam1;
    Diff_pattern = AF_beam2 - AF_beam1;
    mono_ratio = Diff_pattern ./ Sum_pattern;

    % 두 빔의 교차 영역인 인접 스캔 경계구간 내부만 1차함수 형태
    window_mask = (theta >= beam1_center) & (theta <= beam2_center);
    mono_ratio(~window_mask) = NaN;

    % 구간 오차 사선 플롯(0되는 지점)
    plot(theta, mono_ratio, 'LineWidth', 2.5); hold on;

    % 중간값 영점 위치에 수직 점선 표시 확인용
    xline(mid_angle, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.5);
end

grid on;
yline(0, 'k-', 'LineWidth', 1); % Y=0 기준선
xlabel('Angle (Degrees)', 'FontSize', 12);
ylabel('Monopulse Ratio (AF_2 - AF_1) / (AF_2 + AF_1)', 'FontSize', 12);
title('Adjacent Beam Pair Monopulse curves (Centered at Mid-Angles)', 'FontSize', 14);
xlim([-70 70]); ylim([-1 1]);
