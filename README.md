[README.md](https://github.com/user-attachments/files/30980160/README.md)
# Dolph-Chebyshev 기반 스마트 안테나 빔포밍 설계

배열 이론 정립 → 인접 빔 모노펄스 추적 → 고장 소자 Least Square 보상 → 비균일 배열(NLA) 통계적 복원력 분석

김선환 · 광운대학교 전자융합공학과 · 스마트안테나 (4학년 1학기)

## 프로젝트 개요

네 개 프로젝트는 독립된 실습이 아니라 하나의 설계 축을 따라 점진적으로 확장된 연속 프로젝트다.

| Part | 제목 | 핵심 내용 | 코드 |
|---|---|---|---|
| I | Dolph-Chebyshev Beamforming 설계 | DFS 샘플링으로 가중치를 직접 도출해 단일 빔의 주엽-부엽 trade-off 제어 | `프젝1. 체비쇼프/dolph_chebyshev.m` |
| II | 인접 빔 페어 모노펄스 | Part I 가중치를 다중 조향 빔으로 확장, 인접 빔 합·차 패턴으로 표적 편이 각도 추적 | `프젝2. 모노펄스/dolph_chebyshev_monopulse.m` |
| III | 고장 소자 제거 및 Least Square 보상  | 소자 고장 시 Least Square로 빔 패턴 복원 가능 조건을 정량 검증 | `프젝3. 고장소자 5개/simulate_custom_5_failures.m` |
| IV | 비균일 배열(NLA)에서의 LS 복원력 통계 분석 | 소자 간격 오류 개수(n=1~14)에 따른 열화·복원 추이를 몬테카를로로 정량화  | `프젝4. NULA/ls_synthesis_1.py`, `ls_synthesis_5.py`, `ls_synthesis_n.py` |


**구현 환경**: MATLAB (Part I~III, Signal Processing Toolbox 미사용) / Python·NumPy (Part IV)

이론 설계 → 시스템 확장 → 고장 대응 → 통계적 일반화로 이어지는 하나의 엔지니어링 흐름이다.

---

## Part I. Dolph-Chebyshev Beamforming 기초 설계

안테나 배열 설계의 주엽 폭-부엽 수준(SLL) trade-off를 Dolph-Chebyshev 방식으로 해결한다. 체비쇼프 다항식을 이산 푸리에 샘플링(DFS)으로 가중치화하고 IDFT로 물리적 가중치를 도출했다.

- **설계 파이프라인**: SLL → 전압비 R → 임계값 z0 → DFS 샘플링(N개 지점) → IDFT → 가중치 w(n)
- **구현**: `dolph_chebyshev.m` — 소자 수 N, SLL을 입력받아 Array Factor를 계산·플로팅. 중심 주파수 10GHz, 소자 간격 λ/2.
- **검증 조건**: N = 10, 100, 1000 (SLL 고정 30dB) / SLL = 30 → 10dB (N 고정)
- **결과**: 설계 목표 SLL을 전 각도에서 오차 없이 만족하는 Equal-Ripple 패턴 구현 (`체비쇼프 n다르게.png`, `체비쇼프 sll다르게.png`)

## Part II. Chebyshev 기반 인접 빔 페어 모노펄스

단일 탐색 빔으로는 표적이 빔 중심/가장자리 중 어디에 있는지 구별할 수 없다는 문제를, 인접 빔 쌍의 합·차 패턴(모노펄스 비, S-curve)으로 해결했다.

- **구현**: `dolph_chebyshev_monopulse.m` — N=21, SLL 30dB, 중심 주파수 1GHz 기준으로 -60°~+60°를 10° 간격으로 스캔하는 13개 조향 빔 격자와 12개 인접 빔 페어의 모노펄스 곡선을 동시에 계산
- **모노펄스 판별식**: M = (AF_idx+1 − AF_idx) / (AF_idx+1 + AF_idx), 경계 구간 외 NaN 처리로 발산 방지
- **결과**: 전 구간에서 발산 없는 완벽한 선형 Zero-Crossing 곡선 확인 (`result_01_baseline_N21_SLL30_scan10.png`, `result_06_N21_SLL30_scan10.png`)

| 파라미터 변화 실험 | 결과 이미지 |
|---|---|
| N = 10 / 11 / 21 (SLL 30dB 고정) | `result_02_N10_SLL30.png`, `result_03_N11_SLL30.png`, `result_04_N21_SLL30.png` |
| SLL = 10 / 50dB (N=21 고정) | `result_05_N21_SLL10.png`, `result_07_N21_SLL50.png` |
| 스캔 간격 = 5° / 20° | `result_08_scan5.png`, `result_09_scan20.png` |

| 파라미터 | 관찰된 Trade-off |
|---|---|
| 소자 수 N (10/11/21) | N↑ → 주엽 예리해지나 선형 추적 가능 폭은 좁아짐 |
| SLL (10/30/50dB) | 10dB: 비선형 왜곡 발생 / 50dB: 감도 저하 / 30dB가 최적 절충점 |
| 스캔 간격 (5/10/20°) | 5°: 신뢰도↑·연산부하↑ / 20°: 연산부하↓·탐지 사각지대 위험↑ |

## Part III. 고장 소자 제거 및 Least Square 보상

N=21 배열(중심 주파수 3GHz, 목표 SLL -25dB)에서 임의 5개 소자 고장을 모사하고 Least Square로 가중치를 복원했다. 복소 신호 공간에서 단순 Transpose와 Hermitian Transpose의 차이를 이론·실험적으로 비교했다.

- **구현**: `simulate_custom_5_failures.m` — `randperm(21,5)`로 고장 소자를 무작위 추출하고, 정상(BP_D)/무보상(BP_fault)/Hermitian LS 복원(BP_C) 세 패턴을 SLL·MSE와 함께 비교 출력
- **복원식**: w_LS = (A_fault^H A_fault)^-1 A_fault^H b_D (H = Conjugate Transpose)
- **핵심 발견**: 단순 Transpose는 복소 기하 거리 정보를 소실해 복원이 붕괴(SLL -3~0dB)되는 반면, Hermitian Transpose는 실용적 수준(SLL -19.80dB)으로 복원됨

| Case | 제거 소자 위치 | SLL | 비고 | 결과 이미지 |
|---|---|---|---|---|
| 1a/1b | 끝단 대칭/비대칭 | -24.70dB | 양호한 복원 | `case_1_original.png` |
| 2a/2b | 중앙 연속/분산 | -11.45~-13.80dB | 복원 실패 | `case_2_original.png` |
| 3a/3b | 좌/우측 편중 | -24.70dB | LS가 위상중심 자동 보정 | `case_3_original.png` |
| 4a/4b | 분산(중심 생존/유실) | ≈-12dB / ≈-6dB | 중심 생존 여부가 결정적 | `case_4_original.png` |

**결론**: 복원 가능성은 결함 개수·분산도보다 "최고 가중치를 지닌 위상 중심 소자의 생존 여부"가 지배적이다. 외곽 소자 결함은 유효 개구면 축소만으로 복원 가능하나, 중심 소자 결함은 원천적으로 복원 불가능하다.

## Part IV. 비균일 배열(NLA)에서의 LS 복원력 통계 분석

기말고사 문제(N=15, SLL -20dB 배열에서 소자 간격 0.2λ~0.7λ 무작위 비균일화, 학번 기반 조향각 15° 적용)에서 다룬 도구를 그대로 가져오되, **몬테카를로 반복 실험으로 통계적 신뢰도까지 확장**한 별도 분석이다. 기말고사 제출본(`스마트_안테나_기말고사_제출.docx`)은 단일 시행·특정 조향각 기준이었던 반면, 이 코드는 반복 횟수와 오류 소자 개수를 체계적으로 변화시켜 LS의 복원력 자체를 정량화했다.

- **`ls_synthesis_1.py`**: 소자 1개 간격만 무작위(0.2λ~0.7λ)로 흔들었을 때, LS 미적용/적용 시 부엽(PSLL, 목표 -20dB)을 1000회 몬테카를로로 비교 (`ls_pattern_comparison_v2.png`, `ls_boxplot_v2.png`)
- **`ls_synthesis_5.py`**: 5개 소자 간격을 동시에 무작위화한 경우로 확장, 동일하게 1000회 반복 (`ls_pattern_comparison_5_v2.png`, `ls_boxplot_5_v2.png`)
- **`ls_synthesis_n.py`**: 오류 간격 개수 n을 1개부터 14개(전체)까지 증가시키며 각 n마다 500회씩 반복(`LS_Trend_Analysis.md`, `ls_trend_n_v2.png`). 목표 SLL을 -25dB로 더 엄격하게 설정

**핵심 결과** (`LS_Trend_Analysis.md` 기준, N=15, -25dB 타겟, n=14 전 구간 붕괴 조건):
- 기존 가중치 그대로 사용 시(No LS): 부엽이 평균 **-14.7dB**까지 열화되어 목표 달성 완전 실패
- LS로 가중치 재계산 시: n=14의 극한 조건에서도 부엽을 평균 **-21.6dB ~ -22.2dB**로 방어
- 두 방식의 격차는 n=14에서 약 **7dB**로, 손상이 심할수록 LS 최적화의 가치가 커짐을 확인

**원인 분석** (보고서 기재, 정성적 추정):
1. 그레이팅 로브 가설 — 간격이 λ/2를 넘어 최대 0.7λ까지 허용되며 불필요한 에너지가 가시영역에 유입
2. 비대칭 기하구조 가설 — 비균일·비대칭 배열에서 LS 오차 최소화 과정 중 에너지가 사이드 방향으로 분산

---

## 사용된 이론 및 기술 요소

| 구분 | 내용 |
|---|---|
| 배열 안테나 이론 | Dolph-Chebyshev Beamforming, Array Factor, SLL-HPBW Trade-off, 개구면(Effective Aperture), 테이퍼링, 그레이팅 로브 |
| 신호처리 수학 | Chebyshev 다항식, DFS(공간 주파수 샘플링), IDFT 기반 가중치 도출, Hermitian/단순 Transpose 비교 |
| 레이더 추적 이론 | 인접 빔 페어 모노펄스, 합·차 패턴, 모노펄스 비(S-curve), Zero-Crossing 판별 |
| 최적화 및 검증 기법 | Least-Squares, 투영법(Projection Method), MSE 기반 오차 평가, 몬테카를로 시뮬레이션(500/1000회), Boxplot 기반 통계 검증 |
| 구현 환경 | MATLAB(Part I~III, Signal Processing Toolbox 미사용) / Python·NumPy·Matplotlib(Part IV) |

## 파일 구조

```
smart-antenna-beamforming/
├── README.md
├── 프젝1. 체비쇼프/
│   ├── dolph_chebyshev.m
│   ├── 체비쇼프 n다르게.png
│   └── 체비쇼프 sll다르게.png
├── 프젝2. 모노펄스/
│   ├── dolph_chebyshev_monopulse.m
│   └── result_01~09_*.png (9개 파라미터 변화 실험 결과)
├── 프젝3. 고장소자 5개/
│   ├── simulate_custom_5_failures.m
│   └── case_1~4_original.png (8개 Case 비교)
└── 프젝4. NULA/
    ├── ls_synthesis_1.py
    ├── ls_synthesis_5.py
    ├── ls_synthesis_n.py
    ├── LS_Trend_Analysis.md
    └── ls_pattern_comparison*.png, ls_boxplot*.png, ls_trend_n_v2.png
```

