import numpy as np
import matplotlib.pyplot as plt
import os

out_dir = r"C:\Users\kimse\.gemini\antigravity\brain\9142a3ae-c2a8-4506-88ee-82c0717651e9\scratch"

N = 15
SLL_dB = 25
R = 10**(SLL_dB / 20.0)
z0 = np.cosh(np.arccosh(R) / (N - 1))
W_k = np.zeros(N, dtype=complex)
for k in range(N):
    u_k = (2 * np.pi * k) / N
    x_k = z0 * np.cos(u_k / 2)
    if np.abs(x_k) <= 1: W_k[k] = np.cos((N - 1) * np.arccos(x_k))
    else: W_k[k] = np.cosh((N - 1) * np.arccosh(np.complex128(x_k)))

w_ula = np.zeros(N, dtype=complex)
for n in range(N):
    sum_val = 0
    for k in range(N): sum_val += W_k[k] * np.exp(1j * 2 * np.pi * k * (n - (N+1)/2) / N)
    w_ula[n] = np.real(sum_val)
w_ula = w_ula / np.max(w_ula)

theta = np.linspace(-90, 90, 721)
theta_rad = np.radians(theta)
x_ula = np.arange(N) * 0.5
x_ula -= np.mean(x_ula)
A_ula = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_ula))
F_des = A_ula @ w_ula

n_values = np.arange(1, 15)
mean_psll_no_ls = []
mean_psll_ls = []

np.random.seed(42)

for num_pert in n_values:
    psll_no_ls = []
    psll_ls = []
    for _ in range(500): # 500 trials per n point
        d = np.ones(N-1) * 0.5
        idx = np.random.choice(N-1, num_pert, replace=False)
        for i in idx:
            d[i] = np.random.uniform(0.2, 0.7)
        
        x_pert = np.zeros(N)
        for i in range(1, N): x_pert[i] = x_pert[i-1] + d[i-1]
        x_pert -= np.mean(x_pert)
        
        A_pert = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_pert))
        
        # No LS
        F_no_ls = A_pert @ (w_ula / np.max(w_ula))
        mag_no_ls = 20 * np.log10(np.abs(F_no_ls) / np.max(np.abs(F_no_ls)))
        
        # LS Synthesis
        w_ls, _, _, _ = np.linalg.lstsq(A_pert, F_des, rcond=None)
        F_ls = A_pert @ (w_ls / np.max(np.abs(w_ls)))
        mag_ls = 20 * np.log10(np.abs(F_ls) / np.max(np.abs(F_ls)))
        
        mask = np.abs(theta) > 15
        psll_no_ls.append(np.max(mag_no_ls[mask]))
        psll_ls.append(np.max(mag_ls[mask]))
        
    mean_psll_no_ls.append(np.mean(psll_no_ls))
    mean_psll_ls.append(np.mean(psll_ls))
    print(f"n={num_pert}: No LS {mean_psll_no_ls[-1]:.2f}dB, LS {mean_psll_ls[-1]:.2f}dB")

plt.figure(figsize=(10, 6))
plt.plot(n_values, mean_psll_no_ls, color='#e74c3c', marker='o', linewidth=2, label='Original Weights (No LS)')
plt.plot(n_values, mean_psll_ls, color='#3498db', marker='s', linewidth=2, label='Optimized Weights (LS)')
plt.fill_between(n_values, mean_psll_no_ls, mean_psll_ls, color='gray', alpha=0.1, label='Improvement Gap')
plt.axhline(-25, color='k', linestyle=':', linewidth=2, label='Target SLL (-25dB)')
plt.title('Mean Sidelobe Level vs Number of Perturbed Spacings (n)', fontsize=14, fontweight='bold')
plt.xlabel('Number of Random Spacings (n out of 14)', fontsize=12)
plt.ylabel('Mean Peak Sidelobe Level (dB)', fontsize=12)
plt.xticks(n_values)
plt.legend(loc='lower left')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'ls_trend_n.png'), dpi=150)
plt.close()
