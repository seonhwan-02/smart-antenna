import numpy as np
import matplotlib.pyplot as plt
import os

out_dir = r"C:\Users\kimse\.gemini\antigravity\brain\9142a3ae-c2a8-4506-88ee-82c0717651e9\scratch"

N = 15
SLL_dB = 20
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

theta = np.linspace(-90, 90, 1801)
theta_rad = np.radians(theta)
x_ula = np.arange(N) * 0.5
x_ula -= np.mean(x_ula)
A_ula = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_ula))
F_des = A_ula @ w_ula

np.random.seed(42)
psll_no_ls = []
psll_ls = []

num_perturbed = 5

for _ in range(1000):
    d = np.ones(N-1) * 0.5
    idx = np.random.choice(N-1, num_perturbed, replace=False)
    for i in idx:
        d[i] = np.random.uniform(0.2, 0.7)
    
    x_pert = np.zeros(N)
    for i in range(1, N): x_pert[i] = x_pert[i-1] + d[i-1]
    x_pert -= np.mean(x_pert)
    
    A_pert = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_pert))
    
    # Method 1: No LS
    F_no_ls = A_pert @ (w_ula / np.max(w_ula))
    mag_no_ls = 20 * np.log10(np.abs(F_no_ls) / np.max(np.abs(F_no_ls)))
    
    # Method 2: LS Synthesis
    w_ls, _, _, _ = np.linalg.lstsq(A_pert, F_des, rcond=None)
    F_ls = A_pert @ (w_ls / np.max(np.abs(w_ls)))
    mag_ls = 20 * np.log10(np.abs(F_ls) / np.max(np.abs(F_ls)))
    
    mask = np.abs(theta) > 15
    psll_no_ls.append(np.max(mag_no_ls[mask]))
    psll_ls.append(np.max(mag_ls[mask]))

mean_no_ls = np.mean(psll_no_ls)
mean_ls = np.mean(psll_ls)
print(f"5 elements -> Mean PSLL -> No LS: {mean_no_ls:.2f}dB | With LS: {mean_ls:.2f}dB")

# Example worst case for 5 elements
d = np.ones(N-1) * 0.5
d[2] = 0.7; d[5] = 0.7; d[7] = 0.2; d[10] = 0.7; d[12] = 0.2
x_pert = np.zeros(N)
for i in range(1, N): x_pert[i] = x_pert[i-1] + d[i-1]
x_pert -= np.mean(x_pert)
A_pert = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_pert))

F_no_ls = A_pert @ (w_ula / np.max(w_ula))
mag_no_ls = 20 * np.log10(np.abs(F_no_ls) / np.max(np.abs(F_no_ls)))

w_ls, _, _, _ = np.linalg.lstsq(A_pert, F_des, rcond=None)
F_ls = A_pert @ (w_ls / np.max(np.abs(w_ls)))
mag_ls = 20 * np.log10(np.abs(F_ls) / np.max(np.abs(F_ls)))

plt.figure(figsize=(10, 5))
plt.plot(theta, mag_no_ls, color='#e74c3c', linestyle='--', linewidth=1.5, label=f'No LS (PSLL: {np.max(mag_no_ls[np.abs(theta)>15]):.2f}dB)')
plt.plot(theta, mag_ls, color='#3498db', linestyle='-', linewidth=1.5, label=f'With LS (PSLL: {np.max(mag_ls[np.abs(theta)>15]):.2f}dB)')
plt.axhline(-20, color='k', linestyle=':', label='Target SLL (-20dB)')
plt.ylim([-35, 0]); plt.xlim([-90, 90])
plt.title('Pattern Restoration using LS Synthesis (5 Spacings Perturbed)')
plt.xlabel('Angle (Degrees)')
plt.ylabel('Magnitude (dB)')
plt.legend()
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'ls_pattern_comparison_5.png'), dpi=150)
plt.close()

plt.figure(figsize=(8, 6))
box = plt.boxplot([psll_no_ls, psll_ls], patch_artist=True, labels=['Original Weights (No LS)', 'Optimized Weights (LS)'])
colors = ['#e74c3c', '#3498db']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
plt.title('Statistical SLL Improvement with LS (5 Perturbed, 1000 trials)')
plt.ylabel('Peak Sidelobe Level (dB)')
plt.axhline(-20, color='black', linestyle=':', linewidth=2, label='Target SLL (-20dB)')
plt.legend(loc='upper right')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'ls_boxplot_5.png'), dpi=150)
plt.close()
