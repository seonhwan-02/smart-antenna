import numpy as np
import matplotlib.pyplot as plt
import os

out_dir = r"C:\Users\kimse\.gemini\antigravity\brain\9142a3ae-c2a8-4506-88ee-82c0717651e9\scratch"
os.makedirs(out_dir, exist_ok=True)

# 1. ULA Chebyshev Pattern (Desired)
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

# 2. Monte Carlo: 1 Spacing Perturbed
np.random.seed(42)
psll_no_ls = []
psll_ls = []

print("Running Monte Carlo with LS Optimization...")
for _ in range(1000):
    d = np.ones(N-1) * 0.5
    idx = np.random.randint(0, N-1)
    d[idx] = np.random.uniform(0.2, 0.7)
    
    x_pert = np.zeros(N)
    for i in range(1, N): x_pert[i] = x_pert[i-1] + d[i-1]
    x_pert -= np.mean(x_pert)
    
    A_pert = np.exp(1j * 2 * np.pi * np.outer(np.sin(theta_rad), x_pert))
    
    # Method 1: No LS (Original weights)
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
print(f"Mean PSLL -> No LS: {mean_no_ls:.2f}dB | With LS: {mean_ls:.2f}dB")

# Example Pattern Comparison
d = np.ones(N-1) * 0.5
d[7] = 0.7 # Worst case roughly (center spacing perturbed)
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
plt.ylim([-40, 0]); plt.xlim([-90, 90])
plt.title('Pattern Restoration using Least Squares Synthesis (1 Spacing Perturbed)')
plt.xlabel('Angle (Degrees)')
plt.ylabel('Magnitude (dB)')
plt.legend()
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'ls_pattern_comparison.png'), dpi=150)
plt.close()

plt.figure(figsize=(8, 6))
box = plt.boxplot([psll_no_ls, psll_ls], patch_artist=True, labels=['Original Weights (No LS)', 'Optimized Weights (LS)'])
colors = ['#e74c3c', '#3498db']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
plt.title('Statistical SLL Improvement with LS (1000 Monte Carlo trials)')
plt.ylabel('Peak Sidelobe Level (dB)')
plt.axhline(-20, color='black', linestyle=':', linewidth=2, label='Target SLL (-20dB)')
plt.legend(loc='upper right')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'ls_boxplot.png'), dpi=150)
plt.close()
