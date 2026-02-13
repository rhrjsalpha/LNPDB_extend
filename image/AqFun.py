import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Data
# -----------------------
x = np.linspace(0, 10, 600)

mean = 0.75*np.sin(0.75*x) + 0.25*np.sin(1.6*x + 0.8)
std  = 0.18 + 0.22*(np.sin(0.55*x + 1.2)**2)

upper = mean + std
lower = mean - std

# -----------------------
# Acquisition (UCB example)
# -----------------------
kappa = 1.2
acq = mean + kappa*std

acq_norm = (acq - acq.min()) / (acq.max() - acq.min() + 1e-12)

imax = np.argmax(acq_norm)
x_max = x[imax]
y_max = acq_norm[imax]

# -----------------------
# Figure
# -----------------------
fig = plt.figure(figsize=(8, 5), dpi=200)
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.02)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

# ---- 목적함수 (파란 계열)
mean_color = "#1f3b73"
unc_color  = "#4a90e2"

ax1.fill_between(x, lower, upper, color=unc_color, alpha=0.25, label="Uncertainty")
ax1.plot(x, mean, color=mean_color, linewidth=3, label="Mean")

ax1.legend(frameon=False, loc="upper right", fontsize=14)

ax1.set_xticks([]); ax1.set_yticks([])
for s in ax1.spines.values():
    s.set_visible(False)

# ---- Acquisition (초록 계열)
aq_line  = "#1b5e20"
aq_fill  = "#66bb6a"

ax2.fill_between(x, 0, acq_norm, color=aq_fill, alpha=0.25)
ax2.plot(x, acq_norm, color=aq_line, linewidth=3)

# acquisition max (빨간 삼각형)
ax2.scatter([x_max], [y_max], marker="v", s=180,
            color="#c62828", zorder=5)

ax2.text(x_max + 0.2, y_max + 0.05,
         "acquisition max",
         fontsize=16,
         color="#333333",
         va="bottom")

ax2.set_ylim(0, 1.05)
ax2.set_xticks([]); ax2.set_yticks([])
for s in ax2.spines.values():
    s.set_visible(False)

plt.tight_layout()
plt.show()
