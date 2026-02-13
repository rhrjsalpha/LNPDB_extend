import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Data (예시용)
# -----------------------
x = np.linspace(0, 10, 600)

# mean: 파형 형태(예시)
mean = 0.75*np.sin(0.75*x) + 0.25*np.sin(1.6*x + 0.8)

# uncertainty: 위치에 따라 폭이 변하도록(예시)
std = 0.18 + 0.22*(np.sin(0.55*x + 1.2)**2)

upper = mean + std
lower = mean - std

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(8, 4), dpi=200)

ax.fill_between(x, lower, upper, alpha=0.25, label="Uncertainty")
ax.plot(x, mean, linewidth=3, label="Mean")

# 스타일: 축 제거, 여백 최소
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)

ax.legend(frameon=False, loc="upper right", fontsize=14)
plt.tight_layout()

# 저장(원하면)
# plt.savefig("gp_mean_uncertainty.png", dpi=600, bbox_inches="tight", transparent=True)

plt.show()
