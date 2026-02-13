import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors

# ----------------------------
# 채도 낮추는 함수
# ----------------------------
def desaturate_cmap(cmap_name="turbo", saturation_scale=0.6):
    base = plt.get_cmap(cmap_name)
    colors = base(np.linspace(0, 1, 256))

    new_colors = []
    for r, g, b, a in colors:
        h, s, v = mcolors.rgb_to_hsv([r, g, b])
        s *= saturation_scale   # 채도 감소
        r2, g2, b2 = mcolors.hsv_to_rgb([h, s, v])
        new_colors.append([r2, g2, b2, a])

    return ListedColormap(new_colors)

# ----------------------------
# 설정
# ----------------------------
num_molecules = 20
latent_dim = 40

top_rows = 12
left_cols = 10

np.random.seed(0)
data = np.random.rand(num_molecules, latent_dim)

top_block = data[:top_rows, :left_cols]
bottom_block = data[-1:, :left_cols]
right_top = data[:top_rows, -1:]
right_bottom = data[-1:, -1:]

# 👇 채도 낮춘 colormap
cmap = desaturate_cmap("turbo", saturation_scale=0.55)

# ----------------------------
# GridSpec
# ----------------------------
fig = plt.figure(figsize=(4, 6))
gs = GridSpec(
    2, 2,
    width_ratios=[left_cols, 1],
    height_ratios=[top_rows, 1],
    wspace=0.2,
    hspace=0.18
)

ax1 = fig.add_subplot(gs[0, 0])
sns.heatmap(top_block, cmap=cmap, cbar=False,
            xticklabels=False,
            yticklabels=[str(i+1) for i in range(top_rows)],
            ax=ax1)

ax2 = fig.add_subplot(gs[1, 0])
sns.heatmap(bottom_block, cmap=cmap, cbar=False,
            xticklabels=False,
            yticklabels=["N"],
            ax=ax2)

ax3 = fig.add_subplot(gs[0, 1])
sns.heatmap(right_top, cmap=cmap, cbar=False,
            xticklabels=False, yticklabels=False, ax=ax3)

ax4 = fig.add_subplot(gs[1, 1])
sns.heatmap(right_bottom, cmap=cmap, cbar=False,
            xticklabels=False, yticklabels=False, ax=ax4)

for ax in [ax1, ax2, ax3, ax4]:
    for spine in ax.spines.values():
        spine.set_visible(False)

plt.show()
