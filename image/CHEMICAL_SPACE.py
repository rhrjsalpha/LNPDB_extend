import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap

# ----------------------------
# 설정
# ----------------------------
num_molecules = 100000
latent_dim = 128

np.random.seed(42)
latent_vectors = np.random.rand(num_molecules, latent_dim)

# 예시 property (예: ΔG)
property_values = np.random.normal(loc=-8, scale=1.5, size=num_molecules)

# ----------------------------
# UMAP 차원 축소
# ----------------------------
reducer = umap.UMAP(
    n_neighbors=15,
    min_dist=0.2,
    metric="cosine",
    random_state=42
)

embedding_2d = reducer.fit_transform(latent_vectors)

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(8,6))

scatter = plt.scatter(
    embedding_2d[:,0],
    embedding_2d[:,1],
    c=property_values,
    cmap="viridis",
    s=40,
    alpha=0.85
)

plt.colorbar(scatter, label="ΔG (kcal/mol)")
plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.title("Chemical Space Projection from Latent Embedding")

plt.tight_layout()
plt.show()
