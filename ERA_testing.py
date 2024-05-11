import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ERA(X, r):
    C = np.dot(X.T, X) / (X.shape[0] - 1)
    U, Sigma, Vt = np.linalg.svd(C)
    U_r = U[:, :r]
    EOFs = np.dot(X, U_r)
    PCs = np.dot(U_r.T, X.T)

    return EOFs, PCs

file_path = "p_lock_offset.txt"
df = pd.read_csv(file_path, delim_whitespace=True)

data_for_ERA = df.drop(columns=['time'])

# Perform ERA analysis
r = 50  # Truncation rank
EOFs, PCs = ERA(data_for_ERA.values, r)

plt.figure(figsize=(12, 6))
for i in range(3):
    plt.subplot(1, 3, i + 1)
    plt.plot(EOFs[:, i], label=f"Mode {i + 1}")
    plt.xlabel("Index")
    plt.ylabel("Amplitude")
    plt.legend()
plt.tight_layout()
plt.show()
