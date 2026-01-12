import pandas as pd
import matplotlib.pyplot as plt

# This script build the coverage scheme plot from the output given by bowtie

# IMPORT THE DATA
cov = pd.read_csv("coverage.txt", sep="\t", header=None, names=["ref", "pos", "cov"])
plt.figure(figsize=(8,3))

seq = "CAAGGGCTGAACCTAGAATTCTTTTTGGATGAGGCTTTGGAAAAGAATTCTAGGTTCAGCCCTCG"
# Positions of each region
strand5p = seq[:22]
loop = seq[22:65-22]  
strand3p = seq[-22:]

plt.figure(figsize=(8,3))

plt.bar(cov["pos"], cov["cov"], color="teal")
plt.xlabel("pre-miRNA position")
plt.ylabel("Coverage (reads)")
plt.title("Coverage Scheme of pre-miRNA", pad = 30)
plt.tight_layout()

y_offset = cov["cov"].max() * 1.1  # Vertical position above the max
x_positions = cov["pos"]

for i, base in enumerate(seq, start=1):
    # Colour by region
    if i <= 22:
        color = "red"      # 5p
    elif i > 43:
        color = "green"    # 3p
    else:
        color = "orange"   # loop
    plt.text(i, y_offset, base, color=color, ha='center', va='bottom', fontsize=8, fontfamily='monospace')

# Colour legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(color='red', label='5p strand (1–22)'),
    Patch(color='orange', label='Loop (23–43)'),
    Patch(color='green', label='3p strand (44–65)')
]
plt.legend(handles=legend_elements, loc='upper right')


plt.show()
plt.savefig("coverage_scheme.png", dpi=300)
