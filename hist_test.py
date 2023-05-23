import numpy as np
import matplotlib.pyplot as plt

a1 = np.array([1,2,3,4,5,6,7,8])
a2 = np.array([[1,2,3,4,5],[3,4,5,6,7]])
np.hstack(a2)

fig, ax = plt.subplots(figsize=(6,4))
ax.hist(a2, rwidth=.9, color='blue')

plt.savefig('../../../plots/hist_test.png', dpi=300, transparent=False)
plt.close()
