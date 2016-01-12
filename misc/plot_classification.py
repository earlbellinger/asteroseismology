import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

np.random.seed(13)

output_dir = 'plots/classification'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

N=100
x = np.random.random(N)
err = np.random.normal(0, 2, N)

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111)
ax.plot(x[err>0], x[err>0]+err[err>0], 'o', ms=20, color='#C87A2F')
ax.plot(x[err<0], x[err<0]+err[err<0], '*', ms=20, color='#23373b')
ax.plot([-0.2,1.2],[-0.2,1.2],'k-')
ax.plot([-0.2,1.2],[-0.25,1.25],'k--')
ax.plot([-0.2,1.2],[-0.15,1.15],'k--')
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.xlim([-0.2, 1.2])
plt.ylim([-0.2, 1.2])
fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'classification.png'), 
      bbox_inches='tight')


fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111)
y = x + err
xor = (x < 0.5) & (y < 0.5) | (x > 0.5) & (y > 0.5)
ax.plot(x[xor], y[xor], 'o', ms=20, color='#C87A2F')
ax.plot(x[np.logical_not(xor)], 
        y[np.logical_not(xor)], '*', ms=20, color='#23373b')
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.xlim([-0.2, 1.2])
plt.ylim([-0.2, 1.2])
fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'classification-xor.png'), 
      bbox_inches='tight')



