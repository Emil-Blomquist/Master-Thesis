import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))

n = 5
dt = 0.02
t = np.linspace(0, (n - 1)*dt, n) + 0.5*dt

print(t)

histogramColor = [1, 0, 0]

for i in range(0, n):
  Hi = 2/dt**0.5 * ( (i + 1)**0.5 - i**0.5 )
  hist, = plt.plot([i*dt, i*dt, (i + 1)*dt, (i + 1)*dt], [0, Hi, Hi, 0], color=histogramColor)


x = np.linspace(0.001, (n + 1)*dt, 1000)
func, = plt.plot(x, x**-0.5, ':')
true, = plt.plot(t, t**-0.5, 'o')


plt.legend([func, true, hist], [r'$ 1/\sqrt{x}$', 'True value', 'Histogram approximation'], fontsize=18)
e = 0.001
plt.xlim([0 - e, n*dt + e])
plt.ylim([0, 20])
plt.tight_layout()
plt.savefig('plots/plot.pdf')
plt.show() 

