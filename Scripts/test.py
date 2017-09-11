import numpy as np
import matplotlib.pyplot as plt


phi = np.linspace(0, 2*np.pi, 100)
x = np.cos(phi)
y = np.sin(phi)
plt.plot(x, y, label=r'$\Sigma$')


X = np.linspace(-1, 1, 10)
Y = X + np.random.rand(10)*0.1
E = np.random.rand(10)*0.2
plt.errorbar(X, Y, yerr=E, marker='.', ls=":", label='$Z_0^\mathrm{DMC} $')



plt.title(r'$ p = {0} $'.format(0), fontsize=18)
plt.ylabel('$ E_0 $', fontsize=18)
plt.xlabel(r'$ \alpha $', fontsize=18)
plt.legend(loc=2, fontsize=18)

plt.xlim([-2, 2])

plt.tight_layout()
plt.savefig('plot.pdf')
plt.show()