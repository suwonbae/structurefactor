import numpy as np
import matplotlib.pyplot as plt

#plt.rc('font', size=13, family="Times New Roman")
#plt.rcParams["text.usetex"] = True

class Sq:

    def __init__(self, f_name):
        self.f_name = f_name
        res = np.loadtxt(f_name)
        self.calculate(res)

    def calculate(self, res, **kwargs):

        q = res[:,:3]
        modulus = np.linalg.norm(q, axis=1)
        res = res[np.argsort(modulus)]

        q = res[:,:3]
        modulus = np.linalg.norm(q, axis=1)

        val = res[:,3]

        average = []
        denom = 0
        modulus_tmp = 0
        sq = 0
        for i in range(len(res)-1):
            denom += 1
            modulus_tmp += modulus[i]
            sq += val[i]

            arc_A = modulus[i]
            arc_B = modulus[i+1]

            if ( np.abs(arc_B**2 - arc_A**2) > 0.001 ):
                average.append([modulus_tmp/denom, sq/denom])

                denom = 0
                modulus_tmp = 0
                sq = 0

        self.res = np.asarray(average)
        
lists = {}

freqs = [0, 1, 2, 3, 4]

for ind, freq in enumerate(freqs):
    lists.update({ind: Sq('sq_{}.txt'.format(freq))})
    #a.update({ind: Sq('sq.txt'.format(freq))})

res = np.empty((len(freqs), len(lists[0].res[:,0])))

for ind, freq in enumerate(freqs):
    res[ind,:] = lists[ind].res[:,1]

colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(a)))

fig, ax = plt.subplots(figsize=(4,3))
#ax.plot(a[0].res[:,0], a[0].res[:,1], lw=1)
ax.plot(lists[0].res[:,0], res.mean(axis=0), 'o', ms=0.2)
ax.set_xlim(0, 20)
ax.set_xlabel(r'$k\sigma$')
ax.set_ylabel(r'$S(k)$')
ax.set_ylim(0, 3)
ax.set_yticks([0, 1, 2, 3])
plt.tight_layout(pad=1,h_pad=None,w_pad=None,rect=None)
plt.savefig('sq.png', dpi=200)


'''
freqs = [*np.linspace(0, 1000, 11, dtype=int),
        *np.linspace(2000, 10000, 9, dtype=int),
        *np.linspace(20000, 100000, 9, dtype=int),
        *np.linspace(200000, 1000000, 9, dtype=int)]
        '''
