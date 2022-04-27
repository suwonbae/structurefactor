import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

class Sq:

    def __init__(self, f_name):
        self.f_name = f_name
        self.res = np.loadtxt(f_name)

    def calculate(self, **kwargs):

        q = self.res[:,0]
        val = self.res[:,1]

        if 'q_num' in kwargs:
            q_num = kwargs['q_num']
        else:
            delta = 0.05
            q_num = int((np.max(q))/delta)

        q_interest = np.linspace(0.3, np.max(q), q_num)

        bin_width = q_interest[1] - q_interest[0]
        self.q = q_interest[1:] - bin_width/2
        self.I = stats.binned_statistic(q, val, statistic='mean', bins=q_interest).statistic


sq = Sq('sq.txt')
sq.calculate()

'''
#for multiple xyz to calculate ensemble average
sq_list = {}

freqs = [0, 1, 2, 3, 4]

for ind, freq in enumerate(freqs):
    sq_list.update({ind: Sq('sq_{}.txt'.format(freq))})
    sq_list[ind].calculate()

res = np.empty((len(freqs), len(sq_list[0].q)))

for ind, freq in enumerate(freqs):
    res[ind,:] = a[ind].I
'''

fig, ax = plt.subplots(figsize=(4,3))
ax.plot(sq.q, sq.I, color='black', lw=1)
'''
#for multiple xyz to calculate ensemble average
for ind in sq_list.keys():
    ax.plot(sq_list[ind].q, sq_list[ind].I, 'o', ms=0.5, color='black', alpha=0.2)
ax.plot(sq_list[ind].q, res.mean(axis=0), lw=1)
'''
ax.set_xlim(0, 20)
ax.set_xlabel(r'$k\sigma$')
ax.set_ylabel(r'$S(k)$')
ax.set_ylim(0, 2)
plt.tight_layout(pad=1,h_pad=None,w_pad=None,rect=None)
plt.savefig('sq.png', dpi=500)

