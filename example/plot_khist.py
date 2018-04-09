
import numpy
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == '__main__':

    thin = 5
    skip = 50000

    k = numpy.loadtxt('khist.txt')

    fig, ax = P.subplots()

    divider = make_axes_locatable(ax)
    err = divider.append_axes('right', 1.2, pad = 0.25, sharey = ax)
    err.yaxis.tick_right()

    kmin = 0
    kmax = 50

    hist, edges = numpy.histogram(k, bins = kmax - kmin + 1, range = (kmin, kmax))
    centres = (edges[1:] + edges[:-1])/2.0

    it = numpy.arange(len(k)) * thin + skip
    
    ax.plot(it, k, 'k-')
    ax.set_ylim(kmin, kmax)

    ax.set_xlabel('Iteration')
    ax.set_ylabel('No. Coefficients')

    err.plot(hist, centres, 'k-')
    err.axhline(y = numpy.mean(k), linestyle = 'dashed', color = 'red')
    err.set_xticks([])

    P.show()
    

                                  
