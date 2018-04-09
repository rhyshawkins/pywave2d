
import numpy
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == '__main__':

    thin = 5
    skip = 50000

    ps = numpy.loadtxt('priorscale.txt')

    fig, ax = P.subplots()

    divider = make_axes_locatable(ax)
    err = divider.append_axes('right', 1.2, pad = 0.25, sharey = ax)
    err.yaxis.tick_right()

    bins = 250
    scalemin = 0.0
    scalemax = 0.25

    hist, edges = numpy.histogram(ps, bins = bins, range = (scalemin, scalemax))
    centres = (edges[1:] + edges[:-1])/2.0

    it = numpy.arange(len(ps)) * thin + skip
    
    ax.plot(it, ps, 'k-')
    ax.set_ylim(scalemin, scalemax)

    ax.set_xlabel('Interation')
    ax.set_ylabel('Prior Width')

    err.plot(hist, centres, 'k-')
    err.axhline(y = numpy.mean(ps), linestyle = 'dashed', color = 'red')
    err.axhline(y = 0.1, linestyle = 'dashed', color = 'blue')
    err.set_xticks([])

    P.show()
    

                                  
