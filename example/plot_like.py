
import numpy
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == '__main__':

    thin = 5
    skip = 50000
    ndata = 64
    chi = float(ndata)/2.0

    #
    # Likelihood (with chi^2 target in blue)
    #
    
    like = numpy.loadtxt('like.txt')
    escale = numpy.loadtxt('errscale.txt')
    
    fig, ax = P.subplots(2, 1)
    fig.set_tight_layout(True)

    divider = make_axes_locatable(ax[0])
    err = divider.append_axes('right', 1.2, pad = 0.25, sharey = ax[0])
    err.yaxis.tick_right()

    bins = 250
    likemin = 0.0
    likemax = 50.0

    hist, edges = numpy.histogram(like, bins = bins, range = (likemin, likemax))
    centres = (edges[1:] + edges[:-1])/2.0

    it = numpy.arange(len(like)) * thin + skip
    
    ax[0].plot(it, like, 'k-')
    ax[0].axhline(y = chi, linestyle = 'dashed', color = 'blue')
    
    ax[0].set_ylim(likemin, likemax)

    ax[0].set_xlabel('Interation')
    ax[0].set_ylabel('Likelihood')

    err.plot(hist, centres, 'k-')
    err.axhline(y = numpy.mean(like), linestyle = 'dashed', color = 'red')
    err.axhline(y = chi, linestyle = 'dashed', color = 'blue')
    err.set_xticks([])

    #
    # Hierarchical Error Scale (with 1.0 target in blue)
    #

    divider = make_axes_locatable(ax[1])
    err = divider.append_axes('right', 1.2, pad = 0.25, sharey = ax[1])
    err.yaxis.tick_right()

    bins = 250
    escalemin = 0.5
    escalemax = 1.5

    hist, edges = numpy.histogram(escale, bins = bins, range = (escalemin, escalemax))
    centres = (edges[1:] + edges[:-1])/2.0

    it = numpy.arange(len(escale)) * thin + skip
    
    ax[1].plot(it, escale, 'k-')
    ax[1].axhline(y = 1.0, linestyle = 'dashed', color = 'blue')
    
    ax[1].set_ylim(escalemin, escalemax)

    ax[1].set_xlabel('Interation')
    ax[1].set_ylabel('Hierarchical Error Scale')

    err.plot(hist, centres, 'k-')
    err.axhline(y = numpy.mean(escale), linestyle = 'dashed', color = 'red')
    err.axhline(y = 1.0, linestyle = 'dashed', color = 'blue')
    err.set_xticks([])
    

    P.show()
    

                                  
