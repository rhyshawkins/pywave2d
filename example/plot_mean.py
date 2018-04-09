
import numpy
import matplotlib.pyplot as P

if __name__ == '__main__':

    mean = numpy.loadtxt('mean.txt')
    truemodel = numpy.loadtxt('image_true.txt')

    fig, ax = P.subplots(1, 2)

    vmin = min([numpy.min(mean), numpy.min(truemodel)])
    vmax = max([numpy.max(mean), numpy.max(truemodel)])
    cmap = 'inferno'
    
    ax[0].imshow(mean, vmin = vmin, vmax = vmax, cmap = cmap)
    ax[1].imshow(truemodel, vmin = vmin, vmax = vmax, cmap = cmap)

    P.show()
    
    
