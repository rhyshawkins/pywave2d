
import glob
import sys

sys.path.append('..')
sys.path.extend(glob.glob('../build/lib*'))
                
import numpy

import pywave2d
import chainhistory

if __name__ == '__main__':

    history = chainhistory.ChainHistory()

    history.load('chain_history.dat')

    degreex = 5
    degreey = 5
    kmax = 100
    skip = 50000
    thin = 5

    #
    # Create the Trans-d Tree Wavelet Object
    #
    wave2d = pywave2d.pyWave2D(degreex, degreey, pywave2d.WAVELET_CDF97, kmax)

    width = wave2d.get_width()
    height = wave2d.get_height()

    mean = numpy.zeros((width, height))
    meann = 0

    likehist = []
    normhist = []
    errscalehist = []
    priorscalehist = []
    khist = []

    coeff_min = -1.0
    coeff_max = 1.0
    coeff_bins = 200
    coeff_hist, edges = numpy.histogram([], bins = coeff_bins, range = (coeff_min, coeff_max));
    
    for i, model in enumerate(history):

        (like, norm, errscale, priorscale), indices, values = model

        if (i > skip) and ((i - skip) % thin == 0):

            wave2d.set_model_indices(indices)
            wave2d.set_model_values(values)
            wave2d.model_commit()

            meann = meann + 1
            
            image = wave2d.get_image()
            delta = image - mean
            mean = mean + delta/float(meann)

            likehist.append(like)
            normhist.append(norm)
            errscalehist.append(errscale)
            priorscalehist.append(priorscale)
            khist.append(len(indices))

            coeff_hist_model, _ = numpy.histogram(values, bins = 200, range = (coeff_min, coeff_max))

            coeff_hist = coeff_hist + coeff_hist_model

        if ((i + 1) % 10000 == 0):
            print i + 1, '...'
            

    numpy.savetxt('mean.txt', mean)

    numpy.savetxt('like.txt', numpy.array(likehist))
    numpy.savetxt('norm.txt', numpy.array(normhist))
    numpy.savetxt('errscale.txt', numpy.array(errscalehist))
    numpy.savetxt('priorscale.txt', numpy.array(priorscalehist))
    numpy.savetxt('khist.txt', numpy.array(khist))

    centres = (edges[1:] + edges[:-1])/2
    S = numpy.sum(coeff_hist)
    dx = (coeff_max - coeff_min)/float(coeff_bins)
    N = len(centres)
    ncoeff_hist = numpy.array(coeff_hist, dtype = float)/(dx * S)
    numpy.savetxt('coeffhist.txt', numpy.hstack((centres.reshape((N,1)), ncoeff_hist.reshape((N,1)))))
    
    #
    # Save final model
    #
    f = open('final_model_post.txt', 'w')
    for (i, v) in zip(indices, values):
        f.write('%d %15.9f\n' % (i, v))
    f.close()
