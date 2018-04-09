

import numpy
import matplotlib.pyplot as P

def mklaplace(sigma):

    def laplace(x):

        return numpy.exp(-numpy.abs(x)/sigma)/(2.0 * sigma)

    return laplace

if __name__ == '__main__':

    ch = numpy.loadtxt('coeffhist.txt')

    fig, ax = P.subplots()

    #
    # Ensemble distribution
    #
    ax.plot(ch[:,0], ch[:,1], 'k-')

    #
    # Expected distribution
    #
    sigma = 0.1
    target = mklaplace(sigma)(ch[:,0])

    ax.plot(ch[:,0], target, 'r:')

    ax.set_xlabel('Coeff. Magnitude')
    ax.set_ylabel('Probability')

    P.show()
    

                                  
