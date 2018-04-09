
import glob
import sys
sys.path.append('..')
sys.path.extend(glob.glob('../build/lib*'))

print sys.path

import numpy

import pywave2d
import chainhistory

def likelihood(wave2d, xmin, xmax, ymin, ymax, data, sigma):

    rows, cols = data.shape

    image = wave2d.get_image()#numpy.array(, copy = True)
    width, height = image.shape

    like = 0.0
    
    for r in range(rows):
        x, y, z = data[r, :]

        i = int((x - xmin)/(xmax - xmin) * float(width))
        if (i < 0 or i >= width):
            raise Exception('x out of range: %f' % x)

        j = int((y - ymin)/(ymax - ymin) * float(height))
        if (j < 0 or j >= height):
            raise Exception('y out of range: %f' % y)

        pred_z = image[j, i]

        dz = z - pred_z

        like = like + dz*dz/(2.0 * sigma*sigma)

    norm = float(rows) * numpy.log(sigma)
    
    return (like, norm)

DELTA_INIT=0
DELTA_ADD=1
DELTA_DEL=2
DELTA_CHANGE=3
    
def summary(name, propose, accept):

    print name, ':',
    
    for p, a in zip(propose, accept):

        r = 0.0
        if p > 0:
            r = float(a)/float(p) * 100.0
            
        print '%7.3f ' % r,

    print ''


if __name__ == '__main__':

    kmax = 100
    iterations = 100000
    verbosity = 1000

    #
    # Configure/Load the data
    #
    xmin = -10.0
    xmax = 10.0
    ymin = -10.0
    ymax = 10.0
    
    data = numpy.loadtxt('data.txt', skiprows = 1)

    sigma = 0.1

    #
    # Create the Trans-d Tree Wavelet Object
    #
    wave2d = pywave2d.pyWave2D(5, 5, pywave2d.WAVELET_CDF97, kmax)

    #
    # Initialise the model (parameter is mean of image)
    #
    wave2d.initialise(0.0)
    
    #
    # Load prior/proposal file
    #
    wave2d.load_priorproposal('priorproposal.txt')

    #
    # Initial likelihood
    #
    current_like, current_norm = likelihood(wave2d, xmin, xmax, ymin, ymax, data, sigma)
    print 'Initial likelihood: %f %f' % (current_like, current_norm)

    wave2d.model_retrieve()
    current_indices = list(wave2d.get_model_indices())
    current_values = list(wave2d.get_model_values())

    history = chainhistory.ChainHistory()

    current_err_scale = 1.0
    current_prior_scale = wave2d.get_prior_scale()
    
    history.add_initial(current_like, current_norm, current_err_scale, current_prior_scale, current_indices, current_values)

    hierarchical_error_std = 0.1
    ph = [0]
    ah = [0]

    hierarchical_prior_std = 0.1
    pp = [0]
    ap = [0]
    
    for i in xrange(iterations):

        #
        # Perform a Birth/Death/Value proposal. Invalid proposals
        # are possible (prior violation, death with 1 parameter,
        # birth with kmax parameters etc) and are immediately
        # rejected.
        #
        ppratio = wave2d.propose()

        valid = ppratio > pywave2d.pyWave2D.INVALID_PROPOSAL
        if (valid):

            proposed_like, proposed_norm = likelihood(wave2d, xmin, xmax, ymin, ymax, data, current_err_scale * sigma)
            
            u = numpy.log(numpy.random.uniform())

            if (u < ((current_like + current_norm) - (proposed_like + proposed_norm) +
                     ppratio)):

                wave2d.accept()
                current_like = proposed_like
                current_norm = proposed_norm
                
            else:

                wave2d.reject()

        #
        # Perform a hierarchical scale proposal
        #
        if (hierarchical_error_std > 0.0):

            ph[0] = ph[0] + 1
            proposed_err_scale = current_err_scale +numpy.random.normal() * hierarchical_error_std
            if proposed_err_scale > 0.0 and proposed_err_scale < 5.0:

                proposed_like, proposed_norm = likelihood(wave2d, xmin, xmax, ymin, ymax, data, proposed_err_scale * sigma)

                u = numpy.log(numpy.random.uniform())

                if (u < ((current_like + current_norm) - (proposed_like + proposed_norm))):

                    current_err_scale = proposed_err_scale
                    current_like = proposed_like
                    current_norm = proposed_norm
                
                    ah[0] = ah[0] + 1
                


        #
        # Perform a hierarchical prior proposal
        #
        if (hierarchical_prior_std > 0.0):

            pp[0] = pp[0] + 1
            proposed_prior_scale = current_prior_scale * numpy.exp(numpy.random.normal() * hierarchical_prior_std)

            priorratio = wave2d.propose_prior_scale(proposed_prior_scale)
            if priorratio > pywave2d.pyWave2D.INVALID_PROPOSAL:
                
                u = numpy.log(numpy.random.uniform())

                if (u < priorratio):

                    current_prior_scale = proposed_prior_scale
                    wave2d.set_prior_scale(proposed_prior_scale)
                
                    ap[0] = ap[0] + 1

        

        #
        # Update status
        #
        if (i + 1) % verbosity == 0:

            print '%6d: Like %10.6f (%10.6f) Err %10.6f Prior %10.6f k %d' % (i + 1,
                                                                              current_like,
                                                                              current_norm,
                                                                              current_err_scale,
                                                                              current_prior_scale,
                                                                              wave2d.get_k())
            summary('Birth', wave2d.pb_counts(), wave2d.ab_counts())
            summary('Death', wave2d.pd_counts(), wave2d.ad_counts())
            summary('Value',  wave2d.pv_counts(), wave2d.av_counts())
            if hierarchical_error_std > 0.0:
                summary('Err  ', ph, ah)
            if hierarchical_prior_std > 0.0:
                summary('Prior', pp, ap)


        #
        # Get the model for saving to ensemble
        #
        wave2d.model_retrieve()
        indices = wave2d.get_model_indices()
        values = wave2d.get_model_values()

        history.add_delta(current_like, current_norm, current_err_scale, current_prior_scale, current_indices, current_values, indices, values)

        current_indices = list(indices)
        current_values = list(values)

#        print wave2d.get_logprior()
        
    #
    # Save final model
    #
    f = open('final_model.txt', 'w')
    for (i, v) in zip(current_indices, current_values):
        f.write('%d %15.9f\n' % (i, v))
    f.close()

    #
    # Save chain history
    #
    history.save('chain_history.dat')
    
    
