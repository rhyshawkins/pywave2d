
import math
import numpy

NSAMPLES = 60
XMIN = -10.0
XMAX = 10.0
YMIN = -10.0
YMAX = 10.0

WIDTH = 32
HEIGHT = 32

def g(xn, yn):

    x = (xn - XMIN)/(XMAX - XMIN) * 2.0 - 1.0
    y = (yn - YMIN)/(YMAX - YMIN) * 2.0 - 1.0

    return math.exp(-0.5 * (x*x + y*y)/(0.25*0.25))

if __name__ == '__main__':


    xp = numpy.random.uniform(low = XMIN, high = XMAX, size = (NSAMPLES,))
    yp = numpy.random.uniform(low = YMIN, high = YMAX, size = (NSAMPLES,))
    
    noise = 0.1

    f = open('data.txt', 'w')
    f.write('%d\n' % NSAMPLES)
    for x, y in zip(xp, yp):
        z = g(x, y) + numpy.random.normal() * noise
        f.write('%15.9f %15.9f %15.9f\n' % (x, y, z))

    f.close()

    f = open('data_true.txt', 'w')
    f.write('%d\n' % NSAMPLES)
    for x, y in zip(xp, yp):
        z = g(x, y)
        f.write('%15.9f %15.9f %15.9f\n' % (x, y, z))

    f.close()

    img = numpy.zeros((WIDTH, HEIGHT))
    for j in range(HEIGHT):

        y = YMIN + (YMAX - YMIN) * (float(j) + 0.5)/float(HEIGHT)
        
        for i in range(WIDTH):
            x = XMIN + (XMAX - XMIN) * (float(i) + 0.5)/float(WIDTH)
            img[i, j] = g(x, y)

    numpy.savetxt('image_true.txt', img)
    

    

    
