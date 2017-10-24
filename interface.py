#!/usr/bin/env python2
#@LOADMPIFORPY@
import sys
import math as m
sys.path.append("/home/preuss10/Compiled/lib/python2.7/site-packages/")
import Sherpa

def calcMatEl (p_1, p_2, p_3, p_4, p_s):
    # Add this to the execution arguments to prevent Sherpa from starting the cross section integration
    sys.argv.append('INIT_ONLY=2')

    Generator=Sherpa.Sherpa()
    
    try:
        Generator.InitializeTheRun(len(sys.argv),sys.argv)
        Process=Sherpa.MEProcess(Generator)

        # Incoming flavors must be added first!
        Process.AddInFlav(1);
        Process.AddInFlav(-1);
        Process.AddOutFlav(1);
        Process.AddOutFlav(-1);
        Process.AddOutFlav(21);
        Process.Initialize();

        # First argument corresponds to particle index:
        # index 0 correspons to particle added first, index 1 is the particle added second, and so on...
        #Process.SetMomentum(0, 50.,0.,0.,50.)
        #Process.SetMomentum(1, 50.,0.,0.,-50.)
        #Process.SetMomentum(2, 50.,0.,50.,0.)
        #Process.SetMomentum(3, 50.,0.,-50.,0.)
        #print '\nSquared ME: ', Process.CSMatrixElement()

        # Momentum setting via list of floats
        Process.SetMomenta([p_1,
                            p_2,
                            p_3,
                            p_4,
                            p_s])
        MatEl = Process.CSMatrixElement()
        print '\nSquared ME: ', MatEl
        
        # Random momenta
        #E_cms = 500.0
        #tp = Process.TestPoint(E_cms)
        #print '\nRandom test point: ', tp[0], tp[1], tp[2], tp[3]
        #print 'Squared ME: ', Process.CSMatrixElement(), '\n'

    except Sherpa.Exception as exc:
        print exc
        exit(1)
    
    return MatEl;
