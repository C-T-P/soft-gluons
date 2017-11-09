#!/usr/bin/env python2
#@LOADMPIFORPY@
import sys
import math as m
sys.path.append("/home/preuss10/Compiled/lib/python2.7/site-packages/")
import Sherpa

def calcMatEl (process, p_1, p_2, p_3, p_4, p_s):
    # Add this to the execution arguments to prevent Sherpa from starting the cross section integration
    sys.argv.append('INIT_ONLY=2')

    Generator=Sherpa.Sherpa()
    
    if process[0] == 1:
        if process[1] == 1:
            pid_1 = 1
            pid_2 = -1
            pid_3 = 1 
            pid_4 = -1
            
        elif process[1] == 2:
            pid_1 = 1
            pid_2 = -1
            pid_3 = 2
            pid_4 = -2
            
        elif process[1] == 3:
            pid_1 = 2
            pid_2 = -1
            pid_3 =  2
            pid_4 = -1
            
    elif process[0] == 2:
        if process[1] == 1:
            pid_1 = 1
            pid_2 = 1
            pid_3 = 1
            pid_4 = 1
            
        elif process[1] == 2:
            pid_1 = 1
            pid_2 = 2
            pid_3 = 1
            pid_4 = 2
        
    elif process[0] == 3:
        pid_1 = 1
        pid_2 = 21
        pid_3 = 1
        pid_4 = 21
            
    elif process[0] == 4:
        if process[1] == 1:
            pid_1 = 1
            pid_2 = -1
            pid_3 = 21 
            pid_4 = 21
            
        elif process[1] == 2:
            pid_1 = 21
            pid_2 = 21
            pid_3 = 1
            pid_4 = -1
            
    elif process[0] == 5:
        pid_1 = 21
        pid_2 = 21
        pid_3 = 21 
        pid_4 = 21
        
    else:
        pid_1 = 0
        pid_2 = 0
        pid_3 = 0
        pid_4 = 0
    
    try:
        Generator.InitializeTheRun(len(sys.argv),sys.argv)
        Process=Sherpa.MEProcess(Generator)

        # Incoming flavors must be added first!
        Process.AddInFlav(pid_1);
        Process.AddInFlav(pid_2);
        Process.AddOutFlav(pid_3);
        Process.AddOutFlav(pid_4);
        Process.AddOutFlav(21);
        Process.Initialize();

        # Momentum setting via list of floats
        Process.SetMomenta([p_1,
                            p_2,
                            p_3,
                            p_4,
                            p_s])
        MatEl = Process.CSMatrixElement()
        print '\nSquared ME: ', MatEl

    except Sherpa.Exception as exc:
        print exc
        exit(1)
    
    return MatEl;
