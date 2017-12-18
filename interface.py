#!/usr/bin/env python2
#@LOADMPIFORPY@
import sys
import math as m
sys.path.append("/home/preuss10/Compiled/lib/python2.7/site-packages/")
import Sherpa
import numpy as np

def calcMatEl (process, K, lambda_s, E):
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
        title = "$q\\bar{q} \\to q\\bar{q}(g)$"
            
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
        title = "$qq \\to qq(g)$"
        
    elif process[0] == 3:
        pid_1 = 1
        pid_2 = 21
        pid_3 = 1
        pid_4 = 21
        
        title = "$qg \\to qg(g)$"
            
    elif process[0] == 4:
        if process[1] == 1:
            pid_1 = 1
            pid_2 = -1
            pid_3 = 21 
            pid_4 = 21
            title = "$q\\bar{q} \\to gg(g)$"
            
        elif process[1] == 2:
            pid_1 = 21
            pid_2 = 21
            pid_3 = 1
            pid_4 = -1
            title = "$gg \\to q\\bar{q}(g)$"
            
    elif process[0] == 5:
        pid_1 = 21
        pid_2 = 21
        pid_3 = 21 
        pid_4 = 21
        title = "$gg \\to gg(g)$"
        
    else:
        title = "not defined"   
        
    # Debugging only
    #with open('Run.dat', 'w') as f:
        #f.write("(run){\n\tEVENTS %d\n\tSHERPA_LDADD SherpaMain;\n\tSCALES VAR{sqr(91.18)}\n}(run)" % (0))
        #f.write("\n\n(beam){\n\tBEAM_1  %d; BEAM_ENERGY_1  %.1f;\n\tBEAM_2  %d; BEAM_ENERGY_2  %.1f;\n}(beam)" % (pid_1, E, pid_2, E))
        #f.write("\n\n(isr){\n\tPDF_LIBRARY None;\n}(isr)")
        #f.write("\n\n(processes){\n\tProcess %d %d -> %d %d;\n\tME_Generator Amegic;\n\tOrder (2,0);\n\tEnd process;\n}(processes)" % (pid_1, pid_2, pid_3, pid_4))    
    
    # 2 -> 2 + 1 Process
    with open('Run.dat', 'w') as f:
        f.write("(run){\n\tEVENTS %d\n\tSHERPA_LDADD SherpaMain;\n\tSCALES VAR{sqr(91.18)}\n}(run)" % (0))
        f.write("\n\n(beam){\n\tBEAM_1  %d; BEAM_ENERGY_1  %.1f;\n\tBEAM_2  %d; BEAM_ENERGY_2  %.1f;\n}(beam)" % (pid_1, E, pid_2, E))
        f.write("\n\n(isr){\n\tPDF_LIBRARY None;\n}(isr)")
        f.write("\n\n(processes){\n\tProcess %d %d -> %d %d %d;\n\tME_Generator Amegic;\n\tOrder (3,0);\n\tEnd process;\n}(processes)" % (pid_1, pid_2, pid_3, pid_4, 21))
    
    try:
        Generator.InitializeTheRun(len(sys.argv),sys.argv)
        Process=Sherpa.MEProcess(Generator)

        # Incoming flavours must be added first!
        Process.AddInFlav(pid_1);
        Process.AddInFlav(pid_2);
        Process.AddOutFlav(pid_3);
        Process.AddOutFlav(pid_4);
        Process.AddOutFlav(21);
        Process.Initialize();

        # Momentum setting
        phi_43 = m.pi * 3./2.
        phi_44 = m.pi * 5./2.
        phi_H = m.pi/10.
        phi_s = m.pi/7.
        MatEl = np.zeros((K.size,lambda_s.size))
        
        for i in range(0, K.size):
            N = K[i]
            j = 0
            for l_s in lambda_s:
                k_s = 2.*E*l_s
                p_1 = [E,0.,0.,E]
                p_2 = [E,0.,0.,-E]    
                p_3 = [E,E*m.cos(phi_43+N*phi_H),E*m.sin(phi_43+N*phi_H),0]
                p_4 = [E,E*m.cos(phi_44+N*phi_H),E*m.sin(phi_44+N*phi_H),0]
                p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]
                
                Process.SetMomenta([p_1,
                                    p_2,
                                    p_3,
                                    p_4,
                                    p_s
                                    ])
                MatEl[i,j] = Process.CSMatrixElement()
                j = j+1
        

    except Sherpa.Exception as exc:
        print exc
        exit(1)
    
    return title, MatEl;
