import math as m
import cmath as cm
import numpy as np
import matplotlib.pyplot as plt
import interface as sh
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def writeRunCard(process, E):
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
            title = "$q\\bar{q} \\to gg(g)"
            
        elif process[1] == 2:
            pid_1 = 21
            pid_2 = 21
            pid_3 = 1
            pid_4 = -1
            title = "$qq \\to q\\bar{q}(g)"
            
    elif process[0] == 5:
        pid_1 = 21
        pid_2 = 21
        pid_3 = 21 
        pid_4 = 21
        title = "$gg \\to gg(g)$"
        
    else:
        title = "not defined"   
    
    with open('Run.dat', 'w') as f:
        f.write("(run){\n\tEVENTS %d\n\tSHERPA_LDADD SherpaMain;\n\tSCALES VAR{sqr(91.18)}\n}(run)" % (0))
        f.write("\n\n(beam){\n\tBEAM_1  %d; BEAM_ENERGY_1  %.1f;\n\tBEAM_2  %d; BEAM_ENERGY_2  %.1f;\n}(beam)" % (pid_1, E, pid_2, E))
        f.write("\n\n(isr){\n\tPDF_LIBRARY None;\n}(isr)")
        f.write("\n\n(processes){\n\tProcess %d %d -> %d %d %d;\n\tME_Generator Amegic;\n\tOrder (3,0);\n\tEnd process;\n}(processes)" % (pid_1, pid_2, pid_3, pid_4, 21))
        
    return title;

def minkprod(a, b):
    metric = [[1., 0., 0., 0.],[0.,-1., 0., 0.],[0., 0., -1., 0.], [0., 0., 0., -1.]]
    return np.matmul(a, np.matmul(metric, b));

def calcTrace (process, p_1, p_2, p_3, p_4, p_s):
    N_C = 3.
    C_F = 4./3.
    C_A = 3.
    hard_p = np.array([p_1, p_2, p_3, p_4])
    
    s = minkprod(np.add(p_1, p_2), np.add(p_1, p_2))
    t = minkprod(np.subtract(p_1, p_3), np.subtract(p_1, p_3))
    u = minkprod(np.subtract(p_1, p_4), np.subtract(p_1, p_4))
    
    # qqbar -> qqbar
    if process[0] == 1:
        # qqbar -> qqbar
        if process[1] == 1:
            chi_1 = (pow(t,2)+pow(u,2))/pow(s,2)
            chi_2 = N_C*pow(u,2)/(s*t)-chi_1
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1/pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2)-2/N_C*pow(u,2)/(s*t)
                    
        # qqbar -> q'qbar'
        elif process[1] == 2:
            chi_1 = (pow(t,2)+pow(u,2))/pow(s,2)
            chi_2 = -chi_1
            chi_3 = 1/pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2)
                    
        # qqbar' -> qqbar'
        elif process[1] == 3:
            chi_1 = 0.
            chi_2 = 0.
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)
    
        H = 2./pow(N_C,2)*np.array([[pow(C_F,2)/pow(N_C,2)*chi_1, C_F/pow(N_C, 2)*chi_2],
                    [C_F/pow(N_C,2)*chi_2, chi_3]])
        M = np.array([[pow(N_C,2), 0.],
                    [0., 1./4.*(pow(N_C,2)-1)]])
    
        C_12 = np.array([[C_F, -C_F/(2.*N_C)+1./4.*(1.-1./pow(N_C,2))],
                        [0., -1./(2.*N_C)]])
        C_13 = np.array([[0., 1./4.*(1.-1./pow(N_C,2))],
                        [1., -1./N_C]])
        C_14 = np.array([[0., C_F/(2.*N_C)],
                        [1., -1./(2.*N_C)+C_F]])
        C_34 = C_12
        C_24 = C_13
        C_23 = C_14
        Gamma = np.zeros((2,2))
        
    # qq -> qq    
    elif process[0] == 2:
        # qq -> qq
        if process[1] == 1:
            chi_1 = (pow(t,2)+pow(s,2))/pow(u,2)
            chi_2 = N_C*pow(s,2)/(t*u)-chi_1
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1/pow(N_C,2)*(pow(s,2)+pow(t,2))/pow(u,2)-2/N_C*pow(s,2)/(t*u)
        
        # qq' -> qq'
        elif process[1] == 2:
            chi_1 = 0.
            chi_2 = 0.
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)
            
        H = np.array([[2.*pow(C_F,2)/pow(N_C, 4)*chi_1, 2.*C_F/pow(N_C, 4)*chi_2],
                    [2.*C_F/pow(N_C, 4)*chi_2, 2./pow(N_C, 2)*chi_3]])
        M = np.array([[pow(N_C,2), 0.],
                    [0., 1./4.*(pow(N_C,2)-1)]])
        
        C_12 = np.array([[0., (pow(N_C,2)-1)/pow(N_C,2)],[4., -4./n]])
        C_13 = np.array([[2.*(pow(N_C,2)-1)/N_C, -(N_C+1)/pow(N_C,2)],[0., -2./N_C]])
        C_14 = np.array([[0., (pow(N_C,2)-1)/pow(N_C,2)],[4., 2.*(pow(N_C,2)-2)/N_C]])
        C_34 = C_12
        C_24 = C_13
        C_23 = C_14
        Gamma = np.zeros((2,2))
        
    # qg -> qg
    elif process[0] == 3:
        chi_1 = -(pow(s,2)+pow(u,2))/(s*u)
        chi_2 = (1+2*u/t)*chi_1
        chi_3 = (1-4*s*u/pow(t,2))*chi_1
                    
        prfct = 1./(2.*N_C*(pow(N_C,2)-1))
        H = np.array([[prfct/pow(N_C,2)*chi_1, prfct/N_C*chi_1, prfct/N_C*chi_2],
                [prfct/N_C*chi_1, prfct*chi_1, prfct*chi_2],
                [prfct/N_C*chi_2, prfct*chi_2, prfct*chi_3]])
        M = np.array([[C_F*2*pow(N_C,2),0.,0.],
                [0., C_F*(pow(N_C,2)-4),0.],
                [0.,0.,C_F*pow(N_C,2)]])
        
        C_12 = np.array([[0., 0., -2j],[0., 1.j*N_C/2., -1.j*N_C/2.],[-1.j, -1.j*(pow(N_C,2)-4)/(2.*N_C), -1.j*N_C/2.]])
        C_13 = np.array([[2.*(pow(N_C, 2)-1)/N_C, 0., 0.],[0., -2./N_C, 0.],[0., 0., -2./N_C]])
        C_14 = np.array([[0., 0., -2.j],[0., -N_C/2., -1.j*N_C/2.],[-1.j, -1.j*(pow(N_C, 2)-4)/(2.*N_C), -1.j*N_C/2.]])
        C_24 = np.array([[N_C, 0., 0.],[0., N_C/2., 0.],[0., 0., N_C/2.]])
        C_34 = C_12
        C_23 = C_14
        Gamma = np.zeros((3,3))
        
    # qqbar -> gg and gg -> qqbar
    elif process[0] == 4:
        #qqbar -> gg
        if process[1] == 1:
            prfct = pow(N_C,2)
        
        # gg -> qqbar
        elif process[1] == 2:
            prfct = pow(pow(N_C,2)-1,2)
            
        chi_1 = (pow(t,2)+pow(u,2))/(t*u)
        chi_2 = (1+2*u/s)*chi_1
        chi_3 = (1-4*t*u/pow(s,2))*chi_1
                    
        H = np.array([[prfct/pow(N_C,2)*chi_1, prfct/N_C*chi_1, prfct/N_C*chi_2],
                [prfct/N_C*chi_1, prfct*chi_1, prfct*chi_2],
                [prfct/N_C*chi_2, prfct*chi_2, prfct*chi_3]])
        M = np.array([[C_F*2*pow(N_C,2),0.,0.],
                [0., C_F*(pow(N_C,2)-4),0.],
                [0.,0.,C_F*pow(N_C,2)]])
        Gamma = np.zeros((3,3))
        
    # gg -> gg
    elif process[0] == 5:
        chi_1 = 1-t*u/pow(s,2)-s*t/pow(u,2)+pow(t,2)/(s*u)
        chi_2 = s*t/pow(u,2)-t*u/pow(s,2)+pow(u,2)/(s*t)-pow(s,2)/(t*u)
        chi_3 = 27/4 - 9*(s*u/pow(t,2) +1/4*t*u/pow(s,2)+1/4*s*t/pow(u,2)) + 9/2*(pow(u,2)/(s*t)+pow(s,2)/(t*u)-1/2*pow(t,2)/(s*u))
        
        H = np.array([[9.*chi_1, 9./2.*chi_1, 9./2.*chi_2, 0., -3.*chi_1],
            [9./2.*chi_1, 9./4.*chi_1, 9./4.*chi_2, 0., -3./2.*chi_1],
            [9./2.*chi_2, 9./4.*chi_2, chi_3, 0., -3./2.*chi_2],
            [0., 0., 0., 0., 0.],
            [-3*chi_1, -3./2.*chi_1, -3./2.*chi_2, 0., chi_1]])
        M = np.array([
            [1.,0.,0.,0.,0.],
            [0.,8.,0.,0.,0.],
            [0.,0.,8.,0.,0.],
            [0.,0.,0.,20.,0.],
            [0.,0.,0.,0.,27.]
            ])
        Gamma = np.zeros((5,5))
        
    for i in range(0,4):
        for j in range(i+1, 4):
            if (i == 0 and j == 1):
                C = C_12
            elif (i == 2 and j == 3):
                C = C_34
            elif (i == 0 and j == 2): 
                C = C_13
            elif (i == 1 and j == 3):
                C = C_24
            elif (i == 0 and j == 3): 
                C = C_14
            elif (i == 1 and j == 2):
                C = C_23
            
            Gamma = np.add(Gamma, minkprod(hard_p[i], hard_p[j])/(minkprod(hard_p[i], p_s)*minkprod(hard_p[j], p_s))*C)
    
    trace = np.trace(np.matmul(np.matmul(H, M), Gamma))
    
    return trace;

alpha_s = 0.118
E = 1.e3
phi_43 = m.pi * 3./2.
phi_44 = m.pi * 5./2.
phi_H = m.pi/10.
phi_s = m.pi/7.
R_s = np.array([])
lambda_s = np.array([1.e-6, 5.e-6, 1.e-5, 5.e-5, 1.e-4, 5.e-4, 1.e-3, 5e-3, 1.e-2, 5.e-2, 1.e-1,])

process = [3,1]

plot_title = writeRunCard(process, E)
    
for l_s in lambda_s:
    k_s = 2.*E*l_s
    p_1 = [E,0.,0.,E]
    p_2 = [E,0.,0.,-E]    
    p_3 = [E,E*m.cos(phi_43+phi_H),E*m.sin(phi_43+phi_H),0]
    p_4 = [E,E*m.cos(phi_44+phi_H),E*m.sin(phi_44+phi_H),0]
    p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]
    
    MatEl = sh.calcMatEl(p_1, p_2, p_3, p_4, p_s)
    trace = calcTrace(process, p_1, p_2, p_3, p_4, p_s)
    R_s = np.append(R_s, [alpha_s/m.pi*trace/MatEl])
        
print "lambda_s\tR_s"
for i in range(0, lambda_s.size):
    print lambda_s[i], "\t\t", R_s[i]
        
#Plotting
plt.title(plot_title + '$\mathrm{~at~} E = %.1f \mathrm{~GeV}$' %(E))
plt.xlabel('$\lambda_s$')
plt.ylabel('$R_s$')
plt.xscale('log')
plt.plot(lambda_s, R_s)
plt.show()

    
    
