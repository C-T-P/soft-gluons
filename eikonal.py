import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt
import interface as sh
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

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
    
        C_12 = np.array([[0., (pow(N_C,2)-1)/(4*pow(N_C,2))],
                        [1., (pow(N_C,2)-2)/(2.*N_C)]])
        C_13 = np.array([[(pow(N_C,2)-1)/(2.*N_C), 0.],
                        [0., -1./(2.*N_C)]])
        C_14 = -np.array([[0., (pow(N_C,2)-1)/(4.*pow(N_C,2))],
                        [1., -1./N_C]])
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
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1./pow(N_C,2)*(pow(s,2)+pow(t,2))/pow(u,2)-2./N_C*pow(s,2)/(t*u)
        
        # qq' -> qq'
        elif process[1] == 2:
            chi_1 = 0.
            chi_2 = 0.
            chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)
            
        H = np.array([[2.*pow(C_F,2)/pow(N_C, 4)*chi_1, 2.*C_F/pow(N_C, 4)*chi_2],
                    [2.*C_F/pow(N_C, 4)*chi_2, 2./pow(N_C, 2)*chi_3]])
        M = np.array([[pow(N_C,2), 0.],
                    [0., 1./4.*(pow(N_C,2)-1)]])
        
        C_12 = -np.array([[0., (pow(N_C,2)-1)/(4*pow(N_C,2))],[1., -1./N_C]])
        C_13 = np.array([[(pow(N_C,2)-1)/(2.*N_C), 0.],[0., -1./(2.*N_C)]])
        C_14 = np.array([[0., (pow(N_C,2)-1)/(4.*pow(N_C,2))],[1., (pow(N_C,2)-2)/(2.*N_C)]])
        C_34 = C_12
        C_24 = C_13
        C_23 = C_14
        Gamma = np.zeros((2,2))
        
    # qg -> qg
    elif process[0] == 3:
        chi_1 = -(pow(s,2)+pow(u,2))/(s*u)
        chi_2 = (1.+2.*u/t)*chi_1
        chi_3 = (1.-4.*s*u/pow(t,2))*chi_1
                    
        H = 1./(2.*N_C*(pow(N_C,2)-1))*np.array([[1./pow(N_C,2)*chi_1, 1./N_C*chi_1, 1./N_C*chi_2],
                [1./N_C*chi_1, chi_1, chi_2],
                [1./N_C*chi_2, chi_2, chi_3]])
        M = C_F*np.array([[2*pow(N_C,2),0.,0.],
                [0., (pow(N_C,2)-4),0.],
                [0.,0.,pow(N_C,2)]])
        
        C_12 = np.array([[0., 0., -1./2.],[0., N_C/4., -N_C/4.],[-1., -(pow(N_C,2)-4)/(4.*N_C), N_C/4.]])
        C_13 = np.array([[(pow(N_C, 2)-1)/(2.*N_C), 0., 0.],[0., -1./(2.*N_C), 0.],[0., 0., -1./(2.*N_C)]])
        C_14 = np.array([[0., 0., 1./2.],[0., N_C/4., N_C/4.],[1., (pow(N_C, 2)-4)/(4.*N_C), N_C/4.]])
        C_24 = np.array([[N_C, 0., 0.],[0., N_C/2., 0.],[0., 0., N_C/2.]])
        C_34 = C_12
        C_23 = C_14
        Gamma = np.zeros((3,3))
        
    # qqbar -> gg and gg -> qqbar
    elif process[0] == 4:
        #qqbar -> gg
        if process[1] == 1:
            prfct = 1./(6.*pow(N_C,2)) # extra factor of 1/3 for identical final state ptcls
            
            C_12 = np.array([[(pow(N_C,2)-1)/(2.*N_C), 0., 0.],[0., -1./(2.*N_C), 0.], [0., 0., -1./(2.*N_C)]])
            C_13 = -np.array([[0., 0., 1./2.],[0., -N_C/4., N_C/4.],[1., (pow(N_C,2)-4)/(4.*N_C), -N_C/4]])
            C_14 = np.array([[0., 0., 1./2.],[0., N_C/4., N_C/4.],[1., (pow(N_C,2)-4)/(4.*N_C), N_C/4.]])
            C_34 = np.array([[N_C, 0., 0.,],[0., N_C/2., 0.],[0., 0., N_C/2.]])
            C_24 = C_13
            C_23 = C_14
        
        # gg -> qqbar
        elif process[1] == 2:
            prfct = 1./(2.*pow(pow(N_C,2)-1,2))
            
            C_12 = np.array([[N_C, 0., 0.,],[0., N_C/2., 0.],[0., 0., N_C/2.]])
            C_34 = np.array([[(pow(N_C,2)-1)/(2.*N_C), 0., 0.],[0., -1./(2.*N_C), 0.], [0., 0., -1./(2.*N_C)]])
            C_13 = -np.array([[0., 0., -1./2.],[0., -N_C/4., -N_C/4.],[-1., -(pow(N_C,2)-4)/(4.*N_C), -N_C/4]])
            C_14 = -np.array([[0., 0., 1./2.],[0., -N_C/4., N_C/4.],[1., (pow(N_C,2)-4)/(4.*N_C), -N_C/4.]])
            C_24 = C_13
            C_23 = C_14
            
        chi_1 = (pow(t,2)+pow(u,2))/(t*u)
        chi_2 = (1.+2.*u/s)*chi_1
        chi_3 = (1.-4.*t*u/pow(s,2))*chi_1
                    
        H = prfct*np.array([[1./pow(N_C,2)*chi_1, 1./N_C*chi_1, 1./N_C*chi_2],
                [1./N_C*chi_1, chi_1, chi_2],
                [1./N_C*chi_2, chi_2, chi_3]])
        M = C_F*np.array([[2*pow(N_C,2),0.,0.],
                [0., pow(N_C,2)-4,0.],
                [0.,0.,pow(N_C,2)]])
        
        
        Gamma = np.zeros((3,3))
        
    # gg -> gg
    elif process[0] == 5:
        chi_1 = 1.-t*u/pow(s,2)-s*t/pow(u,2)+pow(t,2)/(s*u)
        chi_2 = s*t/pow(u,2)-t*u/pow(s,2)+pow(u,2)/(s*t)-pow(s,2)/(t*u)
        chi_3 = 27./4 - 9.*(s*u/pow(t,2) +1./4.*t*u/pow(s,2)+1./4.*s*t/pow(u,2)) + 9./2.*(pow(u,2)/(s*t)+pow(s,2)/(t*u)-1./2.*pow(t,2)/(s*u))
        
        H = 1./3.*1/16.*np.array([[9.*chi_1, 9./2.*chi_1, 9./2.*chi_2, 0., -3.*chi_1],
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
        
        C_12 = np.array([
            [0, 0, 3., 0, 0],
            [0, 3./4., 3./4., -3./2., 0.],
            [3./8., 3./4., 3./4., 0., 9./8.],
            [0, 3./5., 0, 3./2., 9./10.],
            [0, 0, 1./3., 2./3., 2.]
            ])
        C_13 = np.array([
            [3., 0, 0, 0, 0],
            [0, 3./2., 0, 0, 0],
            [0, 0, 3./2., 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, -1.]
            ])
        C_14= np.array([
            [0, 0, -3., 0, 0],
            [0, 3./4., -3./4., -3./2., 0],
            [-3./8., -3./4., 3./4., 0, -9./8.],
            [0, -3./5., 0, 3./2., -9./10.],
            [0, 0, -1./3., -2./3., 2.]
            ])
        C_34 = C_12
        C_24 = C_13
        C_23 = C_14
        
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
    #trace = np.trace(np.matmul(H, M))
    
    return trace;

# Definition of constants
alpha_s = 0.118
phi_43 = m.pi * 3./2.
phi_44 = m.pi * 5./2.
phi_H = m.pi/10.
phi_s = m.pi/7.

# Process specification
process = [1,1]
E = 7.e3
K = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])
lambda_s = np.logspace(-2, -1, 70)
R_s = np.zeros((K.size,lambda_s.size))

# calculate Matrix Elements and get title for plot in first variable
plot_title, MatEl = sh.calcMatEl(process, K, lambda_s, E)

# Plot setup
plt.title(plot_title + '$\mathrm{~at~} E_{\mathrm{Beam}} = %.1f \mathrm{~GeV}$' %(E))
plt.xlabel('$\lambda_s$')
plt.ylabel('$R_s$')
plt.xscale('log') 

for i in range(1, K.size):
    N = K[i]
    for j in range(0, lambda_s.size):
        k_s = 2.*E*lambda_s[j]
        p_1 = [E,0.,0.,E]
        p_2 = [E,0.,0.,-E]    
        p_3 = [E,E*m.cos(phi_43+N*phi_H),E*m.sin(phi_43+N*phi_H),0]
        p_4 = [E,E*m.cos(phi_44+N*phi_H),E*m.sin(phi_44+N*phi_H),0]
        p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]
        
        #Ntrace = calcTrace(process, p_1, p_2, p_3, p_4, p_s)
        #R_s[i,j] = Ntrace/MatEl[i,j]*(4*m.pi*alpha_s)**2
        trace = calcTrace(process, p_1, p_2, p_3, p_4, p_s)
        R_s[i,j] = trace/MatEl[i,j]*(4*m.pi*alpha_s)**3
        
    #print R_s[i]
    plt.plot(lambda_s, R_s[i], label="N = %i" %(N))    

plt.legend(bbox_to_anchor=(-0.01, 0), loc='lower left')   
plt.savefig("R_plot_E-"+str(E)+"_process-"+str(process[0])+str(process[1]), dpi=300, orientation='portrait', papertype = 'a4', format='pdf', transparent=False, frameon=True)
#plt.show()