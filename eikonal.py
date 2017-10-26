import math as m
import cmath as cm
import numpy as np
import matplotlib.pyplot as plt
import interface as sh
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def writeRunCard(E):
    pid_1 = 1
    pid_2 = -1
    pid_3 = 1
    pid_4 = -1
    
    with open('Run.dat', 'w') as f:
        f.write("(run){\n\tEVENTS %d\n\tSHERPA_LDADD SherpaMain;\n}(run)" % (0))
        f.write("\n\n(beam){\n\tBEAM_1  %d; BEAM_ENERGY_1  %.1f;\n\tBEAM_2  %d; BEAM_ENERGY_2  %.1f;\n}(beam)" % (pid_1, E, pid_2, E))
        f.write("\n\n(isr){\n\tPDF_LIBRARY None;\n}(isr)")
        f.write("\n\n(processes){\n\tProcess %d %d -> %d %d %d;\n\tME_Generator Amegic;\n\tEnd process;\n}(processes)" % (pid_1, pid_2, pid_3, pid_4, 21))
        
    return "$qq \\to qq(g)$";

def minkprod(a, b):
    metric = [[1., 0., 0., 0.],[0.,-1., 0., 0.],[0., 0., -1., 0.], [0., 0., 0., -1.]]
    return np.matmul(a, np.matmul(metric, b));

def calcTrace (p_1, p_2, p_3, p_4, p_s):
    eps = 1e-10
    N_C = 3.
    C_F = 4./3.
    C_A = 3.
    
    hard_p = 1.e9*np.array([p_1, p_2, p_3, p_4]) # factor of 10^9 due to GeV
        
    s = minkprod(np.add(p_1, p_2), np.add(p_1, p_2))
    t = minkprod(np.subtract(p_1, p_3), np.subtract(p_1, p_3))
    u = minkprod(np.subtract(p_1, p_4), np.subtract(p_1, p_4))
    
    chi_1 = (pow(t,2)+pow(u,2))/pow(s,2)
    chi_2 = N_C*pow(u,2)/(s*t)-chi_1
    chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1./pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2)-2./N_C*pow(u,2)/(s*t)
    
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

    Gamma = np.zeros((2,2))
    for i in range(0,4):
        for j in range(i+1, 4):
            if (i == 0 and j == 1) or (i == 2 and j == 3):
                C = C_12
            elif (i == 0 and j == 2) or (i == 1 and j == 3):
                C = C_13
            elif (i == 0 and j == 3) or (i == 1 and j == 2):
                C = C_14
            else:
                C = np.zeros((2,2))
            
            Gamma += minkprod(hard_p[i], hard_p[j])/(minkprod(hard_p[i], p_s)*minkprod(hard_p[j], p_s))*C
    
    trace = np.trace(np.matmul(np.matmul(H, M), Gamma))
    #trace = np.trace(np.matmul(H,np.matmul(np.conj(Gamma.transpose()),np.matmul(M, Gamma)))) 
    if trace.imag < eps:
        trace = trace.real
    
    return trace;

alpha_s = 0.118
E = 7.e3
phi_43 = m.pi * 3./2.
phi_44 = m.pi * 5./2.
phi_H = m.pi/10.
phi_s = m.pi/7.
R_s = np.array([])
lambda_s = np.array([1.e-6, 5.e-6, 1.e-5, 5.e-5, 1.e-4, 5.e-4, 1.e-3, 5e-3, 1.e-2, 5.e-2, 1.e-1, 5e-1, 1.])

plot_title = writeRunCard(E)

for l_s in lambda_s:
    k_s = 2.*E*l_s
    p_1 = [E,0.,0.,E]
    p_2 = [E,0.,0.,-E]    
    p_3 = [E,E*m.cos(phi_43+phi_H),E*m.sin(phi_43+phi_H),0]
    p_4 = [E,E*m.cos(phi_44+phi_H),E*m.sin(phi_44+phi_H),0]
    p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]
    
    MatEl = sh.calcMatEl(p_1, p_2, p_3, p_4, p_s)
    trace = calcTrace(p_1, p_2, p_3, p_4, p_s)
    R_s = np.append(R_s, [alpha_s/m.pi* trace/MatEl])
        
print "lambda_s\tR_s"
for i in range(0, lambda_s.size):
    print lambda_s[i], "\t\t", R_s[i]
        
#Plotting
plt.title(plot_title)
plt.xlabel('$\lambda_s$')
plt.ylabel('$R_s$')
plt.xscale('log')
plt.plot(lambda_s, R_s)
plt.show()

    
    
