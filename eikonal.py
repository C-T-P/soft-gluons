import math as m
import cmath as cm
import numpy as np
import interface as sh

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
        
    return;

def minkprod(a, b):
    metric = [[1., 0., 0., 0.],[0.,-1., 0., 0.],[0., 0., -1., 0.], [0., 0., 0., -1.]]
    return np.matmul(a, np.matmul(metric, b));

def calcTrace (p_1, p_2, p_3, p_4, p_s):
    eps = 1e-10
    N_C = 3.
    C_F = 4./3.
    C_A = 3.
    alpha_s = 0.118
    
    s = minkprod(np.add(p_1, p_2), np.add(p_1, p_2))
    t = minkprod(np.subtract(p_1, p_3), np.subtract(p_1, p_3))
    u = minkprod(np.subtract(p_1, p_4), np.subtract(p_1, p_4))
    
    #s = (p_1[0]+p_2[0])**2-(p_1[1]+p_2[1])**2-(p_1[2]+p_2[2])**2-(p_1[3]+p_2[3])**2
    #t = (p_1[0]-p_3[0])**2-(p_1[1]-p_3[1])**2-(p_1[2]-p_3[2])**2-(p_1[3]-p_3[3])**2
    #u = (p_1[0]-p_4[0])**2-(p_1[1]-p_4[1])**2-(p_1[2]-p_4[2])**2-(p_1[3]-p_4[3])**2
    #print s, t, u
    
    # qqbar -> qqbar using the phd thesis
    chi_1 = 2.*pow(C_F,2)/pow(N_C,4)*(pow(t,2)+pow(u,2))/pow(s,2)
    chi_2 = 2.*C_F/pow(N_C,3)*(pow(u,2)/(s*t)-1./N_C*chi_1)
    chi_3 = 1./pow(N_C,2)*(2.*(pow(s,2)+pow(u,2))/pow(t,2)+2./pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2)-4./N_C*pow(u,2)/(s*t))
        
    H = pow(alpha_s, 2)*np.array([[chi_1, chi_2], [chi_2, chi_3]])
    
    #H = np.array([[2*pow(C_F,2)/pow(N_C, 4)*chi_1, 2*C_F/pow(N_C, 4)*chi_2],
                #[2*C_F/pow(N_C, 4)*chi_2, 2/pow(N_C, 2)*chi_3]])
    M = np.array([[pow(N_C,2), 0.],
                [0., 1./4.*(pow(N_C,2)-1)]])
        
    Gamma = (minkprod(p_1, p_2)/(minkprod(p_1, p_s)*minkprod(p_2, p_s)) + minkprod(p_3, p_4)/(minkprod(p_3, p_s)*minkprod(p_4, p_s))) *np.array([[C_F, 1./2],[0, -1./(2.*N_C)]]) + (minkprod(p_1, p_3)/(minkprod(p_1, p_s) + minkprod(p_3, p_s)) + minkprod(p_2, p_4)/(minkprod(p_2, p_s)*minkprod(p_4, p_s))) * np.array([[C_F, 1./2],[1./2., -1./(2.*N_C)]]) + ( minkprod(p_1, p_4)/(minkprod(p_1, p_s)*minkprod(p_4, p_s)) + minkprod(p_2, p_3)/(minkprod(p_2, p_s)*minkprod(p_3, p_s)) ) * np.array([[-1./(2.*N_C), 0.],[1./2., C_F]])
    
    print Gamma
    
    T1 = np.array([[1.0, 1./N_C],[0., 2.]])
    T2 = np.array([[1.0, -1./N_C],[0., 1./2.]])
    Gamma = np.matmul(np.matmul(T1, Gamma), T2)
    
    trace = np.trace(np.matmul(H,np.matmul(np.conj(Gamma.transpose()),np.matmul(M, Gamma)))) 
    print trace
    if trace.imag < eps:
        trace = trace.real
    
    return trace;

alpha_s = 0.118
E = 50.0
phi_43 = m.pi * 3./2.
phi_44 = m.pi * 5./2.
phi_H = m.pi/10.
phi_s = m.pi/7.
R_s = np.array([])
lambda_s = np.array([])
maxrange = 10

writeRunCard(E)

for i in range (0,maxrange):
    lambda_s = np.append(lambda_s, pow(10, -i))
    k_s = 2.*E*lambda_s[i]
    p_1 = [E,0.,0.,E]
    p_2 = [E,0.,0.,-E]    
    p_3 = [E,E*m.cos(phi_43+phi_H),E*m.sin(phi_43+phi_H),0]
    p_4 = [E,E*m.cos(phi_44+phi_H),E*m.sin(phi_44+phi_H),0]
    p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]

    MatEl = sh.calcMatEl(p_1, p_2, p_3, p_4, p_s)
    trace = calcTrace(p_1, p_2, p_3, p_4, p_s)
    R_s = np.append(R_s, [alpha_s/m.pi * trace/MatEl])
    
print "lambda_s\tR_s"
for i in range(0,maxrange):
    print lambda_s[i], "\t\t", R_s[i]
    
    
