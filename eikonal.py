import math as m
import cmath as cm
import numpy as np
import interface as sh

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
    
    with open('Run.dat', 'w') as f:
        f.write("(run){\n\tEVENTS %d\n\tSHERPA_LDADD SherpaMain;\n}(run)" % (0))
        f.write("\n\n(beam){\n\tBEAM_1  %d; BEAM_ENERGY_1  %.1f;\n\tBEAM_2  %d; BEAM_ENERGY_2  %.1f;\n}(beam)" % (pid_1, E, pid_2, E))
        f.write("\n\n(isr){\n\tPDF_LIBRARY None;\n}(isr)")
        f.write("\n\n(processes){\n\tProcess %d %d -> %d %d %d;\n\tME_Generator Amegic;\n\tEnd process;\n}(processes)" % (pid_1, pid_2, pid_3, pid_4, 21))
        
    return;

def calcTrace (p_1, p_2, p_3, p_4, process):
    eps = 1e-10
    N_C = 3.
    C_F = 4./3.
    C_A = 3.
    
    s = (p_1[0]+p_2[0])**2-(p_1[1]+p_2[1])**2-(p_1[2]+p_2[2])**2-(p_1[3]+p_2[3])**2
    t = (p_1[0]-p_3[0])**2-(p_1[1]-p_3[1])**2-(p_1[2]-p_3[2])**2-(p_1[3]-p_3[3])**2
    u = (p_1[0]-p_4[0])**2-(p_1[1]-p_4[1])**2-(p_1[2]-p_4[2])**2-(p_1[3]-p_4[3])**2
    
    T = cm.log(-t/s)+m.pi*1j
    U = cm.log(-u/s)+m.pi*1j
    
    chi_1 = 0.
    chi_2 = 0.
    chi_3 = 0.
    prfct = 0.
    R_s = 0.+0.j
    
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
        
        H = np.array([[2*pow(C_F,2)/pow(N_C, 4)*chi_1, 2*C_F/pow(N_C, 4)*chi_2],
                [2*C_F/pow(N_C, 4)*chi_2, 2/pow(N_C, 2)*chi_3]])
        M = np.array([[pow(N_C,2), 0],
                [0, 1/4*(pow(N_C,2)-1)]])
        Gamma = np.array([[2*C_F*T, -C_F/N_C*U],
                [-2*U, -1/N_C*(T-2*U)]])
    
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
        Gamma = np.array([[2*C_F*T, C_F/N_C*U],
                [2*U, -1/N_C*(T-2*U)]])
        
    # qg -> qg
    elif process[0] == 3:
        chi_1 = -(pow(s,2)+pow(u,2))/(s*u)
        chi_2 = (1+2*u/t)*chi_1
        chi_3 = (1-4*s*u/pow(t,2))*chi_1
                    
        prfct = 1./(2.*N_C*(pow(N_C,2)-1))
        H = np.array([[prfct/pow(N_C,2)*chi_1, prfct/N_C*chi_1, prfct/N_C*chi_2],
                [prfct/N_C*chi_1, prfct*chi_1, prfct*chi_2],
                [prfct/N_C*chi_2, prfct*chi_2, prfct*chi_3]])
        Gamma = np.array([[(C_F+C_A)*T, 0., U],
                [0., C_F*T+C_A/2*U, C_A/2*U],
                [2*U, (pow(N_C,2)-4)/(2*N_C)*U, C_F*T+C_A/2*U]])
        M = np.array([[C_F*2*pow(N_C,2),0.,0.],
                [0., C_F*(pow(N_C,2)-4),0.],
                [0.,0.,C_F*pow(N_C,2)]])
        
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
        Gamma = np.array([[0., 0., U-T],
                [0., C_A/2*(T+U), C_A/2.*(U-T)],
                [2.*(U-T), (pow(N_C,2)-4.)/(2*N_C)*(U-T), C_A/2.*(T+U)]])
        M = np.array([[C_F*2*pow(N_C,2),0.,0.],
                [0., C_F*(pow(N_C,2)-4),0.],
                [0.,0.,C_F*pow(N_C,2)]])
        
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
        Gamma = np.array([
            [6*T, 0., -6.*U, 0., 0.],
            [0., 3.*T+3./2.*U, -3./2.*U, -3.*U, 0.],
            [-3./4.*U, -3/2*U, 3.*T+3./2.*U, 0., -9./4.*U],
            [0., -6./5.*U, 0., 3.*U, -9./5.*U],
            [0., 0., -2./3.*U, -4/3*U, -2.*T+4.*U]])
        
    trace = np.trace(np.matmul(H,np.matmul(np.conj(Gamma.transpose()),np.matmul(M, Gamma))))
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
process = [5, 1]
maxrange = 10

writeRunCard(process, E)

for i in range (0,maxrange):
    lambda_s = np.append(lambda_s, 1*pow(10, -i))
    k_s = 2*E*lambda_s[i]
    p_1 = [E,0.,0.,E]
    p_2 = [E,0.,0.,-E]    
    p_3 = [E,E*m.cos(phi_43+phi_H),E*m.sin(phi_43+phi_H),0]
    p_4 = [E,E*m.cos(phi_44+phi_H),E*m.sin(phi_44+phi_H),0]
    p_s = [k_s,k_s*m.cos(phi_s),k_s*m.sin(phi_s),0.]

    MatEl = sh.calcMatEl(p_1, p_2, p_3, p_4, p_s)
    trace = calcTrace(p_1, p_2, p_3, p_4, process)
    R_s = np.append(R_s, [alpha_s/m.pi * trace/MatEl])
    
print "lambda_s\tR_s"
for i in range(0,maxrange):
    print lambda_s[i], "\t\t", R_s[i]
    
    
