#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>

#define N_C 3
#define C_F 4/3
#define C_A 3

void mmult(int dim, double complex mat1[dim][dim], double complex mat2[dim][dim], double complex matResult[dim][dim]);
void hermconj(int dim, double complex m[dim][dim], double complex mresult[dim][dim]);
double complex calc_R_s(double alpha_s, double E_cms, double theta_cms, int process, int subprocess);
double complex mtrace(int dim, double complex mat[dim][dim]);
void printmat(int dim, double complex m[dim][dim]);

int main(int argc, char **argv) {
    double alpha_s = 1e0;
    double theta_cms = M_PI/2; // in rad
    
    double E_cms;
    int process;
    int subprocess;
    if(argv[1] != NULL && argv[2] != NULL && argv[3] != NULL) {
        E_cms = atof(argv[1]);
        process = atoi(argv[2]);
        subprocess = atoi(argv[3]);
    }
    else {
        return 1;
    }
    
    double R_s;
    int delta = 1;
    printf("Process Nr %d.%d\n", process, subprocess);
    for(int i=1; i <= delta; i++) {
        R_s = calc_R_s(alpha_s, E_cms, i*theta_cms/delta, process, subprocess);
        printf("%f\t%f+i%f\n", i*theta_cms/delta, creal(R_s), cimag(R_s));
    }
    return 0;
}
double complex calc_R_s(double alpha_s, double E_cms, double theta_cms, int process, int subprocess) {
    double complex trace_mat_El = 1.0;
    
    double s = 4*pow(E_cms,2), t = -2*pow(E_cms,2)*(1-cos(theta_cms)), u = -2*pow(E_cms,2)*(1+cos(theta_cms));
    
    double complex T = clog(-t/s)+M_PI*I;
    double complex U = clog(-u/s)+M_PI*I;
    double chi_1;
    double chi_2;
    double chi_3;
    double prfct;
    double complex R_s;
        
    switch(process) {
        // qqbar -> qqbar type
        case 1:
            switch(subprocess) {
                // qqbar -> qqbar
                case 1:
                    chi_1 = (pow(t,2)+pow(u,2))/pow(s,2);
                    chi_2 = N_C*pow(u,2)/(s*t)-chi_1;
                    chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1/pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2)-2/N_C*pow(u,2)/(s*t);
                    break;
                    
                // qqbar -> q'qbar'
                case 2:
                    chi_1 = (pow(t,2)+pow(u,2))/pow(s,2);
                    chi_2 = -chi_1;
                    chi_3 = 1/pow(N_C,2)*(pow(t,2)+pow(u,2))/pow(s,2);
                    break;
                    
                // qqbar' -> qqbar'
                case 3:
                    chi_1 = 0;
                    chi_2 = 0;
                    chi_3 = (pow(s,2)+pow(u,2))/pow(t,2);
                    break;
            }
            
            double complex H_qqbar[2][2] = {
                {2*pow(C_F,2)/pow(N_C, 4)*chi_1, 2*C_F/pow(N_C, 4)*chi_2},
                {2*C_F/pow(N_C, 4)*chi_2, 2.0/pow(N_C, 2)*chi_3}
            };
            double complex M_qqbar[2][2] = {
                {pow(N_C,2), 0},
                {0, 1.0/4*(pow(N_C,2)-1)}
            };
            double complex Gamma_qqbar[2][2] = {
                {2*C_F*T, -C_F/N_C*U},
                {-2*U, -1.0/N_C*(T-2*U)}
            };
            double complex Gamma_qqbar_dagger[2][2];
            hermconj(2, Gamma_qqbar, Gamma_qqbar_dagger);
            
            mmult(2, Gamma_qqbar_dagger, M_qqbar, Gamma_qqbar_dagger);
            mmult(2, Gamma_qqbar, Gamma_qqbar_dagger, Gamma_qqbar);
            mmult(2, H_qqbar, Gamma_qqbar, H_qqbar);
            R_s = alpha_s/M_PI*mtrace(2, H_qqbar)/trace_mat_El;
            break;
            
        // qq -> qq
        case 2:
            switch(subprocess) {
                // qq -> qq
                case 1:
                    chi_1 = (pow(t,2)+pow(s,2))/pow(u,2);
                    chi_2 = N_C*pow(s,2)/(t*u)-chi_1;
                    chi_3 = (pow(s,2)+pow(u,2))/pow(t,2)+1/pow(N_C,2)*(pow(s,2)+pow(t,2))/pow(u,2)-2/N_C*pow(s,2)/(t*u);
                    break;
                    
                // qq' -> qq'
                case 2:
                    chi_1 = 0;
                    chi_2 = 0;
                    chi_3 = (pow(s,2)+pow(u,2))/pow(t,2);
                    break;
            }
            
            double complex H_qq[2][2] = {
                {2*pow(C_F,2)/pow(N_C, 4)*chi_1, 2*C_F/pow(N_C, 4)*chi_2},
                {2*C_F/pow(N_C, 4)*chi_2, 2.0/pow(N_C, 2)*chi_3}
            };
            double complex M_qq[2][2] = {
                {pow(N_C,2), 0},
                {0, 1.0/4*(pow(N_C,2)-1)}
            };
            double complex Gamma_qq[2][2] = {
                {2*C_F*T, C_F/N_C*U},
                {2*U, -1.0/N_C*(T-2*U)}
            };
            double complex Gamma_qq_dagger[2][2];
            hermconj(2, Gamma_qq, Gamma_qq_dagger);
            
            mmult(2, Gamma_qq_dagger, M_qq, Gamma_qq_dagger);
            mmult(2, Gamma_qq, Gamma_qq_dagger, Gamma_qq);
            mmult(2, H_qq, Gamma_qq, H_qq);
            R_s = alpha_s/M_PI*mtrace(2, H_qq)/trace_mat_El;
            break;
            
        // qg -> qg
        case 3:
            chi_1 = -(pow(s,2)+pow(u,2))/(s*u);
            chi_2 = (1+2*u/t)*chi_1;
            chi_3 = (1-4*s*u/pow(t,2))*chi_1;
                    
            prfct = 1.0/(2*N_C*(pow(N_C,2)-1));
            double complex H_qg[3][3] = {
                {prfct/pow(N_C,2)*chi_1, prfct/N_C*chi_1, prfct/N_C*chi_2},
                {prfct/N_C*chi_1, prfct*chi_1, prfct*chi_2},
                {prfct/N_C*chi_2, prfct*chi_2, prfct*chi_3}
            };
            double complex M_qg[3][3] = {
                {C_F*2*pow(N_C,2),0,0},
                {0, C_F*(pow(N_C,2)-4),0},
                {0,0,C_F*pow(N_C,2)}
            };
            double complex Gamma_qg[3][3] = {
                {(C_F+C_A)*T, 0, U},
                {0, C_F*T+C_A/2*U, C_A/2*U},
                {2*U, (pow(N_C,2)-4)/(2*N_C)*U, C_F*T+C_A/2*U}
            };
            double complex Gamma_qg_dagger[2][2];
            hermconj(3, Gamma_qg, Gamma_qg_dagger);
                
            mmult(3, Gamma_qg_dagger, M_qg, Gamma_qg_dagger);
            mmult(3, Gamma_qg_dagger, Gamma_qg, Gamma_qg);
            mmult(3, H_qg, Gamma_qg, H_qg);
            R_s = alpha_s/M_PI*mtrace(3, H_qg)/trace_mat_El;
            break;
            
        // qqbar -> gg and gg -> qqbar
        case 4:
            prfct = 0;
            switch(subprocess) {
                // qqbar -> gg
                case 1:
                    prfct = pow(N_C,2);
                    break;
                
                // gg -> qqbar
                case 2:
                    prfct = pow(pow(N_C,2)-1,2);
                    break;
            }
            chi_1 = (pow(t,2)+pow(u,2))/(t*u);
            chi_2 = (1+2*u/s)*chi_1;
            chi_3 = (1-4*t*u/pow(s,2))*chi_1;
                    
            double complex H_qqbargg[3][3] = {
                {prfct/pow(N_C,2)*chi_1, prfct/N_C*chi_1, prfct/N_C*chi_2},
                {prfct/N_C*chi_1, prfct*chi_1, prfct*chi_2},
                {prfct/N_C*chi_2, prfct*chi_2, prfct*chi_3}
            };
            double complex M_qqbargg[3][3] = {
                {C_F*2*pow(N_C,2),0,0},
                {0, C_F*(pow(N_C,2)-4),0},
                {0,0,C_F*pow(N_C,2)}
            };
            double complex Gamma_qqbargg[3][3] = {
                {0, 0, U-T},
                {0, C_A/2*(T+U), C_A/2*(U-T)},
                {2*(U-T), (pow(N_C,2)-4)/(2*N_C)*(U-T), C_A/2*(T+U)}
            };
            double complex Gamma_qqbargg_dagger[3][3];
            hermconj(3, Gamma_qqbargg, Gamma_qqbargg_dagger);
            
            mmult(3, Gamma_qqbargg_dagger, M_qqbargg, Gamma_qqbargg_dagger);
            mmult(3, Gamma_qqbargg_dagger, Gamma_qqbargg, Gamma_qqbargg);
            mmult(3, H_qqbargg, Gamma_qqbargg, H_qqbargg);
            R_s = alpha_s/M_PI*mtrace(3, H_qqbargg)/trace_mat_El;
            break;
            
        // gg -> gg
        case 5:
            chi_1 = 1-t*u/pow(s,2)-s*t/pow(u,2)+pow(t,2)/(s*u);
            chi_2 = s*t/pow(u,2)-t*u/pow(s,2)+pow(u,2)/(s*t)-pow(s,2)/(t*u);
            chi_3 = 27/4 - 9*(s*u/pow(t,2) +1/4*t*u/pow(s,2)+1/4*s*t/pow(u,2)) + 9/2*(pow(u,2)/(s*t)+pow(s,2)/(t*u)-1/2*pow(t,2)/(s*u));
            
            double complex H_gg[5][5] = {
                {9*chi_1, 9.0/2.0*chi_1, 9.0/2*chi_2, 0.0, -3*chi_1},
                {9.0/2*chi_1, 9.0/4*chi_1, 9.0/4*chi_2, 0, -3.0/2*chi_1},
                {9.0/2*chi_2, 9.0/4*chi_2, chi_3, 0, -3.0/2*chi_2},
                {0, 0, 0, 0, 0},
                {-3*chi_1, -3/2*chi_1, -3/2*chi_2, 0, chi_1}
            };
            double complex Gamma_gg[5][5] = {
                {6*T, 0, -6*U, 0, 0},
                {0, 3*T+3.0/2*U, -3.0/2*U, -3*U, 0},
                {-3.0/4*U, -3.0/2*U, 3*T+3.0/2*U, 0, -9.0/4*U},
                {0, -6.0/5*U, 0, 3*U, -9.0/5*U},
                {0, 0, -2.0/3*U, -4.0/3*U, -2*T+4*U},
            };
            double complex Gamma_gg_dagger[5][5];
            hermconj(5, Gamma_gg, Gamma_gg_dagger);
                
            double complex M_gg[5][5] = {
                {1,0,0,0,0},
                {0,8,0,0,0},
                {0,0,8,0,0},
                {0,0,0,20,0},
                {0,0,0,0,27}
            };
            mmult(5, Gamma_gg_dagger, M_gg, Gamma_gg_dagger);
            mmult(5, Gamma_gg_dagger, Gamma_gg, Gamma_gg);
            mmult(5, H_gg, Gamma_gg, H_gg);
            R_s = alpha_s/M_PI*mtrace(5, H_gg)/trace_mat_El;
            break;
            
        default:
            R_s = NAN;
            break;
    }
    
    return R_s;
}
void mmult(int dim, double complex mat1[dim][dim], double complex mat2[dim][dim], double complex mresult[dim][dim]) {
    double complex dummy[dim][dim];
    
    for (int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            dummy[i][j] = 0;
    
    for (int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            for(int k = 0; k < dim; k++)
                dummy[i][j] += mat1[i][k]*mat2[k][j];

    for (int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            mresult[i][j] = dummy[i][j];
}
void hermconj(int dim, double complex m[dim][dim], double complex mresult[dim][dim]) {
    for(int i=0; i<dim; i++)
        for(int j=0; j<dim; j++)
            mresult[i][j] = conj(m[j][i]);
} 
double complex mtrace(int dim, double complex mat[dim][dim]) {
    double complex trace = 0;
    for (int i = 0; i < dim; i++)
        trace += mat[i][i];
    
    return trace;
}
void printmat(int dim, double complex m[dim][dim]) {
    for(int i =0; i<dim; i++) {
        for(int j=0; j<dim; j++) {
            printf("%f+i%f\t", creal(m[i][j]), cimag(m[i][j]));
        }
        printf("\n");
    }
    printf("\n");
}
