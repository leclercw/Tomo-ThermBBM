#ifndef _GENESIS_
#define _GENESIS_

typedef double R;

void gener2_disk_alea(int NB_DIS,R R_DIS,int N_TOT,int H1, int H2, bool * LIST_N, R * LIST_X, R * LIST_Y, R * LIST_R);
bool bincl(R R_CYL,R AA,R BB,R CC, R X1, R Y1, R Z1, R X2, R Y2, R Z2, R XX, R YY, R ZZ);

void gener3_cyl_alea_curve_int_per(int THETAM,int PHIM,int nC, R facf, R eche, int N_TOT,int V_TOT,int H1,int H2,int H3, bool * LIST_N, R* LIST_PS, R* LIST_GA, R* LIST_PH);

void gener3_sphere_alea_per(int nS, R eche, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N, R* LIST_R, R* LIST_X, R* LIST_Y, R* LIST_Z);
void gener3_sphere_alea(int nS, R eche, R buffer, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N, R* LIST_R, R* LIST_X, R* LIST_Y, R* LIST_Z);
void gener3_cyl_inf_per(int nC, R eche, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir);
void gener3_cyl_inf(int nC, R eche, R buffer, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir);

void gener3_cyl_align(int DENSI,int nC, R facf, R eche, R buffer, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir);
void gener3_cyl_align_per(int DENSI,int nC, R facf, R eche, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir);

void gener3_cyl_alea(int DENSI,int nC, R facf, R eche, R buffer,  int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,R * LIST_PS,R * LIST_GA,R * LIST_PH);
void gener3_cyl_alea_per(int DENSI,int nC, R facf, R eche, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,R * LIST_PS,R * LIST_GA,R * LIST_PH);

#endif
