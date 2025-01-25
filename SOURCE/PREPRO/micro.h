#ifndef __MICRO__
#define __MICRO__

typedef double R;

void search_cont(int & N_SEG, bool * LIST_N, int H1, int H2, int IN, int OUT);
void treat_cont(int N_SEG, bool * LIST_N, int ** LIST_Q, int * LIST_V, int H1, int H2, int IN, int OUT);
void search_surf(int & N_QUA, bool * LIST_N, int H1,int H2,int H3, int IN, int OUT);
void treat_surf(int N_QUA, bool * LIST_N, int ** LIST_Q, int * LIST_V, int H1,int H2,int H3, int IN, int OUT);
void treatbm2(int H1, int H2, bool * LIST_N);
void treatbp2(int H1, int H2, bool * LIST_N);
void treatbm2b(int H1, int H2, bool * LIST_N);
void treatbp2b(int H1, int H2, bool * LIST_N);
void detect_agr2(int H1, int H2, int * LIST_A, bool * LIST_N, int N_TOT);
void treat_per2(int H1, int H2, bool * LIST_N);
void treat_per3(int H1, int H2, int H3, bool * LIST_N);
void eval_fv(R & PR_INC, int N_TOT, bool * LIST_N);

#endif
