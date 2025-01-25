#ifndef ___VTK___
#define ___VTK___

typedef double R;

void read_H(int & H_TOT,const char * filename);
void readmap3(int H1, int & NB_SPH, R & L_VER, R & R_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, const char * filename);
void readVTK2(int & H1, int & H2, R & PR_SPHR , bool * LIST_N, const char * filename);
void readVTK2map(int HX,int HY,int & H1, int & H2, R & PR_SPHR , bool * LIST_N, const char * filename);
void readVTK3(int & H1, int & H2, int & H3, R & PR_SPHR , bool * LIST_N, const char * filename);
void expVTK2(bool * LIST_N, int H1, int H2, int N_TOT,R PR_DISR);
void expVTK3(bool * LIST_N, int H1, int H2, int H3, int N_TOT,R PR_SPHR);
void ExpCTR3(R * sig, int H1, int H2, int H3, int N_TOT);
void ExpTXT3(R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_R, int NB_INC, int H1, int H2, int H3);

#endif
