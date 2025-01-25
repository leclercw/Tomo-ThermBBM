#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <sys/resource.h> 
#include <string.h>

#ifndef __CONF__
#define __CONF__



void allocat_memory_elas2D(){

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
//	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }

LIST_R = new R[NB_DIS];
LIST_X = new R[NB_DIS];
LIST_Y = new R[NB_DIS];

LIST_N = new bool[N_TOT];
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}
LIST_A = new int[N_TOT];

LIST_SXX = new R[N_TOT];
LIST_SXY = new R[N_TOT];
LIST_SYY = new R[N_TOT];
LIST_UX = new R[N_TOT];
LIST_UY = new R[N_TOT];

LIST_V= new int [(H1+1)*(H2+1)];
LIST_Q = new int *[2];
for(int i=0;i<2;i++){
LIST_Q[i] = new int [N_TOT];
}


}

void allocat_memory_ther2D(){

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
//	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }

LIST_R = new R[NB_DIS];
LIST_X = new R[NB_DIS];
LIST_Y = new R[NB_DIS];

LIST_N = new bool[N_TOT];
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}
LIST_A = new int[N_TOT];

LIST_FX = new R[N_TOT];
LIST_FY = new R[N_TOT];
LIST_T = new R[N_TOT];

LIST_V= new int [(H1+1)*(H2+1)];
LIST_Q = new int *[2];
for(int i=0;i<2;i++){
LIST_Q[i] = new int [N_TOT];
}


}

void allocat_memory_elas3D(){

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
//	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }

LIST_R = new R[NB_SPH];
LIST_X = new R[NB_SPH];
LIST_Y = new R[NB_SPH];
LIST_Z = new R[NB_SPH];

LIST_N = new bool[N_TOT];
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}
LIST_A = new int[N_TOT];

LIST_PH = new R[N_TOT];
LIST_PS = new R[N_TOT];
LIST_GA = new R[N_TOT];

LIST_SXX = new R[N_TOT];
LIST_SXY = new R[N_TOT];
LIST_SXZ = new R[N_TOT];
LIST_SYY = new R[N_TOT];
LIST_SYZ = new R[N_TOT];
LIST_SZZ = new R[N_TOT];

LIST_UX = new R[N_TOT];
LIST_UY = new R[N_TOT];
LIST_UZ = new R[N_TOT];

LIST_V= new int [(H1+1)*(H2+1)*(H3+1)];
LIST_Q = new int *[4];
for(int i=0;i<4;i++){
LIST_Q[i] = new int [N_TOT];
}


}

void allocat_memory_ther3D(){

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
//	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }

LIST_R = new R[NB_SPH];
LIST_X = new R[NB_SPH];
LIST_Y = new R[NB_SPH];
LIST_Z = new R[NB_SPH];

LIST_PH = new R[N_TOT];
LIST_PS = new R[N_TOT];
LIST_GA = new R[N_TOT];

LIST_N = new bool[N_TOT];
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}
LIST_A = new int[N_TOT];

LIST_FX = new R[N_TOT];
LIST_FY = new R[N_TOT];
LIST_FZ = new R[N_TOT];
LIST_T = new R[N_TOT];

LIST_V= new int [(H1+1)*(H2+1)*(H3+1)];
LIST_Q = new int *[4];
for(int i=0;i<4;i++){
LIST_Q[i] = new int [N_TOT];
}


}


////////////////////////////////////////////////////////////////////////////////

void free_memory_elas2D(){

delete [] LIST_R;
delete [] LIST_X;
delete [] LIST_Y;

delete [] LIST_N;
delete [] LIST_A;

delete [] LIST_SXX;
delete [] LIST_SXY;
delete [] LIST_SYY;
delete [] LIST_UX;
delete [] LIST_UY;

delete [] LIST_V;
for(int i=0;i<2;i++){
delete [] LIST_Q[i];
}

}

void free_memory_ther2D(){

delete [] LIST_R;
delete [] LIST_X;
delete [] LIST_Y;

delete [] LIST_N;
delete [] LIST_A;

delete [] LIST_FX;
delete [] LIST_FY;
delete [] LIST_T;

delete [] LIST_V;
for(int i=0;i<2;i++){
delete [] LIST_Q[i];
}

}

void free_memory_elas3D(){

delete [] LIST_R;
delete [] LIST_X;
delete [] LIST_Y;
delete [] LIST_Z;

delete [] LIST_N;
delete [] LIST_A;

delete [] LIST_SXX;
delete [] LIST_SXY;
delete [] LIST_SXZ;
delete [] LIST_SYY;
delete [] LIST_SYZ;
delete [] LIST_SZZ;
delete [] LIST_UX;
delete [] LIST_UY;
delete [] LIST_UZ;

delete [] LIST_V;
for(int i=0;i<4;i++){
delete [] LIST_Q[i];
}

}

void free_memory_ther3D(){

delete [] LIST_R;
delete [] LIST_X;
delete [] LIST_Y;
delete [] LIST_Z;

delete [] LIST_N;
delete [] LIST_A;

delete [] LIST_FX;
delete [] LIST_FY;
delete [] LIST_FZ;
delete [] LIST_T;

delete [] LIST_V;
for(int i=0;i<4;i++){
delete [] LIST_Q[i];
}

}


























#endif

