#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/resource.h> 

#include "../PREPRO/VTK.h"
#include "../PREPRO/genesis.h"
#include "../PREPRO/micro.h"
#include "../PROCESS/reso_FFT.h"
#include "../PROCESS/ftstring.h"
#include "../PROCESS/ctime.h"

#include "view.h"
#include "conf.h"

typedef double R;
const R Pi=3.14159265;

using namespace std;

int main(int argc, char *argv[])
{
R t1=CPUtime();

/*
// RESOLUTION 1
string str_H1 = argv[1];
from_string<int>(H1,str_H1,dec);

// RESOLUTION 2
string str_H2 = argv[2];
from_string<int>(H2,str_H2,dec);

// RESOLUTION 3
string str_H3 = argv[3];
from_string<int>(H3,str_H3,dec);
*/

V_TOT = H1*H2;
N_TOT = H1*H2*H3;
cH=min(H1,H2);
cH=min(cH,H3);


// Nom du fichier microstructure
string filename = argv[1];
const char * filename1 = filename.c_str();

// CONTRASTE
string str_cont = argv[2];
from_string<R>(cont,str_cont,dec);

// NOMBRE D INCLUSIONS
//string str_NB_SPH = argv[5];
//from_string<int>(NB_SPH,str_NB_SPH,dec);
//string str_NB_FIB = argv[5];
//from_string<int>(NB_FIB,str_NB_FIB,dec);

// RAPPORT DES ECHELLES
//string str_R_ECH = argv[6];
//from_string<R>(R_ECH,str_R_ECH,dec);

// FACTEUR DE FORME
//string str_FACF = argv[7];
//from_string<R>(FACF,str_FACF,dec);

read_H(H1,filename1);
H2=H1;
H3=H1;
V_TOT = H1*H2;
N_TOT = H1*H2*H3;
cH=min(H1,H2);
cH=min(cH,H3);

// ALLOCATION DYNAMIQUE
allocat_memory_ther3D();

// LECTURE DE LA GRILLE
readVTK3(H1,H2,H3,PR_INC,LIST_N,filename1);

// GENERATION DE LA GRILLE

//gener3_cyl_alea_curve_int_per(10.,10.,NB_FIB,FACF,R_ECH,N_TOT,V_TOT,H1,H2,H3,LIST_N,LIST_PS,LIST_GA,LIST_PH);
//gener3_cyl_alea(50,NB_FIB,FACF,R_ECH,0.,N_TOT,V_TOT,H1,H2,H3,LIST_N,LIST_PS,LIST_GA,LIST_PH);
//gener3_cyl_align_per(50,NB_FIB,FACF,R_ECH,N_TOT,V_TOT,H1,H2,H3,LIST_N,1);
//gener3_cyl_inf_per(NB_FIB, R_ECH,5., N_TOT,V_TOT,H1,H2,H3,LIST_N,3);
//gener3_sphere_alea(NB_SPH,R_ECH,10.,5.,N_TOT,V_TOT,H1,H2,H3,LIST_N,LIST_R,LIST_X,LIST_Y,LIST_Z);
treat_per3(H1,H2,H3,LIST_N);

// COMPACITE
eval_fv(PR_INC, N_TOT, LIST_N);
pC=PR_INC;
cout<<"Pourcentage d inclusions après traitement: "<<PR_INC<<endl;

// REPERAGE CONTOURS
int N_QUA;
int IN=1;
int OUT=0;

// Search surf
search_surf(N_QUA,LIST_N,H1,H2,H3,IN,OUT);
nQ=N_QUA;

// Treat surf
treat_surf(N_QUA,LIST_N,LIST_Q,LIST_V,H1,H2,H3,IN,OUT);

expVTK3(LIST_N,H1,H2,H3,N_TOT,PR_INC);


//// HOMOGENEISATION

// Caractéristiques matériaux
// Matrice
R K1  = 1.;
cout <<"Conductivité thermique de la matrice: "<<K1<<endl;

// inclusion 
R K2  = cont;
cout <<"Conductivité thermique de(s) inclusions(s): "<<K2<<endl;
      
//milieu 0
R K0 =-sqrt(K1*K2);

R K1pK0=K1+K0;
R K2pK0=K2+K0;
R iK1=2*K0/(K1-K0);
R iK2=2*K0/(K2-K0);

// Résolution fft3
R flu11[3];
R flu22[3];
R flu33[3];

/*
reso_fft3_ther_temp(LIST_N,LIST_T,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu11,1,mintt,maxtt);
reso_fft3_ther_temp(LIST_N,LIST_T,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu22,2,mintt,maxtt);
reso_fft3_ther_temp(LIST_N,LIST_T,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu33,3,mintt,maxtt);
*/
reso_fft3_ther_flux(LIST_N,LIST_FX,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu11,1,1,minfx,maxfx);
reso_fft3_ther_flux(LIST_N,LIST_FY,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu22,2,2,minfy,maxfy);
reso_fft3_ther_flux(LIST_N,LIST_FY,H1,H2,H3,N_TOT,K0,K1,K2,K1pK0,K2pK0,iK1,iK2,flu33,3,3,minfy,maxfy);

cout<<"Phix :"<<minfx<<", "<<maxfx<<endl;
//cout<<"Temp :"<<mintt<<", "<<maxtt<<endl;

cout<<"Conductivités thermiques estimées :"<<endl;
cout<<"K1 : "<<flu11[0]<<endl;
cout<<"K2 : "<<flu22[1]<<endl;
cout<<"K3 : "<<flu33[2]<<endl;

fstream g;

string *str;
str = new string[6];

str[0] = to_string(H1); str[1] = to_string(PR_INC);str[2] = to_string(cont);
str[3] = to_string(flu11[0]*0.0257); str[4] = to_string(flu22[1]*0.0257);str[5] = to_string(flu33[2]*0.0257);

g.open("resultsmat.txt",fstream::out | fstream::app);
 for(int i=0;i<6;i++)
 {
 g<<str[i]<<'\t';
 }
g<<endl;
g.close();


R t2=CPUtime()-t1;
cout<<"Temps de calcul CPU : "<<t2<<endl;


// VISU

//view_init3(argc,argv);
//view_micro3_surf();
//view_micro3_phix(1,cH/2);
//view_micro3_phiy(1,cH/2);
//view_micro3_phiz(1,cH/2);
//view_micro3_t(1,cH/2);

// LIBERATION MEMOIRE
free_memory_ther3D();

return 1;
}
