#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include "VTK.h"
#include "../PROCESS/ftstring.h"

using namespace std;

void read_H(int & H_TOT,const char * filename){

ifstream f(filename,ios::in);
f >> H_TOT;
f.close();

}    

void readmap3(int H1, int & NB_SPH, R & L_VER, R & R_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, const char * filename){

//cout<<"Lecture du fichier de donnees : "<<filename;
ifstream f(filename,ios::in);

f >> NB_SPH >> L_VER >> R_SPH;

R_SPH=R_SPH*R(H1)/L_VER;

cout<<"NB_SPH : "<<NB_SPH<<", L_VER : "<<L_VER<<", R_SPH : "<<R_SPH<<endl;

for (int i=0;i<NB_SPH;i++){

f >> std::scientific>> (LIST_X[i])  >> std::scientific>> (LIST_Y[i]) >> std::scientific>> (LIST_Z[i]);  

LIST_X[i]=(LIST_X[i])*R(H1)/L_VER;
LIST_Y[i]=(LIST_Y[i])*R(H1)/L_VER;
LIST_Z[i]=(LIST_Z[i])*R(H1)/L_VER;
//cout<<"Fibre "<<i<<" - "<<LIST_X[i]<<", "<<LIST_Y[i]<<", "<<LIST_Z[i]<<", "<<LIST_PH[i]<<", "<<LIST_GA[i]<<", "<<LIST_PS[i]<<endl; 
}

//cout<<"Fin de lecture du fichier de donnees "<<filename<<endl<<endl;
f.close();

}


void readVTK2(int & H1, int & H2, R & PR_SPHR , bool * LIST_N, const char * filename){

PR_SPHR=0.;
  
//cout<<"Lecture du fichier de donnees : "<<filename;
ifstream f(filename,ios::in);

f >> H1;
H2=H1;
for (int i=0;i<(H1*H2);i++){
f >> LIST_N[i];
if(LIST_N[i]){PR_SPHR++;}
}
PR_SPHR/=(H1*H2);

//cout<<"Fin de lecture du fichier de donnees "<<filename<<endl<<endl;
f.close();

cout<<"H1 : "<<H1<<", H2 : "<<H2<<", PR_SPHR :"<<PR_SPHR<<endl;

}  

void readVTK2map(int HX,int HY,int & H1, int & H2, R & PR_SPHR , bool * LIST_N, const char * filename){

PR_SPHR=0.;
  
  
  
//cout<<"Lecture du fichier de donnees : "<<filename;
ifstream f(filename,ios::in);

int i=0;
int tmp;
for (int j=0;j<(HX*HY);j++){
int NX=j%HX;
int NY=j/HX;

f >> tmp;
//cout<<NX<<", "<<NY<<endl;
if((NX<H1)&&(NY<H2)){
LIST_N[i]=bool(tmp/255);
if(LIST_N[i]>0){LIST_N[i]=1;PR_SPHR++;}
i++;
}

}
PR_SPHR/=(H1*H2);
cout<<i<<endl;
//cout<<"Fin de lecture du fichier de donnees "<<filename<<endl<<endl;
f.close();

cout<<"H1 : "<<H1<<", H2 : "<<H2<<", PR_SPHR : "<<PR_SPHR<<endl;

}


void readVTK3(int & H1, int & H2, int & H3, R & PR_SPHR , bool * LIST_N, const char * filename){

PR_SPHR=0.;
  
//cout<<"Lecture du fichier de donnees : "<<filename;
ifstream f(filename,ios::in);

f >> H1;
H2=H1;
H3=H1;
int tmp;
for (int i=0;i<(H1*H2*H3);i++){
f >> tmp;
LIST_N[i]=(tmp>0?0:1);
if(LIST_N[i]){PR_SPHR++;}
}
PR_SPHR/=(H1*H2*H3);

//cout<<"Fin de lecture du fichier de donnees "<<filename<<endl<<endl;
f.close();

cout<<"H1 : "<<H1<<", H2 : "<<H2<<", H3 : "<<H3<<", PR_SPHR :"<<PR_SPHR<<endl;

}  

void expVTK2(bool * LIST_N, int H1, int H2, int N_TOT,R PR_DISR){

R space=1./(R(H1));

fstream g;

string pp=to_string(PR_DISR);
string filename="MAP/maille2_"+pp+".vtk";
const char * filename1=filename.c_str();

g.open(filename1,fstream::out);

g<<"# vtk DataFile Version 3.0"<<endl;
g<<"craft output"<<endl;
g<<"ASCII"<<endl;
g<<"DATASET STRUCTURED_POINTS"<<endl;
g<<"DIMENSIONS "<<to_string(H1)<<" "<<to_string(H2)<<" "<<to_string(1)<<endl; 
g<<"ORIGIN 0.000000 0.000000 0.000000"<<endl;
g<<"SPACING  "<<to_string(space)<<" "<<to_string(space)<<" "<<to_string(space)<<endl;  
g<<"POINT_DATA "<<to_string(N_TOT)<<endl;
g<<"SCALARS scalars float"<<endl;
g<<"LOOKUP_TABLE default"<<endl;
for(int i=0;i<N_TOT;i++){
int mat1;
  
if(LIST_N[i]==0){mat1=0;}else{mat1=1;}
  
g<<mat1<<" ";
if(((i+1)%H1)==0){g<<endl;}
}


g.close();

}

void expVTK3(bool * LIST_N, int H1, int H2, int H3, int N_TOT,R PR_SPHR){

R space=1./(R(H1));

fstream g;

string pp=to_string(PR_SPHR);
string filename="MAP/maille3_"+pp+".vtk";
const char * filename1=filename.c_str();

g.open(filename1,fstream::out);

g<<"# vtk DataFile Version 3.0"<<endl;
g<<"craft output"<<endl;
g<<"ASCII"<<endl;
g<<"DATASET STRUCTURED_POINTS"<<endl;
g<<"DIMENSIONS "<<to_string(H1)<<" "<<to_string(H2)<<" "<<to_string(H3)<<endl; 
g<<"ORIGIN 0.000000 0.000000 0.000000"<<endl;
g<<"SPACING  "<<to_string(space)<<" "<<to_string(space)<<" "<<to_string(space)<<endl;  
g<<"POINT_DATA "<<to_string(N_TOT)<<endl;
g<<"SCALARS scalars float"<<endl;
g<<"LOOKUP_TABLE default"<<endl;
for(int i=0;i<N_TOT;i++){
g<<LIST_N[i]<<" ";
if(((i+1)%H1)==0){g<<endl;}
}

g.close();

}

void ExpCTR3(R * sig, int H1, int H2, int H3, int N_TOT){

string filename="POST/Sig.vtk";
const char * signame =filename.c_str(); 
fstream g;

g.open(signame,fstream::out);
g<<to_string(H1)<<" "<<to_string(H2)<<" "<<to_string(H3)<<endl;
for(int i=0;i<N_TOT;i++){
g<<sig[i]<<" ";
if(((i+1)%H1)==0){g<<endl;}
}


g.close();

} 

inline void ExpTXT3(R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_R, int NB_INC, int H1, int H2, int H3){
	
 //post-traitement
  string filename="MAP/Map3.txt";
  const char * filename1 = filename.c_str();

  fstream h;  
  h.open(filename1,fstream::out);
  h<<NB_INC<<endl;
  h<<LIST_R[0]/H1<<endl;  
  for(int it=0;it<NB_INC;it++){
  h<<LIST_X[it]/H1<<'\t'<<LIST_Y[it]/H2<<'\t'<<LIST_Z[it]/H3<<endl;
  } 
  h<<endl;
  h.close();	
	
}
