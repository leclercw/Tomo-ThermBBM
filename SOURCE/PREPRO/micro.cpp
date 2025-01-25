#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "micro.h"

using namespace std;

// Function search contour
void search_cont(int & N_SEG, bool * LIST_N, int H1, int H2, int IN, int OUT)
{
N_SEG=0;

	for (int i=0;i<H1*H2;i++){
	  
			if(LIST_N[i]==IN){
			int l,c,NUMV,NX,NY,NXP,NYP;

			l = i/H1;
			c = i%H1;

			// Neighbour cubes  
			// Pixel 0
			NX=c;
			NY=l-1;
			NYP=(NY+H2)%H2;
			NUMV=NY*H1+NX;

			if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
			N_SEG++;   
			}
			else if(NY!=NYP){
			N_SEG++;    
			}
			// Pixel 1
			NX=c-1;
			NY=l;
			NXP=(NX+H1)%H1;
			NUMV=NY*H1+NX;

			if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
			N_SEG++;   
			}
			else if(NX!=NXP){
			N_SEG++;    
			}			
			// Pixel 2
			NX=c+1;
			NY=l;
			NXP=(NX+H1)%H1;
			NUMV=NY*H1+NX;

			if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
			N_SEG++;   
			}
			else if(NX!=NXP){
			N_SEG++;    
			}		
			// Pixel 3
			NX=c;
			NY=l+1;
			NYP=(NY+H2)%H2;
			NUMV=NY*H1+NX;

			if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
			N_SEG++;   
			}
			else if(NY!=NYP){
			N_SEG++;    
			}	

     	   }

	}

}

// Function treat contour
void treat_cont(int N_SEG, bool * LIST_N, int ** LIST_Q, int * LIST_V, int H1, int H2, int IN, int OUT)
{
	
for(int i=0;i<((H1+1)*(H2+1));i++){
LIST_V[i]=-1;
}

int eq=0;


	for (int i=0;i<H1*H2;i++){
	  
			if(LIST_N[i]==IN){
			int l,c,NUMV,NX,NY,NXP,NYP,n0,n1,n2,n3;

			l = i/H1;
			c = i%H1;

			n0 = (l+0)*(H1+1)+(c+0);
			n1 = (l+0)*(H1+1)+(c+1); 
			n2 = (l+1)*(H1+1)+(c+0);
			n3 = (l+1)*(H1+1)+(c+1); 

			// Neighbour cubes  
			// Pixel 0
			NX=c;
			NY=l-1;
			NYP=(NY+H2)%H2;
			NUMV=NY*H1+NX;

			if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
			LIST_Q[0][eq]=n0;
            LIST_Q[1][eq]=n1;          
			LIST_V[n0]=1;
			LIST_V[n1]=1;              
			eq++;    
			}
			else if(NY!=NYP){
			LIST_Q[0][eq]=n0;
            LIST_Q[1][eq]=n1;          
			LIST_V[n0]=1;
			LIST_V[n1]=1;              
			eq++;     
			}
			// Pixel 1
			NX=c-1;
			NY=l;
			NXP=(NX+H1)%H1;
			NUMV=NY*H1+NX;

			if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
			LIST_Q[0][eq]=n0;
            LIST_Q[1][eq]=n2;          
			LIST_V[n0]=1;
			LIST_V[n2]=1;              
			eq++;   
			}
			else if(NX!=NXP){
			LIST_Q[0][eq]=n0;
            LIST_Q[1][eq]=n2;          
			LIST_V[n0]=1;
			LIST_V[n2]=1;  
			eq++;    
			}			
			// Pixel 2
			NX=c+1;
			NY=l;
			NXP=(NX+H1)%H1;
			NUMV=NY*H1+NX;

			if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
			LIST_Q[0][eq]=n1;
            LIST_Q[1][eq]=n3;          
			LIST_V[n1]=1;
			LIST_V[n3]=1;              
			eq++;   
			}
			else if(NX!=NXP){
			LIST_Q[0][eq]=n1;
            LIST_Q[1][eq]=n3;          
			LIST_V[n1]=1;
			LIST_V[n3]=1;  
			eq++;       
			}		
			// Pixel 3
			NX=c;
			NY=l+1;
			NYP=(NY+H2)%H2;
			NUMV=NY*H1+NX;

			if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
			LIST_Q[0][eq]=n2;
            LIST_Q[1][eq]=n3;          
			LIST_V[n2]=1;
			LIST_V[n3]=1;              
			eq++;    
			}
			else if(NY!=NYP){
			LIST_Q[0][eq]=n2;
            LIST_Q[1][eq]=n3;          
			LIST_V[n2]=1;
			LIST_V[n3]=1;              
			eq++;    
			}	

     	   }

	}

}

// Function search surface

void search_surf(int & N_QUA, bool * LIST_N, int H1,int H2,int H3, int IN, int OUT)
{
N_QUA=0;

for (int i=0;i<H1*H2*H3;i++){
  
if(LIST_N[i]==IN){
int p,l,c,NUMV,NX,NY,NZ,NXP,NYP,NZP;

p = i/(H1*H2);
l = (i%(H1*H2))/H1;
c = i%H1;

// Neighbour cubes  
// Cube 0
NX=c;
NY=l;
NZ=p-1;
NZP=(NZ+H3)%H3;

NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NZ==NZP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NZ!=NZP){
N_QUA++;    
}

// Cube 1
NX=c-1;
NY=l;
NZ=p;
NXP=(NX+H1)%H1;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NX!=NXP){
N_QUA++;    
}

// Cube 2
NX=c;
NY=l-1;
NZ=p;
NYP=(NY+H2)%H2;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NY!=NYP){
N_QUA++;    
}

// Cube 5
NX=c;
NY=l;
NZ=p+1;
NZP=(NZ+H3)%H3;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NZ==NZP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NZ!=NZP){
N_QUA++;    
}

// Cube 3
NX=c+1;
NY=l;
NZ=p;
NXP=(NX+H1)%H1;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NX!=NXP){
N_QUA++;    
}

// Cube 4
NX=c;
NY=l+1;
NZ=p;
NYP=(NY+H2)%H2;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
N_QUA++;   
}
else if(NY!=NYP){
N_QUA++;    
}



}

}

}

// Function treat surface

void treat_surf(int N_QUA, bool * LIST_N, int ** LIST_Q, int * LIST_V, int H1,int H2,int H3, int IN, int OUT)
{

for(int i=0;i<((H1+1)*(H2+1)*(H3+1));i++){
LIST_V[i]=-1;
}


int eq=0;

// Treat surface
for (int i=0;i<H1*H2*H3;i++){
  
if(LIST_N[i]==IN){
int p,l,c,NUMV,NX,NY,NZ,NXP,NYP,NZP,n0,n1,n2,n3,n4,n5,n6,n7;

p = i/(H1*H2);
l = (i%(H1*H2))/H1;
c = i%H1;

n0 = (p+0)*(H1+1)*(H2+1)+(l+0)*(H1+1)+(c+0);
n1 = (p+0)*(H1+1)*(H2+1)+(l+0)*(H1+1)+(c+1); 
n2 = (p+0)*(H1+1)*(H2+1)+(l+1)*(H1+1)+(c+0);
n3 = (p+0)*(H1+1)*(H2+1)+(l+1)*(H1+1)+(c+1); 
n4 = (p+1)*(H1+1)*(H2+1)+(l+0)*(H1+1)+(c+0);
n5 = (p+1)*(H1+1)*(H2+1)+(l+0)*(H1+1)+(c+1); 
n6 = (p+1)*(H1+1)*(H2+1)+(l+1)*(H1+1)+(c+0);
n7 = (p+1)*(H1+1)*(H2+1)+(l+1)*(H1+1)+(c+1); 

// Neighbour cubes  
// Cube 0
NX=c;
NY=l;
NZ=p-1;
NZP=(NZ+H3)%H3;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NZ==NZP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n0;
LIST_Q[1][eq]=n2;
LIST_Q[2][eq]=n3;
LIST_Q[3][eq]=n1;
LIST_V[n0]=1;
LIST_V[n1]=1;
LIST_V[n2]=1;
LIST_V[n3]=1;
eq++;   
}
else if(NZ!=NZP){
LIST_Q[0][eq]=n0;
LIST_Q[1][eq]=n2;
LIST_Q[2][eq]=n3;
LIST_Q[3][eq]=n1;
LIST_V[n0]=1;
LIST_V[n1]=1;
LIST_V[n2]=1;
LIST_V[n3]=1;
eq++;    
}

// Cube 1
NX=c-1;
NY=l;
NZ=p;
NXP=(NX+H1)%H1;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n0;
LIST_Q[1][eq]=n4;
LIST_Q[2][eq]=n6;
LIST_Q[3][eq]=n2;   
LIST_V[n0]=1;
LIST_V[n2]=1;
LIST_V[n4]=1;
LIST_V[n6]=1;
eq++;   
}
else if(NX!=NXP){
LIST_Q[0][eq]=n0;
LIST_Q[1][eq]=n4;
LIST_Q[2][eq]=n6;
LIST_Q[3][eq]=n2; 
LIST_V[n0]=1;
LIST_V[n2]=1;
LIST_V[n4]=1;
LIST_V[n6]=1;
eq++;    
}

// Cube 2
NX=c;
NY=l-1;
NZ=p;
NYP=(NY+H2)%H2;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n4;
LIST_Q[1][eq]=n0;
LIST_Q[2][eq]=n1;
LIST_Q[3][eq]=n5;
LIST_V[n0]=1;
LIST_V[n1]=1;
LIST_V[n4]=1;
LIST_V[n5]=1;
eq++;   
}
else if(NY!=NYP){
LIST_Q[0][eq]=n4;
LIST_Q[1][eq]=n0;
LIST_Q[2][eq]=n1;
LIST_Q[3][eq]=n5;
LIST_V[n0]=1;
LIST_V[n1]=1;
LIST_V[n4]=1;
LIST_V[n5]=1;
eq++;    
}

// Cube 5
NX=c;
NY=l;
NZ=p+1;
NZP=(NZ+H3)%H3;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NZ==NZP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n4;
LIST_Q[1][eq]=n5;
LIST_Q[2][eq]=n7;
LIST_Q[3][eq]=n6;   
LIST_V[n4]=1;
LIST_V[n5]=1;
LIST_V[n6]=1;
LIST_V[n7]=1;
eq++;   
}
else if(NZ!=NZP){
LIST_Q[0][eq]=n4;
LIST_Q[1][eq]=n5;
LIST_Q[2][eq]=n7;
LIST_Q[3][eq]=n6;
LIST_V[n4]=1;
LIST_V[n5]=1;
LIST_V[n6]=1;
LIST_V[n7]=1;
eq++;    
}

// Cube 3
NX=c+1;
NY=l;
NZ=p;
NXP=(NX+H1)%H1;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NX==NXP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n1;
LIST_Q[1][eq]=n3;
LIST_Q[2][eq]=n7;
LIST_Q[3][eq]=n5;
LIST_V[n1]=1;
LIST_V[n3]=1;
LIST_V[n5]=1;
LIST_V[n7]=1;
eq++;   
}
else if(NX!=NXP){
LIST_Q[0][eq]=n1;
LIST_Q[1][eq]=n3;
LIST_Q[2][eq]=n7;
LIST_Q[3][eq]=n5; 
LIST_V[n1]=1;
LIST_V[n3]=1;
LIST_V[n5]=1;
LIST_V[n7]=1;
eq++;    
}

// Cube 4
NX=c;
NY=l+1;
NZ=p;
NYP=(NY+H2)%H2;
NUMV=NZ*(H1*H2)+NY*H1+NX;

if((NY==NYP)&&(LIST_N[NUMV]==OUT)){
LIST_Q[0][eq]=n3;
LIST_Q[1][eq]=n2;
LIST_Q[2][eq]=n6;
LIST_Q[3][eq]=n7;  
LIST_V[n2]=1;
LIST_V[n3]=1;
LIST_V[n6]=1;
LIST_V[n7]=1;
eq++;   
}
else if(NY!=NYP){
LIST_Q[0][eq]=n3;
LIST_Q[1][eq]=n2;
LIST_Q[2][eq]=n6;
LIST_Q[3][eq]=n7;  
LIST_V[n2]=1;
LIST_V[n3]=1;
LIST_V[n6]=1;
LIST_V[n7]=1;
eq++;    
}

}
}

}


void treatbm2(int H1, int H2, bool * LIST_N){

bool * LIST_N2= new bool[H1*H2];

for(long it=0;it<(H1*H2);it++){
LIST_N2[it]=LIST_N[it];
}

for (int i=0;i<(H1*H2);i++){
  if(LIST_N[i]==0){
  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N2=NXP1+H1*NYP1;
  int N3=NX+H1*NYP1;
  int N4=NXM1+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N6=NXM1+H1*NYM1;
  int N7=NX+H1*NYM1;
  int N8=NXP1+H1*NYM1;

  LIST_N2[N1]=0;
  LIST_N2[N2]=0;
  LIST_N2[N3]=0;
  LIST_N2[N4]=0;
  LIST_N2[N5]=0;
  LIST_N2[N6]=0;
  LIST_N2[N7]=0;
  LIST_N2[N8]=0;


  }
}

for(long it=0;it<(H1*H2);it++){
LIST_N[it]=LIST_N2[it];
}

}
 
void treatbp2(int H1, int H2, bool * LIST_N){

bool * LIST_N2= new bool[H1*H2];

for(long it=0;it<(H1*H2);it++){
LIST_N2[it]=LIST_N[it];
}

for (int i=0;i<(H1*H2);i++){
  if(LIST_N[i]==1){
  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N2=NXP1+H1*NYP1;
  int N3=NX+H1*NYP1;
  int N4=NXM1+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N6=NXM1+H1*NYM1;
  int N7=NX+H1*NYM1;
  int N8=NXP1+H1*NYM1;

  LIST_N2[N1]=1;
  LIST_N2[N2]=1;
  LIST_N2[N3]=1;
  LIST_N2[N4]=1;
  LIST_N2[N5]=1;
  LIST_N2[N6]=1;
  LIST_N2[N7]=1;
  LIST_N2[N8]=1;

  }
}

for(long it=0;it<(H1*H2);it++){
LIST_N[it]=LIST_N2[it];
}

}

void treatbm2b(int H1, int H2, bool * LIST_N){

bool * LIST_N2= new bool[H1*H2];

for(long it=0;it<(H1*H2);it++){
LIST_N2[it]=LIST_N[it];
}

for (int i=0;i<(H1*H2);i++){
  if(LIST_N[i]==0){
  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N3=NX+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N7=NX+H1*NYM1;

  LIST_N2[N1]=0;
  LIST_N2[N3]=0;
  LIST_N2[N5]=0;
  LIST_N2[N7]=0;

  }
}

for(long it=0;it<(H1*H2);it++){
LIST_N[it]=LIST_N2[it];
}

}
 
void treatbp2b(int H1, int H2, bool * LIST_N){

bool * LIST_N2= new bool[H1*H2];

for(long it=0;it<(H1*H2);it++){
LIST_N2[it]=LIST_N[it];
}

for (int i=0;i<(H1*H2);i++){
  if(LIST_N[i]==1){
  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N3=NX+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N7=NX+H1*NYM1;

  LIST_N2[N1]=1;
  LIST_N2[N3]=1;
  LIST_N2[N5]=1;
  LIST_N2[N7]=1;

  }
}

for(long it=0;it<(H1*H2);it++){
LIST_N[it]=LIST_N2[it];
}

}

void detect_agr2(int H1, int H2, int * LIST_A, bool * LIST_N, int N_TOT){

int jt=1;
for(int it=0;it<N_TOT;it++){
  LIST_A[it]=0;	
  if(LIST_N[it]>0){
  LIST_A[it]=jt;
  jt++;
  }
}

int NBA_old,NBA;
NBA_old=0;
NBA=jt;

cout<<endl;

bool boola=1;

while(boola){

boola=0;

for (int i=0;i<(H1*H2);i++){
  if(LIST_A[i]>0){

  int valmin=LIST_A[i];

  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N2=NXP1+H1*NYP1;
  int N3=NX+H1*NYP1;
  int N4=NXM1+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N6=NXM1+H1*NYM1;
  int N7=NX+H1*NYM1;
  int N8=NXP1+H1*NYM1;

  if(LIST_A[N1]>0) valmin=((valmin<LIST_A[N1])?valmin:LIST_A[N1]);
  if(LIST_A[N2]>0) valmin=((valmin<LIST_A[N2])?valmin:LIST_A[N2]);
  if(LIST_A[N3]>0) valmin=((valmin<LIST_A[N3])?valmin:LIST_A[N3]);
  if(LIST_A[N4]>0) valmin=((valmin<LIST_A[N4])?valmin:LIST_A[N4]);
  if(LIST_A[N5]>0) valmin=((valmin<LIST_A[N5])?valmin:LIST_A[N5]);
  if(LIST_A[N6]>0) valmin=((valmin<LIST_A[N6])?valmin:LIST_A[N6]);
  if(LIST_A[N7]>0) valmin=((valmin<LIST_A[N7])?valmin:LIST_A[N7]);
  if(LIST_A[N8]>0) valmin=((valmin<LIST_A[N8])?valmin:LIST_A[N8]);

  if(LIST_A[i]!=valmin) {LIST_A[i]=valmin;boola=1;}

  if((LIST_A[N1]>0)&&(LIST_A[N1]!=valmin)) {LIST_A[N1]=valmin;boola=1;}
  if((LIST_A[N2]>0)&&(LIST_A[N2]!=valmin)) {LIST_A[N2]=valmin;boola=1;}
  if((LIST_A[N3]>0)&&(LIST_A[N3]!=valmin)) {LIST_A[N3]=valmin;boola=1;}
  if((LIST_A[N4]>0)&&(LIST_A[N4]!=valmin)) {LIST_A[N4]=valmin;boola=1;}
  if((LIST_A[N5]>0)&&(LIST_A[N5]!=valmin)) {LIST_A[N5]=valmin;boola=1;}
  if((LIST_A[N6]>0)&&(LIST_A[N6]!=valmin)) {LIST_A[N6]=valmin;boola=1;}
  if((LIST_A[N7]>0)&&(LIST_A[N7]!=valmin)) {LIST_A[N7]=valmin;boola=1;}
  if((LIST_A[N8]>0)&&(LIST_A[N8]!=valmin)) {LIST_A[N8]=valmin;boola=1;}

  }
}


for (int i=(H1*H2)-1;i>=0;i--){
  if(LIST_A[i]>0){

  int valmin=LIST_A[i];

  int NX=i%H1;
  int NY=i/H1;
  int NXP1=(NX+H1+1)%H1;
  int NXM1=(NX+H1-1)%H1; 
  int NYP1=(NY+H2+1)%H2;
  int NYM1=(NY+H2-1)%H2; 

  int N1=NXP1+H1*NY;
  int N2=NXP1+H1*NYP1;
  int N3=NX+H1*NYP1;
  int N4=NXM1+H1*NYP1;
  int N5=NXM1+H1*NY;
  int N6=NXM1+H1*NYM1;
  int N7=NX+H1*NYM1;
  int N8=NXP1+H1*NYM1;

  if(LIST_A[N1]>0) valmin=((valmin<LIST_A[N1])?valmin:LIST_A[N1]);
  if(LIST_A[N2]>0) valmin=((valmin<LIST_A[N2])?valmin:LIST_A[N2]);
  if(LIST_A[N3]>0) valmin=((valmin<LIST_A[N3])?valmin:LIST_A[N3]);
  if(LIST_A[N4]>0) valmin=((valmin<LIST_A[N4])?valmin:LIST_A[N4]);
  if(LIST_A[N5]>0) valmin=((valmin<LIST_A[N5])?valmin:LIST_A[N5]);
  if(LIST_A[N6]>0) valmin=((valmin<LIST_A[N6])?valmin:LIST_A[N6]);
  if(LIST_A[N7]>0) valmin=((valmin<LIST_A[N7])?valmin:LIST_A[N7]);
  if(LIST_A[N8]>0) valmin=((valmin<LIST_A[N8])?valmin:LIST_A[N8]);

  if(LIST_A[i]!=valmin) {LIST_A[i]=valmin;boola=1;}

  if((LIST_A[N1]>0)&&(LIST_A[N1]!=valmin)) {LIST_A[N1]=valmin;boola=1;}
  if((LIST_A[N2]>0)&&(LIST_A[N2]!=valmin)) {LIST_A[N2]=valmin;boola=1;}
  if((LIST_A[N3]>0)&&(LIST_A[N3]!=valmin)) {LIST_A[N3]=valmin;boola=1;}
  if((LIST_A[N4]>0)&&(LIST_A[N4]!=valmin)) {LIST_A[N4]=valmin;boola=1;}
  if((LIST_A[N5]>0)&&(LIST_A[N5]!=valmin)) {LIST_A[N5]=valmin;boola=1;}
  if((LIST_A[N6]>0)&&(LIST_A[N6]!=valmin)) {LIST_A[N6]=valmin;boola=1;}
  if((LIST_A[N7]>0)&&(LIST_A[N7]!=valmin)) {LIST_A[N7]=valmin;boola=1;}
  if((LIST_A[N8]>0)&&(LIST_A[N8]!=valmin)) {LIST_A[N8]=valmin;boola=1;}



  }
}

int valmax=0;
for (int i=0;i<(H1*H2);i++){
if(LIST_A[i]>0) valmax=(valmax>LIST_A[i])?valmax:LIST_A[i];
}

int LIST_TMP[valmax];
for (int i=0;i<valmax;i++){
LIST_TMP[i]=0;
}

for (int i=0;i<(H1*H2);i++){
if(LIST_A[i]>0) LIST_TMP[LIST_A[i]-1]=1;
}

cout<<"Nombre d'agregats potentiel:"<<valmax<<endl;
int j=0;
for (int i=0;i<valmax;i++){
	if(LIST_TMP[i]>0) {
	j++;
	LIST_TMP[i]=j; 
	}
}
cout<<"Nombre d'agregats potentiel:"<<j<<endl<<endl;

for (int i=0;i<(H1*H2);i++){
if(LIST_A[i]>0) LIST_A[i]=LIST_TMP[LIST_A[i]-1];
}

NBA_old=NBA;
NBA=j;
}


cout<<"Nombre final d'agregats:"<<NBA<<endl<<endl;

}


void treat_per2(int H1, int H2, bool * LIST_N)
{  
 // Liste des noeuds à conditions aux limites
int i1=H2-2;
int i2=H1-2;
int i3=H2-2;
int i4=H1-2;
int i5=1;
int i6=1;
int i7=1;
int i8=1;

int P1[i1];
int P2[i2];
int P3[i3];
int P4[i4];
int P5[i5];
int P6[i6];
int P7[i7];
int P8[i8];

for (int it=1;it<(H1-1);it++){
P2[it-1]=it;
P4[it-1]=H1*(H2-1)+it;
}	

for (int it=1;it<(H2-1);it++){	
P1[it-1]=it*H1;
P3[it-1]=(it+1)*H1-1;
}


P5[0]=0;
P6[0]=H1-1;
P7[0]=H1*H2-1;
P8[0]=H1*(H2-1);

// cas bord 1-3
for(int it=0;it<i1;it++){
int a1=LIST_N[P1[it]];
int a3=LIST_N[P3[it]];
LIST_N[P1[it]]=max(a1,a3);
LIST_N[P3[it]]=LIST_N[P1[it]];
}

// cas bord 2-4
for(int it=0;it<i2;it++){
int a2=LIST_N[P2[it]];
int a4=LIST_N[P4[it]];
LIST_N[P2[it]]=max(a2,a4);
LIST_N[P4[it]]=LIST_N[P2[it]];
}

// cas extr 5 6 7 8
for(int it=0;it<i5;it++){
int a5=LIST_N[P5[it]];
int a6=LIST_N[P6[it]];
int a7=LIST_N[P7[it]];
int a8=LIST_N[P8[it]];

LIST_N[P5[it]]=max(max(max(a5,a6),a7),a8);
LIST_N[P6[it]]=LIST_N[P5[it]];
LIST_N[P7[it]]=LIST_N[P5[it]];
LIST_N[P8[it]]=LIST_N[P5[it]];
}

}


void treat_per3(int H1, int H2, int H3, bool * LIST_N)
{  
 // Liste des noeuds à conditions aux limites

//faces 
// 9 => (H1-2)(H2-2)
// 26 => (H1-2)(H2-2)
// 10 => (H2-2)(H3-2)
// 11 => (H1-2)(H3-2)
// 12 => (H2-2)(H3-2)
// 13 => (H1-2)(H3-2)
int i9 = (H1-2)*(H2-2);
int i26= (H1-2)*(H2-2);
int i10= (H2-2)*(H3-2);
int i11= (H1-2)*(H3-2);
int i12= (H2-2)*(H3-2);
int i13= (H1-2)*(H3-2);

//bords
// 1 => (H2-2)
// 2 => (H1-2)
// 3 => (H2-2)
// 4 => (H1-2)
// 14 => (H3-2)
// 15 => (H3-2)
// 16 => (H3-2)
// 17 => (H3-2)
// 18 => (H2-2)
// 19 => (H1-2)
// 20 => (H2-2)
// 21 => (H1-2)
int i1 = (H2-2);
int i2 = (H1-2);
int i3 = (H2-2);
int i4 = (H1-2);
int i14= (H3-2);
int i15= (H3-2);
int i16= (H3-2);
int i17= (H3-2);
int i18= (H2-2);
int i19= (H1-2);
int i20= (H2-2);
int i21= (H1-2);

//extremites
// 5 => 1
// 6 => 1
// 7 => 1
// 8 => 1
// 22 => 1
// 23 => 1
// 24 => 1
// 25 => 1
int i5=1;
int i6=1;
int i7=1;
int i8=1;
int i22=1;
int i23=1;
int i24=1;
int i25=1;

int * P1 = new int[i1];
int * P2 = new int[i2];
int * P3 = new int[i3];
int * P4 = new int[i4];
int * P5 = new int[i5];
int * P6 = new int[i6];
int * P7 = new int[i7];
int * P8 = new int[i8];
int * P9 = new int[i9];
int * P10 = new int[i10];
int * P11 = new int[i11];
int * P12 = new int[i12];
int * P13 = new int[i13];
int * P14 = new int[i14];
int * P15 = new int[i15];
int * P16 = new int[i16];
int * P17 = new int[i17];
int * P18 = new int[i18];
int * P19 = new int[i19];
int * P20 = new int[i20];
int * P21 = new int[i21];
int * P22 = new int[i22];
int * P23 = new int[i23];
int * P24 = new int[i24];
int * P25 = new int[i25];
int * P26 = new int[i26];

for (int it=1;it<(H1-1);it++){
P2[it-1]=it;
P4[it-1]=H1*(H2-1)+it;
P19[it-1]=H1*H2*(H3-1)+it;
P21[it-1]=H1*H2*(H3-1)+H1*(H2-1)+it;	

}

for (int it=1;it<(H2-1);it++){
P1[it-1]=it*H1;
P3[it-1]=(it+1)*H1-1;	
P18[it-1]=H1*H2*(H3-1)+it*H1;
P20[it-1]=H1*H2*(H3-1)+(it+1)*H1-1;
	for (int jt=1;jt<(H1-1);jt++){
	P9[(it-1)*(H1-2)+(jt-1)]=jt+H1*it;
	P26[(it-1)*(H1-2)+(jt-1)]=jt+H1*it+H1*H2*(H3-1);
	}
}

for (int it=1;it<(H3-1);it++){
P14[it-1]=it*H1*H2;
P15[it-1]=it*H1*H2+H1-1;  
P16[it-1]=it*H1*H2+H1*H2-1;  
P17[it-1]=it*H1*H2+H1*(H2-1);	
	for (int jt=1;jt<(H2-1);jt++){
	P10[(it-1)*(H2-2)+(jt-1)]=jt*H1+it*H1*H2;
    P12[(it-1)*(H2-2)+(jt-1)]=jt*H1+H1-1+it*H1*H2;	
    }
	for (int jt=1;jt<(H1-1);jt++){
	P11[(it-1)*(H1-2)+(jt-1)]=jt+it*H1*H2;
	P13[(it-1)*(H1-2)+(jt-1)]=jt+it*H1*H2+H1*(H2-1);  
    }    
}

P5[0]=0;
P6[0]=H1-1;
P7[0]=H1*H2-1;
P8[0]=H1*(H2-1);
P22[0]=H1*H2*(H3-1);
P23[0]=H1*H2*(H3-1)+H1-1;
P24[0]=H1*H2*(H3-1)+H1*H2-1;
P25[0]=H1*H2*(H3-1)+H1*(H2-1);

// cas face 9-26
for(int it=0;it<i9;it++){
int a1=LIST_N[P9[it]];
int a2=LIST_N[P26[it]];
LIST_N[P9[it]]=max(a1,a2);
LIST_N[P26[it]]=LIST_N[P9[it]];
}

// cas face 10-12
for(int it=0;it<i10;it++){
int a1=LIST_N[P10[it]];
int a2=LIST_N[P12[it]];
LIST_N[P10[it]]=max(a1,a2);
LIST_N[P12[it]]=LIST_N[P10[it]];
}

// cas face 11-13
for(int it=0;it<i11;it++){
int a1=LIST_N[P11[it]];
int a2=LIST_N[P13[it]];
LIST_N[P11[it]]=max(a1,a2);
LIST_N[P13[it]]=LIST_N[P11[it]];
}

// cas bord 1 3 18 20
for(int it=0;it<i1;it++){
int a1=LIST_N[P1[it]];
int a2=LIST_N[P3[it]];
int a3=LIST_N[P18[it]];
int a4=LIST_N[P20[it]];

LIST_N[P1[it]]=max(max(max(a1,a2),a3),a4);
LIST_N[P3[it]]=LIST_N[P1[it]];
LIST_N[P18[it]]=LIST_N[P1[it]];
LIST_N[P20[it]]=LIST_N[P1[it]];
}

// cas bord 2 4 19 21
for(int it=0;it<i2;it++){
int a1=LIST_N[P2[it]];
int a2=LIST_N[P4[it]];
int a3=LIST_N[P19[it]];
int a4=LIST_N[P21[it]];

LIST_N[P2[it]]=max(max(max(a1,a2),a3),a4);
LIST_N[P4[it]]=LIST_N[P2[it]];
LIST_N[P19[it]]=LIST_N[P2[it]];
LIST_N[P21[it]]=LIST_N[P2[it]];
}

// cas bord 14 15 16 17
for(int it=0;it<i14;it++){
int a1=LIST_N[P14[it]];
int a2=LIST_N[P15[it]];
int a3=LIST_N[P16[it]];
int a4=LIST_N[P17[it]];

LIST_N[P14[it]]=max(max(max(a1,a2),a3),a4);
LIST_N[P15[it]]=LIST_N[P14[it]];
LIST_N[P16[it]]=LIST_N[P14[it]];
LIST_N[P17[it]]=LIST_N[P14[it]];
}

// cas bord 5 6 7 8 22 23 24 25
for(int it=0;it<i5;it++){
int a1=LIST_N[P5[it]];
int a2=LIST_N[P6[it]];
int a3=LIST_N[P7[it]];
int a4=LIST_N[P8[it]];
int a5=LIST_N[P22[it]];
int a6=LIST_N[P23[it]];
int a7=LIST_N[P24[it]];
int a8=LIST_N[P25[it]];

LIST_N[P5[it]]=max(max(max(max(max(max(max(a1,a2),a3),a4),a5),a6),a7),a8);
LIST_N[P6[it]]=LIST_N[P5[it]];
LIST_N[P7[it]]=LIST_N[P5[it]];
LIST_N[P8[it]]=LIST_N[P5[it]];
LIST_N[P22[it]]=LIST_N[P5[it]];
LIST_N[P23[it]]=LIST_N[P5[it]];
LIST_N[P24[it]]=LIST_N[P5[it]];
LIST_N[P25[it]]=LIST_N[P5[it]];

}

delete [] P1;delete [] P2;delete [] P3;delete [] P4;
delete [] P5;delete [] P6;delete [] P7;delete [] P8;
delete [] P9;delete [] P10;delete [] P11;delete [] P12;
delete [] P13;delete [] P14;delete [] P15;delete [] P16;
delete [] P17;delete [] P18;delete [] P19;delete [] P20;
delete [] P21;delete [] P22;delete [] P23;delete [] P24;
delete [] P25;delete [] P26;

}  

void eval_fv(R & PR_INC, int N_TOT, bool * LIST_N){
	
PR_INC=0.;
for(int it=0;it<N_TOT;it++){
if(LIST_N[it]) PR_INC++;
}
PR_INC/=N_TOT;
	
}
