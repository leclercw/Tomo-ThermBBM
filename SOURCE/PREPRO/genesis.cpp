#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <list>

#include "genesis.h"

const R Pi=3.14159265;

using namespace std;

void gener2_disk_alea(int NB_DIS,R R_DIS,int N_TOT,int H1, int H2, bool * LIST_N, R * LIST_X, R * LIST_Y, R * LIST_R)
{  

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
      /*

R XX[NB_DIS];
R YY[NB_DIS];

XX[0]=6.67;YY[0]=5.51;
XX[1]=6.73;YY[1]=17.43;
XX[2]=0.86;YY[2]=6.71;
XX[3]=24.06;YY[3]=27.05;
XX[4]=16.92;YY[4]=16.73;
XX[5]=24.61;YY[5]=20.69;*/

 for(int it=0;it<NB_DIS;it++){
R X_GRA = ((R) (rand()%(H1*100)))/R(100);
R Y_GRA = ((R) (rand()%(H2*100)))/R(100);
LIST_X[it]=X_GRA;
LIST_Y[it]=Y_GRA;

LIST_R[it]=R_DIS;
cout<<"it:"<<it<<" "<<LIST_X[it]<<" "<<LIST_Y[it]<<endl;

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);

R DISTX = (X_GRA-(NX+0.5))*(X_GRA-(NX+0.5))+(Y_GRA-(NY+0.5))*(Y_GRA-(NY+0.5));

if(DISTX<=(R_DIS*R_DIS)){
  int NUM = NY*H1+NX;
  LIST_N[NUM]=1; 

int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    
    R BXPLUS = 0;
    R BXMOIN = 0;
    R BYPLUS = 0;
    R BYMOIN = 0;
    R BXI    = 0;
    R BYI    = 0;
    
    if(NXPLUS>=H1){
    BXPLUS = -1; 
    }
    if(NXMOIN<0){
    BXMOIN = 1; 
    }
    if(NYPLUS>=H2){
    BYPLUS = -1; 
    }
    if(NYMOIN<0){
    BYMOIN = 1; 
    }  
    
    NXPLUS = (NXPLUS+H1)%H1;
    NXMOIN = (NXMOIN+H1)%H1;
    NYPLUS = (NYPLUS+H2)%H2;
    NYMOIN = (NYMOIN+H2)%H2;    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI = NY-BVOIS+j;        
    BYI=0;
    
    if(NYI>=H2){
    BYI = -1; 
    }
    if(NYI<0){
    BYI = 1; 
    }
    
    NYI = (NYI+H2)%H2; 
    R DISTXP = (X_GRA+(BXPLUS*H1)-(NXPLUS+0.5))*(X_GRA+(BXPLUS*H1)-(NXPLUS+0.5))+(Y_GRA+(BYI*H2)-(NYI+0.5))*(Y_GRA+(BYI*H2)-(NYI+0.5));
    R DISTXM = (X_GRA+(BXMOIN*H1)-(NXMOIN+0.5))*(X_GRA+(BXMOIN*H1)-(NXMOIN+0.5))+(Y_GRA+(BYI*H2)-(NYI+0.5))*(Y_GRA+(BYI*H2)-(NYI+0.5));

    if(DISTXP<=(R_DIS*R_DIS)){
    int NELXP = NYI*H1+NXPLUS;
    LIST_N[NELXP]=1;
    }

    if(DISTXM<=(R_DIS*R_DIS)){
    int NELXM = NYI*H1+NXMOIN;
    LIST_N[NELXM]=1;
    }

    if((j>0)&&(j<(2*BVOIS))){
    int NXI=NX-BVOIS+j;
    BXI=0;
    if(NXI>=H1){
    BXI = -1; 
    }
    if(NXI<0){
    BXI = 1; 
    }
     
    NXI = (NXI+H1)%H1; 
    R DISTYP = (X_GRA+(BXI*H1)-(NXI+0.5))*(X_GRA+(BXI*H1)-(NXI+0.5))+(Y_GRA+(BYPLUS*H2)-(NYPLUS+0.5))*(Y_GRA+(BYPLUS*H2)-(NYPLUS+0.5));
    R DISTYM = (X_GRA+(BXI*H1)-(NXI+0.5))*(X_GRA+(BXI*H1)-(NXI+0.5))+(Y_GRA+(BYMOIN*H2)-(NYMOIN+0.5))*(Y_GRA+(BYMOIN*H2)-(NYMOIN+0.5));

    if(DISTYP<=(R_DIS*R_DIS)){
      int NELYP = NYPLUS*H1+NXI;
      LIST_N[NELYP]=1;

    }

    if(DISTYM<=(R_DIS*R_DIS)){
      int NELYM = NYMOIN*H1+NXI;
      LIST_N[NELYM]=1;

    } 
                                 
    }       
    
    }
        
    }while(BVOIS<=(R_DIS+2));
                         
}

}

   
}   




bool bincl(R rC,R AA,R BB,R CC, R X1, R Y1, R Z1, R X2, R Y2, R Z2, R XX, R YY, R ZZ){

bool boolinc=0;

R EQUAC1 = (XX-X1)*(XX-X1)+(YY-Y1)*(YY-Y1)+(ZZ-Z1)*(ZZ-Z1);
R EQUAC2 = (AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1))*(AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1));
R EQUAC  = EQUAC1-EQUAC2;

R EQUAP1 = -(AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1));
R EQUAP2 = (AA*(XX-X2)+BB*(YY-Y2)+CC*(ZZ-Z2));

if((EQUAC<=(rC*rC)+1e-8)&&(EQUAP1<=1e-8)&&(EQUAP2<=1e-8)){
boolinc=1;  
}
 
return boolinc;
}



void gener3_cyl_alea_curve_int_per(int THETAM,int PHIM,int nC, R facf, R eche, int N_TOT,int V_TOT,int H1,int H2,int H3, bool * LIST_N, R* LIST_PS, R* LIST_GA, R* LIST_PH)
{  
 
R lC=R(H1)/eche; 
R rC=lC/(2.*facf);
R dC=2.*rC;

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}   
  
R Pi=3.14159265;  
  
// Densité de points définissant la courbure du cylindre
int DENSI2 = floor(lC/(2*dC));

// Frontières des fibres
R ** TBX= new R * [nC];
R ** TBY= new R * [nC];
R ** TBZ= new R * [nC];
R ** TB_ANGT = new R * [nC];
R ** TB_ANGP = new R * [nC];

for(int it=0;it<nC;it++){
  
    // PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
    
///////////////////////////////////////
//A+---------------+E--------------+B//
///////////////////////////////////////

//////////////////////////////////
//                      B       //
//                     /        //
//	             /          //         
//	           /            //
//	         /              //
//	       /                //
//	    +E                  //
//	   /                    //
//      /                       //
//   / ) theta                  //
//A - - ) - - - -               //
//////////////////////////////////

// Etape 1.1 : point E de départ

R X_E = ( rand()/(double)RAND_MAX ) * R(H1);
R Y_E = ( rand()/(double)RAND_MAX ) * R(H2);
R Z_E = ( rand()/(double)RAND_MAX ) * R(H3);

TBX[it]= new R [2*DENSI2+1];
TBY[it]= new R [2*DENSI2+1];
TBZ[it]= new R [2*DENSI2+1];

TB_ANGT[it] = new R[2*DENSI2];
TB_ANGP[it] = new R[2*DENSI2];

int ORIT   = rand()%360;
int ORIP   = rand()%180;

TBX[it][DENSI2] = X_E;
TBY[it][DENSI2] = Y_E;
TBZ[it][DENSI2] = Z_E;
TB_ANGT[it][DENSI2] = ORIT;
TB_ANGP[it][DENSI2] = ORIP;


R X_FIN,Y_FIN,Z_FIN,X_REF2,Y_REF2,Z_REF2;
int AT_REF,AP_REF,AT_REF2,AP_REF2;  

// Etape 1.2 : EB

AT_REF     = ORIT;
AP_REF     = ORIP;
X_FIN     = X_E+(dC*cos(Pi/180*ORIT)*sin(Pi/180*ORIP));
Y_FIN     = Y_E+(dC*sin(Pi/180*ORIT)*sin(Pi/180*ORIP));
Z_FIN     = Z_E+(dC*cos(Pi/180*ORIP));

// cout<<"1:"<<(dC*cos(Pi/180*ORIT)*sin(Pi/180*ORIP))<<endl;
// cout<<"2:"<<(dC*sin(Pi/180*ORIT)*sin(Pi/180*ORIP))<<endl;
// cout<<"3:"<<(dC*cos(Pi/180*ORIP))<<endl;

//X_REF1    = X_E;
X_REF2    = X_FIN;
//Y_REF1    = Y_E;
Y_REF2    = Y_FIN;
//Z_REF1    = Z_E;
Z_REF2    = Z_FIN;
TBX[it][DENSI2+1]    = X_FIN;
TBY[it][DENSI2+1]    = Y_FIN;
TBZ[it][DENSI2+1]    = Z_FIN;
int nump=DENSI2+1;


for(int j=1;j<DENSI2;j++){
 nump++;
AT_REF2 = AT_REF;
AP_REF2 = AP_REF;
AT_REF  = (THETAM>0)?((rand()%(2*THETAM))-THETAM):0;  
AP_REF  = (PHIM>0)?((rand()%(2*PHIM))-PHIM):0;  
AT_REF  = AT_REF+AT_REF2;
AP_REF  = AP_REF+AP_REF2;
AT_REF  = (AT_REF>360)?(AT_REF-360):AT_REF;
AT_REF  = (AT_REF<0)?(AT_REF+360):AT_REF;
AP_REF  = (AP_REF>180)?(360-AP_REF):AP_REF;
AP_REF  = (AP_REF<0)?(-AP_REF):AP_REF;
X_FIN     = X_REF2+(dC*cos(Pi/180*AT_REF)*sin(Pi/180*AP_REF));
Y_FIN     = Y_REF2+(dC*sin(Pi/180*AT_REF)*sin(Pi/180*AP_REF));
Z_FIN     = Z_REF2+(dC*cos(Pi/180*AP_REF));
TB_ANGT[it][nump-1] = AT_REF;
TB_ANGP[it][nump-1] = AP_REF;
TBX[it][nump]    = X_FIN;
TBY[it][nump]    = Y_FIN;
TBZ[it][nump]    = Z_FIN;
//X_REF1    = X_REF2;
X_REF2    = X_FIN;
//Y_REF1    = Y_REF2;
Y_REF2    = Y_FIN;
//Z_REF1    = Z_REF2;
Z_REF2    = Z_FIN;
}

// Etape 1.5 : E
nump=DENSI2;
nump--;
AT_REF     = ORIT;
AP_REF     = ORIP;
TB_ANGT[it][nump] = ORIT;
TB_ANGP[it][nump] = ORIP;
X_FIN     = X_E-(dC*cos(Pi/180*ORIT)*sin(Pi/180*ORIP));
Y_FIN     = Y_E-(dC*sin(Pi/180*ORIT)*sin(Pi/180*ORIP));
Z_FIN     = Z_E-(dC*cos(Pi/180*ORIP));
// cout<<"1:"<<(dC*cos(Pi/180*ORIT)*sin(Pi/180*ORIP))<<endl;
// cout<<"2:"<<(dC*sin(Pi/180*ORIT)*sin(Pi/180*ORIP))<<endl;
// cout<<"3:"<<(dC*cos(Pi/180*ORIP))<<endl;
TBX[it][nump]    = X_FIN;
TBY[it][nump]    = Y_FIN;
TBZ[it][nump]    = Z_FIN;
//X_REF1    = X_E;
X_REF2    = X_FIN;
//Y_REF1    = Y_E;
Y_REF2    = Y_FIN;
//Z_REF1    = Z_E;
Z_REF2    = Z_FIN;

// cout<<"nump:"<<nump<<","<<TBX[it][nump]<<","<<TBY[it][nump]<<","<<TBZ[it][nump]<<endl; 


// Etape 1.6 : EA
for(int j=1;j<DENSI2;j++){
nump--;
AT_REF2 = AT_REF;
AP_REF2 = AP_REF;
AT_REF  = (THETAM>0)?((rand()%(2*THETAM))-THETAM):0;  
AP_REF  = (PHIM>0)?((rand()%(2*PHIM))-PHIM):0;  
AT_REF  = AT_REF+AT_REF2;
AP_REF  = AP_REF+AP_REF2;
AT_REF  = (AT_REF>360)?(AT_REF-360):AT_REF;
AT_REF  = (AT_REF<0)?(AT_REF+360):AT_REF;
AP_REF  = (AP_REF>180)?(360-AP_REF):AP_REF;
AP_REF  = (AP_REF<0)?(-AP_REF):AP_REF;
X_FIN   = X_REF2-(dC*cos(Pi/180*AT_REF)*sin(Pi/180*AP_REF));
Y_FIN   = Y_REF2-(dC*sin(Pi/180*AT_REF)*sin(Pi/180*AP_REF));
Z_FIN   = Z_REF2-(dC*cos(Pi/180*AP_REF));
TB_ANGT[it][nump] = AT_REF;
TB_ANGP[it][nump] = AP_REF;
TBX[it][nump] = X_FIN;
TBY[it][nump] = Y_FIN;
TBZ[it][nump] = Z_FIN;
//X_REF1  = X_REF2;
X_REF2  = X_FIN;
//Y_REF1  = Y_REF2;
Y_REF2  = Y_FIN;
//Z_REF1  = Z_REF2;
Z_REF2  = Z_FIN;
}

for(int jt=0;jt<2*DENSI2;jt++){

//cout<<"angles :"<<THET1<<","<<PHI1<<","<<(TB_ANGT[it][jt]*Pi/180)<<","<<(TB_ANGP[it][jt]*Pi/180)<<endl;
R PH_REF  = TB_ANGT[it][jt]*Pi/180;
R PS_REF  = 1.57079-TB_ANGP[it][jt]*Pi/180;
R GA_REF  = 1.57079;

R AA = cos(PH_REF)*cos(PS_REF)-cos(GA_REF)*sin(PS_REF)*sin(PH_REF);
R BB = sin(PH_REF)*cos(PS_REF)+cos(GA_REF)*sin(PS_REF)*cos(PH_REF);
R CC = sin(GA_REF)*sin(PS_REF);

// cout<<TBX[it][jt]<<","<<TBX[it][jt+1]<<","<<TBX[it][jt]+dC*AA<<","<<TBX[it][jt]+dC*AAA<<endl;
// cout<<TBY[it][jt]<<","<<TBY[it][jt+1]<<","<<TBY[it][jt]+dC*BB<<","<<TBY[it][jt]+dC*BBB<<endl;
// cout<<TBZ[it][jt]<<","<<TBZ[it][jt+1]<<","<<TBZ[it][jt]+dC*CC<<","<<TBZ[it][jt]+dC*CCC<<endl;

TBX[it][jt+1]=TBX[it][jt]+dC*AA;
TBY[it][jt+1]=TBY[it][jt]+dC*BB;
TBZ[it][jt+1]=TBZ[it][jt]+dC*CC;

R XX = (TBX[it][jt]+TBX[it][jt+1])/2;
R YY = (TBY[it][jt]+TBY[it][jt+1])/2;
R ZZ = (TBZ[it][jt]+TBZ[it][jt+1])/2;

int NX    = (int) floor(XX);
int NY    = (int) floor(YY);
int NZ    = (int) floor(ZZ);
int NXP   = (NX+H1)%H1;
int NYP   = (NY+H2)%H2;
int NZP   = (NZ+H3)%H3;

//cout<<jt<<","<<TBX[it][jt]<<","<<TBY[it][jt]<<","<<TBZ[it][jt]<<endl;

bool boolxa = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],XX,YY,ZZ);
int NUMA  = NZP*V_TOT+NYP*H1+NXP;
//cout<<boolxa<<endl;

if((boolxa==1)&&(LIST_N[NUMA]==0)){
LIST_N[NUMA]=1;
LIST_PH[NUMA]=PH_REF;
LIST_GA[NUMA]=GA_REF;
LIST_PS[NUMA]=PS_REF;
}
    //couches sup
    
for(int BVOIS=1;BVOIS<=(floor(1.5*rC));BVOIS++)
    {      
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
     
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;      
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI  = NY-BVOIS+j;       
    int NYIP = (NYI+H2)%H2; 
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI  = NZ-BVOIS+k;        
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXPLUS+0.5),(NYI+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXMOIN+0.5),(NYI+0.5),(NZI+0.5));

	int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
        int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==0)){
        LIST_N[NELXP]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==0)){
        LIST_N[NELXM]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI  = NX-BVOIS+l+1;        
        int NXIP = (NXI+H1)%H1;  
	
	
 	bool boolxp = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXI+0.5),(NYI+0.5),(NZPLUS+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXI+0.5),(NYI+0.5),(NZMOIN+0.5));

        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==0)){
        LIST_N[NELXP]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==0)){
        LIST_N[NELXM]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }        
        
        }  //END  
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI  = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1; 
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;      
        int NZIP = (NZI+H3)%H3;

 	bool boolxp = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXI+0.5),(NYPLUS+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,TBX[it][jt],TBY[it][jt],TBZ[it][jt],TBX[it][jt+1],TBY[it][jt+1],TBZ[it][jt+1],(NXI+0.5),(NYMOIN+0.5),(NZI+0.5));

        int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;	        
	int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
		
        if((boolxp==1)&&(LIST_N[NELXP]==0)){
        LIST_N[NELXP]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==0)){
        LIST_N[NELXM]=1;
	LIST_PH[NUMA]=PH_REF;
	LIST_GA[NUMA]=GA_REF;
	LIST_PS[NUMA]=PS_REF;
        }     
           
        }

      } //END
       
    }

}
  
}


for(int it=0;it<nC;it++){
delete [] TBX[it];
delete [] TBY[it];
delete [] TBZ[it];
delete [] TB_ANGT[it];
delete [] TB_ANGP[it];
}

delete [] TBX;
delete [] TBY;
delete [] TBZ;
delete [] TB_ANGT;
delete [] TB_ANGP;

}


void gener3_sphere_alea_per(int nS, R eche, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N, R* LIST_R, R* LIST_X, R* LIST_Y, R* LIST_Z)
{  

R rS=H1/(2.*eche);

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
       
int it=0;
while(it<nS){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * R(H1);
R Y_GRA = ( rand()/(double)RAND_MAX ) * R(H2);
R Z_GRA = ( rand()/(double)RAND_MAX ) * R(H3);

LIST_X[it]=X_GRA;
LIST_Y[it]=Y_GRA;
LIST_Z[it]=Z_GRA;
LIST_R[it]=rS;

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-H3-LIST_Z[jt])*(LIST_Z[it]-H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;




	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;




	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-H2-LIST_Y[jt])*(LIST_Y[it]-H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-H1-LIST_X[jt])*(LIST_X[it]-H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	distij=(LIST_X[it]+H1-LIST_X[jt])*(LIST_X[it]+H1-LIST_X[jt])+(LIST_Y[it]+H2-LIST_Y[jt])*(LIST_Y[it]+H2-LIST_Y[jt])+(LIST_Z[it]+H3-LIST_Z[jt])*(LIST_Z[it]+H3-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;



	}


if(boolc==0){
// Test intersection voxel

int BVOIS=0;

    do
    {  
    BVOIS++;        
    
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI = NY-BVOIS+j;  
    int NYIP = (NYI+H2)%H2;  
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI = NZ-BVOIS+k;  
	int NZIP = (NZI+H3)%H3;  

        R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
        R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
      
	int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;

        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
	
        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI = NX-BVOIS+l+1;     
	int NXIP = (NXI+H1)%H1;  
        
        R DISTXP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
        R DISTXM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }  
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;

        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }  
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }        
        
        }    
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1;  
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;     
	int NZIP = (NZI+H3)%H3;  
	
        R DISTXP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
        R DISTXM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    
	int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;
	
        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){        
	boolc=1;  
        }  
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){        
        list_tmp.push_back(NELXP);
        }   
        
        int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
	
        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }  
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

    }
                                    
        
    }while(BVOIS<=(2*rS));


        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		 LIST_N[*it]=1;
	       }	 
	 
       } 

}

}



R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de spheres :"<< phi_inc <<" / "<<(nS*4./3*Pi*rS*rS*rS)/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par sphere :"<< int(phi_inc*N_TOT/nS) <<endl; 

}  


void gener3_sphere_alea(int nS, R eche, R buffer, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N, R* LIST_R, R* LIST_X, R* LIST_Y, R* LIST_Z)
{  

R rS=H1/(2.*eche);
buffer=R(H1)*buffer/100.+rS;

if(((R(H1)-2.*buffer)<0)||((R(H2)-2.*buffer)<0)||((R(H3)-2.*buffer)<0)){
cout<<"Spheres trop grandes !!!"<<endl;
exit(0);
}

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
       
int it=0;
while(it<nS){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * (R(H1)-2.*buffer) + buffer;
R Y_GRA = ( rand()/(double)RAND_MAX ) * (R(H2)-2.*buffer) + buffer;
R Z_GRA = ( rand()/(double)RAND_MAX ) * (R(H3)-2.*buffer) + buffer;

LIST_X[it]=X_GRA;
LIST_Y[it]=Y_GRA;
LIST_Z[it]=Z_GRA;
LIST_R[it]=rS;

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(LIST_X[it]-LIST_X[jt])*(LIST_X[it]-LIST_X[jt])+(LIST_Y[it]-LIST_Y[jt])*(LIST_Y[it]-LIST_Y[jt])+(LIST_Z[it]-LIST_Z[jt])*(LIST_Z[it]-LIST_Z[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100)) boolc=1;

	}


  if(boolc==0){
  // Test intersection voxel

  int BVOIS=0;

    do
    {  
    BVOIS++;        
    
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI = NY-BVOIS+j;  
    int NYIP = (NYI+H2)%H2;  
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI = NZ-BVOIS+k;  
	int NZIP = (NZI+H3)%H3;  

        R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
        R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
      
	int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;

        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
	
        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI = NX-BVOIS+l+1;     
	int NXIP = (NXI+H1)%H1;  
        
        R DISTXP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
        R DISTXM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }  
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;

        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }  
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }        
        
        }    
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1;  
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;     
	int NZIP = (NZI+H3)%H3;  
	
        R DISTXP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
        R DISTXM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    
	int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;
	
        if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==1)){        
	boolc=1;  
        }  
        else if((DISTXP<=(rS*rS))&&(LIST_N[NELXP]==0)){        
        list_tmp.push_back(NELXP);
        }   
        
        int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
	
        if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }  
        else if((DISTXM<=(rS*rS))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

    }
                                    
        
    }while(BVOIS<=(2*rS));



        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		 LIST_N[*it]=1;
	       }	 
	 
       } 

  }

}



R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de spheres :"<< phi_inc <<" / "<<(nS*4./3*Pi*rS*rS*rS)/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par sphere :"<< int(phi_inc*N_TOT/nS) <<endl; 

}  







void gener3_cyl_inf(int nC, R eche, R buffer, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir)
{  

R rC;

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;

R XA[nC],YA[nC],ZA[nC];

if(dir==1){
////////////// direction 1 //////////////  
rC=H2/(2.*eche);
buffer=R(H2)*buffer/100.+rC;

if(((R(H2)-2.*buffer)<0)||((R(H3)-2.*buffer)<0)){
cout<<"Fibres trop grandes !!!"<<endl;
exit(0);
}


















int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R Y_GRA = ( rand()/(double)RAND_MAX ) * (R(H2)-2.*buffer) + buffer;
R Z_GRA = ( rand()/(double)RAND_MAX ) * (R(H3)-2.*buffer) + buffer;

YA[it]=Y_GRA;
ZA[it]=Z_GRA;

int NX    = 0;
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NZI = NZ-BVOIS+j;        
    int NZIP = (NZI+H3)%H3;  


    R DISTXP = (Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    R DISTXM = (Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));


	int NELXP = NZIP*V_TOT+NYPLUSP*H1+0;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+NYMOINP*H1+0;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NYI=NY-BVOIS+j;
            int NYIP = (NYI+H2)%H2;  
	     
	    R DISTYP = (Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
	    R DISTYM = (Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));

		int NELYP = NZPLUSP*V_TOT+NYIP*H1+0;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = NZMOINP*V_TOT+NYIP*H1+0;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H1;jj++){
 		LIST_N[*it+jj]=1;
                }
	       }	 
	 
       } 

  }

}





R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H1))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 





}else if(dir==2){    
////////////// direction 2 //////////////  
rC=H1/(2.*eche);
buffer=R(H1)*buffer/100.+rC;

if(((R(H1)-2.*buffer)<0)||((R(H3)-2.*buffer)<0)){
cout<<"Fibres trop grandes !!!"<<endl;
exit(0);
}










int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * (R(H1)-2.*buffer) + buffer;
R Z_GRA = ( rand()/(double)RAND_MAX ) * (R(H3)-2.*buffer) + buffer;

XA[it]=X_GRA;
ZA[it]=Z_GRA;

int NX    = (int) floor(X_GRA);
int NY    = 0;
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NZI = NZ-BVOIS+j;        
    int NZIP = (NZI+H3)%H3;  


    R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));


	int NELXP = NZIP*V_TOT+0*H1+NXPLUSP;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+0*H1+NXMOINP;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NXI=NX-BVOIS+j;
            int NXIP = (NXI+H1)%H1;  
	     
	    R DISTYP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
	    R DISTYM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));

		int NELYP = NZPLUSP*V_TOT+0*H1+NXIP;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = NZMOINP*V_TOT+0*H1+NXIP;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H2;jj++){
 		LIST_N[*it+jj*H1]=1;
                }
	       }	 
	 
       } 

  }

}




R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H2))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 





}else if(dir==3){    
////////////// direction 3 //////////////  
rC=H1/(2.*eche);
buffer=R(H1)*buffer/100.+rC;

if(((R(H1)-2.*buffer)<0)||((R(H2)-2.*buffer)<0)){
cout<<"Fibres trop grandes !!!"<<endl;
exit(0);
}








int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * (R(H1)-2.*buffer) + buffer;
R Y_GRA = ( rand()/(double)RAND_MAX ) * (R(H2)-2.*buffer) + buffer;

XA[it]=X_GRA;
YA[it]=Y_GRA;

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = 0;

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI = NY-BVOIS+j;        
    int NYIP = (NYI+H2)%H2;  


    R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5));
    R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5));


	int NELXP = 0*V_TOT+NYIP*H1+NXPLUSP;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = 0*V_TOT+NYIP*H1+NXMOINP;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NXI=NX-BVOIS+j;
            int NXIP = (NXI+H1)%H1;  
	     
	    R DISTYP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5));
	    R DISTYM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5));

		int NELYP = 0*V_TOT+NYPLUSP*H1+NXIP;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = 0*V_TOT+NYMOINP*H1+NXIP;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H3;jj++){
 		LIST_N[*it+jj*(H1*H2)]=1;
                }
	       }	 
	 
       } 

  }

}







R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H3))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 
}






}  








void gener3_cyl_inf_per(int nC, R eche, R tol, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir)
{  

R rC;

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;

R XA[nC],YA[nC],ZA[nC];

if(dir==1){
////////////// direction 1 //////////////  
rC=H2/(2.*eche);


















int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R Y_GRA = ( rand()/(double)RAND_MAX ) * R(H2) ;
R Z_GRA = ( rand()/(double)RAND_MAX ) * R(H3) ;

YA[it]=Y_GRA;
ZA[it]=Z_GRA;

int NX    = 0;
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){


	distij=(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NZI = NZ-BVOIS+j;        
    int NZIP = (NZI+H3)%H3;  


    R DISTXP = (Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    R DISTXM = (Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));


	int NELXP = NZIP*V_TOT+NYPLUSP*H1+0;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+NYMOINP*H1+0;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NYI=NY-BVOIS+j;
            int NYIP = (NYI+H2)%H2;  
	     
	    R DISTYP = (Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
	    R DISTYM = (Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));

		int NELYP = NZPLUSP*V_TOT+NYIP*H1+0;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = NZMOINP*V_TOT+NYIP*H1+0;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H1;jj++){
 		LIST_N[*it+jj]=1;
                }
	       }	 
	 
       } 

  }

}





R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H1))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 





}else if(dir==2){    
////////////// direction 2 //////////////  
rC=H1/(2.*eche);










int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * R(H1);
R Z_GRA = ( rand()/(double)RAND_MAX ) * R(H3);

XA[it]=X_GRA;
ZA[it]=Z_GRA;

int NX    = (int) floor(X_GRA);
int NY    = 0;
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(ZA[it]-H3-ZA[jt])*(ZA[it]-H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(ZA[it]+H3-ZA[jt])*(ZA[it]+H3-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NZI = NZ-BVOIS+j;        
    int NZIP = (NZI+H3)%H3;  


    R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));
    R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Z_GRA-(NZI+0.5))*(Z_GRA-(NZI+0.5));


	int NELXP = NZIP*V_TOT+0*H1+NXPLUSP;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = NZIP*V_TOT+0*H1+NXMOINP;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NXI=NX-BVOIS+j;
            int NXIP = (NXI+H1)%H1;  
	     
	    R DISTYP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Z_GRA-(NZPLUS+0.5))*(Z_GRA-(NZPLUS+0.5));
	    R DISTYM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Z_GRA-(NZMOIN+0.5))*(Z_GRA-(NZMOIN+0.5));

		int NELYP = NZPLUSP*V_TOT+0*H1+NXIP;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = NZMOINP*V_TOT+0*H1+NXIP;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H2;jj++){
 		LIST_N[*it+jj*H1]=1;
                }
	       }	 
	 
       } 

  }

}




R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H2))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 





}else if(dir==3){    
////////////// direction 3 //////////////  
rC=H1/(2.*eche);






int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * R(H1);
R Y_GRA = ( rand()/(double)RAND_MAX ) * R(H2);

XA[it]=X_GRA;
YA[it]=Y_GRA;

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = 0;

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique
        int jt;
	R distij;
	for(jt=0;jt<it;jt++){


	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-H1-XA[jt])*(XA[it]-H1-XA[jt])+(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(YA[it]-H2-YA[jt])*(YA[it]-H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;

	distij=(XA[it]+H1-XA[jt])*(XA[it]+H1-XA[jt])+(YA[it]+H2-YA[jt])*(YA[it]+H2-YA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rC*(1.+tol/100)) boolc=1;


	}


  if(boolc==0){
  // Test intersection voxel


int BVOIS=0;

    do
    {
  
    BVOIS++;
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;   
    
    
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI = NY-BVOIS+j;        
    int NYIP = (NYI+H2)%H2;  


    R DISTXP = (X_GRA-(NXPLUS+0.5))*(X_GRA-(NXPLUS+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5));
    R DISTXM = (X_GRA-(NXMOIN+0.5))*(X_GRA-(NXMOIN+0.5))+(Y_GRA-(NYI+0.5))*(Y_GRA-(NYI+0.5));


	int NELXP = 0*V_TOT+NYIP*H1+NXPLUSP;

        if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==1)){
	boolc=1;  
        }   	
        else if((DISTXP<=(rC*rC))&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   

        int NELXM = 0*V_TOT+NYIP*H1+NXMOINP;
	
        if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==1)){
	boolc=1;  
        }    
        else if((DISTXM<=(rC*rC))&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }  


	    if((j>0)&&(j<(2*BVOIS))){
	    int NXI=NX-BVOIS+j;
            int NXIP = (NXI+H1)%H1;  
	     
	    R DISTYP = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYPLUS+0.5))*(Y_GRA-(NYPLUS+0.5));
	    R DISTYM = (X_GRA-(NXI+0.5))*(X_GRA-(NXI+0.5))+(Y_GRA-(NYMOIN+0.5))*(Y_GRA-(NYMOIN+0.5));

		int NELYP = 0*V_TOT+NYPLUSP*H1+NXIP;

		if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==1)){
		boolc=1;  
		}   	
		else if((DISTYP<=(rC*rC))&&(LIST_N[NELYP]==0)){
		list_tmp.push_back(NELYP);
		}   

		int NELYM = 0*V_TOT+NYMOINP*H1+NXIP;

		if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==1)){
		boolc=1;  
		}    
		else if((DISTYM<=(rC*rC))&&(LIST_N[NELYM]==0)){
		list_tmp.push_back(NELYM);
		}  

		                         
	    }       
    
    }
        
    }while(BVOIS<=(2*rC));




        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	it++;
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		for(int jj=1;jj<H3;jj++){
 		LIST_N[*it+jj*(H1*H2)]=1;
                }
	       }	 
	 
       } 

  }

}







R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*R(H3))/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 
}






}  























void gener3_cyl_align(int DENSI,int nC, R facf, R eche, R buffer, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir)
{

R bufferx,buffery,bufferz,lC;
if(dir==1){
lC=R(H1)/eche; 
bufferx=R(H1)*buffer/100.+R(H1)/(2.*eche);
buffery=R(H2)*buffer/100.+R(H2)/(2.*eche*facf);
bufferz=R(H3)*buffer/100.+R(H3)/(2.*eche*facf);
}else if(dir==2){
lC=R(H2)/eche; 
bufferx=R(H1)*buffer/100.+R(H1)/(2.*eche*facf);
buffery=R(H2)*buffer/100.+R(H2)/(2.*eche);
bufferz=R(H3)*buffer/100.+R(H3)/(2.*eche*facf);
}else if(dir==3){
lC=R(H3)/eche; 
bufferx=R(H1)*buffer/100.+R(H1)/(2.*eche*facf);
buffery=R(H2)*buffer/100.+R(H2)/(2.*eche*facf);
bufferz=R(H3)*buffer/100.+R(H3)/(2.*eche);
}else{
lC=R(H1)/eche; 
bufferx=R(H1)*buffer/100.+R(H1)/(2.*eche);
buffery=R(H2)*buffer/100.+R(H2)/(2.*eche*facf);
bufferz=R(H3)*buffer/100.+R(H3)/(2.*eche*facf);
}

if((R(H1)-2.*bufferx)<0||(R(H2)-2.*buffery)<0||(R(H3)-2.*bufferz)<0){
cout<<"Fibres trop longues !!!"<<endl;
exit(0);
}

R rC=lC/(2.*facf);



for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;


int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * (R(H1)-2.*bufferx) + bufferx;
R Y_GRA = ( rand()/(double)RAND_MAX ) * (R(H2)-2.*buffery) + buffery;
R Z_GRA = ( rand()/(double)RAND_MAX ) * (R(H3)-2.*bufferz) + bufferz;

R XX1,YY1,ZZ1,XX2,YY2,ZZ2;
R AA,BB,CC;

if(dir==1){
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}else if(dir==2){
XX1=X_GRA;
YY1=Y_GRA-0.5*lC;
ZZ1=Z_GRA;
XX2=X_GRA;
YY2=Y_GRA+0.5*lC;
ZZ2=Z_GRA;
AA=0.;
BB=1.;
CC=0.;
}else if(dir==3){
XX1=X_GRA;
YY1=Y_GRA;
ZZ1=Z_GRA-0.5*lC;
XX2=X_GRA;
YY2=Y_GRA;
ZZ2=Z_GRA+0.5*lC;
AA=0.;
BB=0.;
CC=1.;
}else{
cout<<"Mauvaise direction !! Dir. x par défaut"<<endl;
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique

for(int jt=0;jt<=DENSI;jt++){

R XX=XX1+jt*lC*AA/DENSI;
R YY=YY1+jt*lC*BB/DENSI;
R ZZ=ZZ1+jt*lC*CC/DENSI;

int NX    = (int) floor(XX);
int NY    = (int) floor(YY);
int NZ    = (int) floor(ZZ);
int NXP   = (NX+H1)%H1;
int NYP   = (NY+H2)%H2;
int NZP   = (NZ+H3)%H3;

bool boolxa = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NX+0.5),(NY+0.5),(NZ+0.5));
int NUMA  = NZP*V_TOT+NYP*H1+NXP;

if((boolxa==1)&&(LIST_N[NUMA]==1))
{
//cout<<it<<","<<jt<<endl;
boolc=1;  
}
else if((boolxa==1)&&(LIST_N[NUMA]==0)){
list_tmp.push_back(NUMA);
}

    //couches sup
    
for(int BVOIS=1;BVOIS<=(floor(1.5*rC));BVOIS++)
    {      
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
     
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;      
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI  = NY-BVOIS+j;       
    int NYIP = (NYI+H2)%H2; 
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI  = NZ-BVOIS+k;        
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXPLUS+0.5),(NYI+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXMOIN+0.5),(NYI+0.5),(NZI+0.5));

	int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
    int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI  = NX-BVOIS+l+1;        
        int NXIP = (NXI+H1)%H1;     
        
	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZPLUS+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZMOIN+0.5));

        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }     
        
        }  //END  
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI  = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1; 
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;      
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYPLUS+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYMOIN+0.5),(NZI+0.5));

        int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;	        
	int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
		
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;	  
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

      } //END
       

       
       
    }
 
  }


        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		 LIST_N[*it]=1;
	       }	
	it++; 	 
       } 



} //fin while it

R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*lC)/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 




} //fin subroutine
















void gener3_cyl_align_per(int DENSI,int nC, R facf, R eche, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,int dir)
{

R lC;
if(dir==1){
lC=R(H1)/eche; 
}else if(dir==2){
lC=R(H2)/eche; 
}else if(dir==3){
lC=R(H3)/eche; 
}else{
lC=R(H1)/eche; 
}

R rC=lC/(2.*facf);

for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
 
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;


int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * R(H1);
R Y_GRA = ( rand()/(double)RAND_MAX ) * R(H2);
R Z_GRA = ( rand()/(double)RAND_MAX ) * R(H3);

R XX1,YY1,ZZ1,XX2,YY2,ZZ2;
R AA,BB,CC;

if(dir==1){
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}else if(dir==2){
XX1=X_GRA;
YY1=Y_GRA-0.5*lC;
ZZ1=Z_GRA;
XX2=X_GRA;
YY2=Y_GRA+0.5*lC;
ZZ2=Z_GRA;
AA=0.;
BB=1.;
CC=0.;
}else if(dir==3){
XX1=X_GRA;
YY1=Y_GRA;
ZZ1=Z_GRA-0.5*lC;
XX2=X_GRA;
YY2=Y_GRA;
ZZ2=Z_GRA+0.5*lC;
AA=0.;
BB=0.;
CC=1.;
}else{
cout<<"Mauvaise direction !! Dir. x par défaut"<<endl;
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}

int NX    = (int) floor(X_GRA);
int NY    = (int) floor(Y_GRA);
int NZ    = (int) floor(Z_GRA);

if(X_GRA==H1){NX=H1-1;}
if(Y_GRA==H2){NY=H2-1;}
if(Z_GRA==H3){NZ=H3-1;}

int NUMA    = NZ*V_TOT+NY*H1+NX;
list_tmp.push_back(NUMA);

// Test intersection géométrique

for(int jt=0;jt<=DENSI;jt++){

R XX=XX1+jt*lC*AA/DENSI;
R YY=YY1+jt*lC*BB/DENSI;
R ZZ=ZZ1+jt*lC*CC/DENSI;

int NX    = (int) floor(XX);
int NY    = (int) floor(YY);
int NZ    = (int) floor(ZZ);
int NXP   = (NX+H1)%H1;
int NYP   = (NY+H2)%H2;
int NZP   = (NZ+H3)%H3;

bool boolxa = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NX+0.5),(NY+0.5),(NZ+0.5));
int NUMA  = NZP*V_TOT+NYP*H1+NXP;

if((boolxa==1)&&(LIST_N[NUMA]==1))
{
//cout<<it<<","<<jt<<endl;
boolc=1;  
}
else if((boolxa==1)&&(LIST_N[NUMA]==0)){
list_tmp.push_back(NUMA);
}

    //couches sup
    
for(int BVOIS=1;BVOIS<=(floor(1.5*rC));BVOIS++)
    {      
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
     
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;      
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI  = NY-BVOIS+j;       
    int NYIP = (NYI+H2)%H2; 
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI  = NZ-BVOIS+k;        
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXPLUS+0.5),(NYI+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXMOIN+0.5),(NYI+0.5),(NZI+0.5));

	int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
    int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI  = NX-BVOIS+l+1;        
        int NXIP = (NXI+H1)%H1;     
        
	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZPLUS+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZMOIN+0.5));

        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }     
        
        }  //END  
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI  = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1; 
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;      
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYPLUS+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYMOIN+0.5),(NZI+0.5));

        int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;	        
	int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
		
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;	  
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

      } //END
       

       
       
    }
 
  }


        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	       for (std::list<int>::iterator it=list_tmp.begin(); it!=list_tmp.end(); ++it){	 
		 LIST_N[*it]=1;
	       }	
	it++; 	 
       } 



} //fin while it

R phi_inc=0.;
for(it=0;it<N_TOT;it++){
if(LIST_N[it]) phi_inc=phi_inc+1.;
}
phi_inc/=N_TOT;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC*lC)/(R(H1)*R(H2)*R(H3))<<endl; 
cout<<"Nombre de voxels par fibre :"<< int(phi_inc*N_TOT/nC) <<endl; 




} //fin subroutine




void gener3_cyl_alea_per(int DENSI,int nC, R facf, R eche, int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,R * LIST_PS,R * LIST_GA,R * LIST_PH)  {
  
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
  
R Pi=3.14159265;
  
R lC=R(H1)/eche; 
R rC=lC/(2.*facf);
  
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
    
int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * R(H1);
R Y_GRA = ( rand()/(double)RAND_MAX ) * R(H2);
R Z_GRA = ( rand()/(double)RAND_MAX ) * R(H3);

R PH_REF  = rand()%360;    
R GA_REF  = ((R) (rand()%100))/R(100);
GA_REF    = acos(GA_REF);
R PS_REF  = rand()%360;  

PH_REF    = PH_REF*Pi/180; 
PS_REF    = PS_REF*Pi/180; 

R AA = cos(PH_REF)*cos(PS_REF)-cos(GA_REF)*sin(PS_REF)*sin(PH_REF);
R BB = sin(PH_REF)*cos(PS_REF)+cos(GA_REF)*sin(PS_REF)*cos(PH_REF);
R CC = sin(GA_REF)*sin(PS_REF);

R XX1=X_GRA-0.5*lC*AA;
R YY1=Y_GRA-0.5*lC*BB;
R ZZ1=Z_GRA-0.5*lC*CC;
R XX2=X_GRA+0.5*lC*AA;
R YY2=Y_GRA+0.5*lC*BB;
R ZZ2=Z_GRA+0.5*lC*CC;


for(int jt=0;jt<=DENSI;jt++){

R XX=XX1+jt*lC*AA/DENSI;
R YY=YY1+jt*lC*BB/DENSI;
R ZZ=ZZ1+jt*lC*CC/DENSI;

int NX    = (int) floor(XX);
int NY    = (int) floor(YY);
int NZ    = (int) floor(ZZ);
int NXP   = (NX+H1)%H1;
int NYP   = (NY+H2)%H2;
int NZP   = (NZ+H3)%H3;

bool boolxa = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NX+0.5),(NY+0.5),(NZ+0.5));
int NUMA  = NZP*V_TOT+NYP*H1+NXP;

if((boolxa==1)&&(LIST_N[NUMA]==1))
{
//cout<<it<<","<<jt<<endl;
boolc=1;  
}
else if((boolxa==1)&&(LIST_N[NUMA]==0)){
list_tmp.push_back(NUMA);
}

    //couches sup
    
for(int BVOIS=1;BVOIS<=(floor(1.5*rC));BVOIS++)
    {      
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
     
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;      
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI  = NY-BVOIS+j;       
    int NYIP = (NYI+H2)%H2; 
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI  = NZ-BVOIS+k;        
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXPLUS+0.5),(NYI+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXMOIN+0.5),(NYI+0.5),(NZI+0.5));

	int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
    int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI  = NX-BVOIS+l+1;        
        int NXIP = (NXI+H1)%H1;     
        
	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZPLUS+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZMOIN+0.5));

        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }     
        
        }  //END  
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI  = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1; 
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;      
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYPLUS+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYMOIN+0.5),(NZI+0.5));

        int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;	        
	int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
		
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;	  
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

      } //END
       

       
       
    }
 
  }
  
        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	 it++;
	//cout<<"dime:"<<list_tmp.size()<<endl;
       for (std::list<int>::iterator iit=list_tmp.begin(); iit!=list_tmp.end(); ++iit){	 
	 LIST_N[*iit]=1;
	 LIST_PH[*iit]=PH_REF;
	 LIST_PS[*iit]=PS_REF;
	 LIST_GA[*iit]=GA_REF;	 
//	 cout<<"PH_REF:"<<PH_REF<<","<<"PS_REF:"<<PS_REF<<","<<"GA_REF:"<<GA_REF<<endl;
       }	 
	 
       } 

}

}








void gener3_cyl_alea(int DENSI,int nC, R facf, R eche, R buffer,  int N_TOT, int V_TOT, int H1, int H2, int H3, bool * LIST_N,R * LIST_PS,R * LIST_GA,R * LIST_PH)  {
  
for(int it=0;it<N_TOT;it++){
LIST_N[it]=0;
}  
  
R Pi=3.14159265;
  
R lC=R(H1)/eche; 
R rC=lC/(2.*facf);
buffer=R(H1)*buffer/100.+lC/2.;

if(((R(H1)-2.*buffer)<0)||((R(H2)-2.*buffer)<0)||((R(H3)-2.*buffer)<0)){
cout<<"Fibres trop grandes !!!"<<endl;
exit(0);
}
  
// PARAMETRES ALEATOIRES
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    srand(tv.tv_usec) ;
    
int it=0;
while(it<nC){
  
list<int> list_tmp; 
bool boolc=0;

R X_GRA = ( rand()/(double)RAND_MAX ) * (R(H1)-2.*buffer) + buffer;
R Y_GRA = ( rand()/(double)RAND_MAX ) * (R(H2)-2.*buffer) + buffer;
R Z_GRA = ( rand()/(double)RAND_MAX ) * (R(H3)-2.*buffer) + buffer;

R PH_REF  = rand()%360;    
R GA_REF  = ((R) (rand()%100))/R(100);
GA_REF    = acos(GA_REF);
R PS_REF  = rand()%360;  

PH_REF    = PH_REF*Pi/180; 
PS_REF    = PS_REF*Pi/180; 

R AA = cos(PH_REF)*cos(PS_REF)-cos(GA_REF)*sin(PS_REF)*sin(PH_REF);
R BB = sin(PH_REF)*cos(PS_REF)+cos(GA_REF)*sin(PS_REF)*cos(PH_REF);
R CC = sin(GA_REF)*sin(PS_REF);

R XX1=X_GRA-0.5*lC*AA;
R YY1=Y_GRA-0.5*lC*BB;
R ZZ1=Z_GRA-0.5*lC*CC;
R XX2=X_GRA+0.5*lC*AA;
R YY2=Y_GRA+0.5*lC*BB;
R ZZ2=Z_GRA+0.5*lC*CC;


for(int jt=0;jt<=DENSI;jt++){

R XX=XX1+jt*lC*AA/DENSI;
R YY=YY1+jt*lC*BB/DENSI;
R ZZ=ZZ1+jt*lC*CC/DENSI;

int NX    = (int) floor(XX);
int NY    = (int) floor(YY);
int NZ    = (int) floor(ZZ);
int NXP   = (NX+H1)%H1;
int NYP   = (NY+H2)%H2;
int NZP   = (NZ+H3)%H3;

bool boolxa = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NX+0.5),(NY+0.5),(NZ+0.5));
int NUMA  = NZP*V_TOT+NYP*H1+NXP;

if((boolxa==1)&&(LIST_N[NUMA]==1))
{
//cout<<it<<","<<jt<<endl;
boolc=1;  
}
else if((boolxa==1)&&(LIST_N[NUMA]==0)){
list_tmp.push_back(NUMA);
}

    //couches sup
    
for(int BVOIS=1;BVOIS<=(floor(1.5*rC));BVOIS++)
    {      
        
    int NXPLUS = NX+(BVOIS);
    int NXMOIN = NX-(BVOIS);
    int NYPLUS = NY+(BVOIS);
    int NYMOIN = NY-(BVOIS);
    int NZPLUS = NZ+(BVOIS);
    int NZMOIN = NZ-(BVOIS);
     
    int NXPLUSP = (NXPLUS+H1)%H1;
    int NXMOINP = (NXMOIN+H1)%H1;
    int NYPLUSP = (NYPLUS+H2)%H2;
    int NYMOINP = (NYMOIN+H2)%H2;    
    int NZPLUSP = (NZPLUS+H3)%H3;
    int NZMOINP = (NZMOIN+H3)%H3;      
        
    for(int j=0;j<(2*BVOIS+1);j++){
    int NYI  = NY-BVOIS+j;       
    int NYIP = (NYI+H2)%H2; 
    
        for(int k=0;k<(2*BVOIS+1);k++){
        int NZI  = NZ-BVOIS+k;        
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXPLUS+0.5),(NYI+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXMOIN+0.5),(NYI+0.5),(NZI+0.5));

	int NELXM = NZIP*V_TOT+NYIP*H1+NXMOINP;
    int NELXP = NZIP*V_TOT+NYIP*H1+NXPLUSP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }          
        
        }

        for(int l=0;l<(2*BVOIS-1);l++){
        int NXI  = NX-BVOIS+l+1;        
        int NXIP = (NXI+H1)%H1;     
        
	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZPLUS+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYI+0.5),(NZMOIN+0.5));

        int NELXM = NZMOINP*V_TOT+NYIP*H1+NXIP;
        int NELXP = NZPLUSP*V_TOT+NYIP*H1+NXIP;
	
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }     
        
        }  //END  
    }
    
    for(int j=0;j<(2*BVOIS-1);j++){
    int NXI  = NX-BVOIS+j+1;        
    int NXIP = (NXI+H1)%H1; 
    
        for(int k=0;k<(2*BVOIS-1);k++){
        int NZI = NZ-BVOIS+k+1;      
        int NZIP = (NZI+H3)%H3;

	bool boolxp = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYPLUS+0.5),(NZI+0.5));
	bool boolxm = bincl(rC,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,(NXI+0.5),(NYMOIN+0.5),(NZI+0.5));

        int NELXP = NZIP*V_TOT+NYPLUSP*H1+NXIP;	        
	int NELXM = NZIP*V_TOT+NYMOINP*H1+NXIP;
		
        if((boolxp==1)&&(LIST_N[NELXP]==1))
	{
	//  cout<<it<<","<<jt<<endl;	  
	boolc=1;  
	}
	else if((boolxp==1)&&(LIST_N[NELXP]==0)){
        list_tmp.push_back(NELXP);
        }   
        
        if((boolxm==1)&&(LIST_N[NELXM]==1))
	{
	//  cout<<it<<","<<jt<<endl;
	boolc=1;  
	}
	else if((boolxm==1)&&(LIST_N[NELXM]==0)){
        list_tmp.push_back(NELXM);
        }   
           
        }

      } //END
       

       
       
    }
 
  }
  
        // Traitement de la liste des éléments à ajouter
       
       if(boolc==0){
	list_tmp.unique();
	 it++;
	//cout<<"dime:"<<list_tmp.size()<<endl;
       for (std::list<int>::iterator iit=list_tmp.begin(); iit!=list_tmp.end(); ++iit){	 
	 LIST_N[*iit]=1;
	 LIST_PH[*iit]=PH_REF;
	 LIST_PS[*iit]=PS_REF;
	 LIST_GA[*iit]=GA_REF;	 
//	 cout<<"PH_REF:"<<PH_REF<<","<<"PS_REF:"<<PS_REF<<","<<"GA_REF:"<<GA_REF<<endl;
       }	 
	 
       } 

}

}










