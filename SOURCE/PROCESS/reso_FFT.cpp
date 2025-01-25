#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "reso_FFT.h"
#include "FFT.h"

const R Pi=3.14159265;
const bool DEBUG = 0;

using namespace std;


void reso_fft2_ther(bool * LIST_N, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir)
{
	
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=2;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;

R KK;

if(N1%2==0) {IP1=2;} else {IP1=1;}

ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2;
	grad_co[j]=grad_co[0]+j*ND1*N2;
	flux[j]=flux[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(flux,N1,N2,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
			mf[0]+=flux[0][it];
            mf[1]+=flux[1][it];
            
			err_flux1=max(flux[0][it],err_flux1);		
			err_flux1=max(flux[1][it],err_flux1);		
			
		}
	}

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(flux,N1,N2,NDIM,-1);

if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2;

for(i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
    ic=ic1+ic2;

	R NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*flux[0][indi]+NY*flux[1][indi])*(NX*flux[0][indi]+NY*flux[1][indi]);
	errf+=((NX*flux[0][indr]+NY*flux[1][indr])*(NX*flux[0][indr]+NY*flux[1][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_fluxx+=errf;	
   }
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;
*/
}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
		}
	}
	
/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(flux,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1]<<", "<<flux[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

for(i2=0;i2<N2;i2++){

	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	    R cy=cos(2*Pi*NY/N2);
	    R sy=sin(2*Pi*NY/N2);	
	    R cyy=2*(1-cy);
	    		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	R NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
		if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  
	    R cx=cos(2*Pi*NX/N1);
	    R sx=sin(2*Pi*NX/N1);	
	    R cxx=2*(1-cx);
	    R cxy=((cx-1)*(cy-1)+sx*sy);
	    R sxy=(sx*(cy-1)-sy*(cx-1));
	    
	    R norm2=cxx+cyy; 
	    norm2*=K0;	    
	    
	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
	    grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;	
	    grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;
	    grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
	    
		grad_co[0][indr]=grad_co[0][indr]/norm2;
		grad_co[0][indi]=grad_co[0][indi]/norm2;
		grad_co[1][indr]=grad_co[1][indr]/norm2;
		grad_co[1][indi]=grad_co[1][indi]/norm2;    
		}
		else if((i1!=0)||(i2!=0))
		{
		grad_co[0][indr]=flux[0][indr]/K0;
		grad_co[0][indi]=flux[0][indi]/K0;
		grad_co[1][indr]=flux[1][indr]/K0;
		grad_co[1][indi]=flux[1][indi]/K0;
		}
		else
		{
		grad_co[0][indr]=gti[0];
		grad_co[0][indi]=0.;
		grad_co[1][indr]=gti[1];
		grad_co[1][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][2*i+j*ND1]<<", "<<grad_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(grad_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			
			err_co+=pow(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),2);
			
		}
	}
	
err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];

fft2(grad_co,N1,N2,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}

void reso_fft2_ther_flux(bool * LIST_N, R * LIST_S, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir, int dir2, R & fmin, R & fmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=2;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;

R KK;

if(N1%2==0) {IP1=2;} else {IP1=1;}

ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2;
	grad_co[j]=grad_co[0]+j*ND1*N2;
	flux[j]=flux[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(flux,N1,N2,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
			mf[0]+=flux[0][it];
            mf[1]+=flux[1][it];
            
			err_flux1=max(flux[0][it],err_flux1);		
			err_flux1=max(flux[1][it],err_flux1);		
			
		}
	}

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(flux,N1,N2,NDIM,-1);

if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2;

for(i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
    ic=ic1+ic2;

	R NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*flux[0][indi]+NY*flux[1][indi])*(NX*flux[0][indi]+NY*flux[1][indi]);
	errf+=((NX*flux[0][indr]+NY*flux[1][indr])*(NX*flux[0][indr]+NY*flux[1][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_fluxx+=errf;	
   }
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;
*/
}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
		}
	}
	
/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(flux,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1]<<", "<<flux[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

for(i2=0;i2<N2;i2++){

	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	    R cy=cos(2*Pi*NY/N2);
	    R sy=sin(2*Pi*NY/N2);	
	    R cyy=2*(1-cy);
	    		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	R NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
		if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  
	    R cx=cos(2*Pi*NX/N1);
	    R sx=sin(2*Pi*NX/N1);	
	    R cxx=2*(1-cx);
	    R cxy=((cx-1)*(cy-1)+sx*sy);
	    R sxy=(sx*(cy-1)-sy*(cx-1));
	    
	    R norm2=cxx+cyy; 
	    norm2*=K0;	    
	    
	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
	    grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;	
	    grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;
	    grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
	    
		grad_co[0][indr]=grad_co[0][indr]/norm2;
		grad_co[0][indi]=grad_co[0][indi]/norm2;
		grad_co[1][indr]=grad_co[1][indr]/norm2;
		grad_co[1][indi]=grad_co[1][indi]/norm2;    
		}
		else if((i1!=0)||(i2!=0))
		{
		grad_co[0][indr]=flux[0][indr]/K0;
		grad_co[0][indi]=flux[0][indi]/K0;
		grad_co[1][indr]=flux[1][indr]/K0;
		grad_co[1][indi]=flux[1][indi]/K0;
		}
		else
		{
		grad_co[0][indr]=gti[0];
		grad_co[0][indi]=0.;
		grad_co[1][indr]=gti[1];
		grad_co[1][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][2*i+j*ND1]<<", "<<grad_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(grad_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			
			err_co+=pow(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),2);
			
		}
	}
	
err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];

fmin=0.;
fmax=0.;

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
		
           if(dir2==1){LIST_S[itc]=KK*grad[0][it];}
           else{LIST_S[itc]=KK*grad[1][it];}
           
           fmin=(fmin<LIST_S[itc])?fmin:LIST_S[itc];
           fmax=(fmax>LIST_S[itc])?fmax:LIST_S[itc];
           
		}
	}

fft2(grad_co,N1,N2,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}

void reso_fft2_ther_temp(bool * LIST_N, R * LIST_T, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir, R & tmin, R & tmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=2;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;

R KK;

if(N1%2==0) {IP1=2;} else {IP1=1;}

ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2;
	grad_co[j]=grad_co[0]+j*ND1*N2;
	flux[j]=flux[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(flux,N1,N2,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
			mf[0]+=flux[0][it];
            mf[1]+=flux[1][it];
            
			err_flux1=max(flux[0][it],err_flux1);		
			err_flux1=max(flux[1][it],err_flux1);		
			
		}
	}

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(flux,N1,N2,NDIM,-1);

if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2;

for(i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
    ic=ic1+ic2;

	R NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*flux[0][indi]+NY*flux[1][indi])*(NX*flux[0][indi]+NY*flux[1][indi]);
	errf+=((NX*flux[0][indr]+NY*flux[1][indr])*(NX*flux[0][indr]+NY*flux[1][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_fluxx+=errf;	
   }
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;
*/
}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			
		}
	}
	
/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(flux,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1]<<", "<<flux[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

for(i2=0;i2<N2;i2++){

	R NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	    R cy=cos(2*Pi*NY/N2);
	    R sy=sin(2*Pi*NY/N2);	
	    R cyy=2*(1-cy);
	    		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	R NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
		if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  
	    R cx=cos(2*Pi*NX/N1);
	    R sx=sin(2*Pi*NX/N1);	
	    R cxx=2*(1-cx);
	    R cxy=((cx-1)*(cy-1)+sx*sy);
	    R sxy=(sx*(cy-1)-sy*(cx-1));
	    
	    R norm2=cxx+cyy; 
	    norm2*=K0;	    
	    
	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
	    grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;	
	    grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;
	    grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
	    
		grad_co[0][indr]=grad_co[0][indr]/norm2;
		grad_co[0][indi]=grad_co[0][indi]/norm2;
		grad_co[1][indr]=grad_co[1][indr]/norm2;
		grad_co[1][indi]=grad_co[1][indi]/norm2;    
		}
		else if((i1!=0)||(i2!=0))
		{
		grad_co[0][indr]=flux[0][indr]/K0;
		grad_co[0][indi]=flux[0][indi]/K0;
		grad_co[1][indr]=flux[1][indr]/K0;
		grad_co[1][indi]=flux[1][indi]/K0;
		}
		else
		{
		grad_co[0][indr]=gti[0];
		grad_co[0][indi]=0.;
		grad_co[1][indr]=gti[1];
		grad_co[1][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][2*i+j*ND1]<<", "<<grad_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(grad_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			
			err_co+=pow(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),2);
			
		}
	}
	
err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1]<<", grad_co :"<<grad_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];

tmin=0.;
tmax=0.;

	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
		
           if(dir==1){LIST_T[itc]=grad[0][it]*(i1+0.5);}
           else{LIST_T[itc]=grad[1][it]*(i2+0.5);}
           
           tmin=(tmin<LIST_T[itc])?tmin:LIST_T[itc];
           tmax=(tmax>LIST_T[itc])?tmax:LIST_T[itc];
           
		}
	}

fft2(grad_co,N1,N2,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}
void reso_fft2_elas(bool * LIST_N, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma,int dir)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=3;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;
R NX,NY;
R c11,c12,c33;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Deformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;
if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];

R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2;
	eps_co[j]=eps_co[0]+j*ND1*N2;
	sig[j]=sig[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(sig,N1,N2,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c33=C2[2][2];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c33=C1[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];	

			err_sig1=max(sig[0][it],err_sig1);
			err_sig1=max(sig[1][it],err_sig1);
			err_sig1=max(sig[2][it],err_sig1);	
			
       }
}

err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(sig,N1,N2,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2;

for(int i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(int i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
        ic=ic1+ic2;

	NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi])*(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi]);
	errf+=((NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr])*(NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_sigg+=errf;	
   }
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]));
cout<<"err_sigg :"<<err_sigg<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c33=C2pC0[2][2];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c33=C1pC0[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
       }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(sig,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1]<<", "<<sig[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R coef1=1./G0;
R coef2=K0/(G0*(K0+G0));

R ccoef1=(K0+G0)/(4*K0*G0);
R ccoef2=(-K0+G0)/(4*K0*G0);
R ccoef3=1./(2*G0);

for(i2=0;i2<N2;i2++){

	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	R NYY=NY*NY;
		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;

    R NXX=NX*NX;
	R norm2=NXX+NYY;
	norm2=1./norm2;
	R norm4=norm2*norm2;
	
	        if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  

	
	        R divs1 = NX*sig[0][indr]+NY*sig[2][indr];
	        R divs2 = NX*sig[2][indr]+NY*sig[1][indr];
	        
	        R dd = NX*divs1+NY*divs2;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[2][indi];
	        R idivs2 = NX*sig[2][indi]+NY*sig[1][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = (NX*iu2+NY*iu1)/2.;	
	        
	        /*

			R gg0=4*coef1*NX*NX/norm2-coef2*NX*NX*NX*NX/norm4;
			R gg1=-coef2*NX*NX*NY*NY/norm4;
			R gg2=4*coef1*NY*NY/norm2-coef2*NY*NY*NY*NY/norm4;
			R gg3=coef1*(NX*NX+NY*NY)/norm2-coef2*NX*NX*NY*NY/norm4;
			R gg4=2*coef1*(NX*NY)/norm2-coef2*NX*NX*NX*NY/norm4;
			R gg5=2*coef1*(NX*NY)/norm2-coef2*NX*NY*NY*NY/norm4;

		  eps_co[0][indr]=gg0*sig[0][indr]+gg1*sig[1][indr]+2*gg4*sig[2][indr];
		  eps_co[1][indr]=gg1*sig[0][indr]+gg2*sig[1][indr]+2*gg5*sig[2][indr];
		  eps_co[2][indr]=gg4*sig[0][indr]+gg5*sig[1][indr]+2*gg3*sig[2][indr]; 
		  eps_co[0][indi]=gg0*sig[0][indi]+gg1*sig[1][indi]+2*gg4*sig[2][indi];
		  eps_co[1][indi]=gg1*sig[0][indi]+gg2*sig[1][indi]+2*gg5*sig[2][indi];
		  eps_co[2][indi]=gg4*sig[0][indi]+gg5*sig[1][indi]+2*gg3*sig[2][indi]; */


		}
		else if((i1!=0)||(i2!=0))
		{
		  eps_co[0][indr]=ccoef1*sig[0][indr]+ccoef2*sig[1][indr];
		  eps_co[1][indr]=ccoef2*sig[0][indr]+ccoef1*sig[1][indr];
		  eps_co[2][indr]=ccoef3*sig[2][indr];
		  eps_co[0][indi]=ccoef1*sig[0][indi]+ccoef2*sig[1][indi];
		  eps_co[1][indi]=ccoef2*sig[0][indi]+ccoef1*sig[1][indi];
		  eps_co[2][indi]=ccoef3*sig[2][indi];
		}
		else
		{
   		  eps_co[0][indr]=di[0];
		  eps_co[0][indi]=0.;
		  eps_co[1][indr]=di[1];
		  eps_co[1][indi]=0.	;
		  eps_co[2][indr]=di[2];
		  eps_co[2][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][2*i+j*ND1]<<", "<<eps_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(eps_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
R ic11,ic12,ic33;
R eps0,eps1,eps2;

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic33=iC2_C0[2][2];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic33=iC1_C0[2][2];		
			}
			
               eps0=eps_co[0][it]-eps[0][it];
               eps1=eps_co[1][it]-eps[1][it];
               eps2=eps_co[2][it]-eps[2][it];
               
                eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*eps1;
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*eps0;
                
                eps[2][it]-=ic33*eps2;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co+=err_co1*err_co1;              
               
                
       }
}


err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];

fft2(eps_co,N1,N2,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}

void reso_fft2_elas_ctr(bool * LIST_N, R * LIST_S, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma, int dir, int dir2, R & cmin, R & cmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=3;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;
R NX,NY;
R c11,c12,c33;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Deformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;
if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];

R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2;
	eps_co[j]=eps_co[0]+j*ND1*N2;
	sig[j]=sig[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(sig,N1,N2,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c33=C2[2][2];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c33=C1[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];	

			err_sig1=max(sig[0][it],err_sig1);
			err_sig1=max(sig[1][it],err_sig1);
			err_sig1=max(sig[2][it],err_sig1);	
			
       }
}

err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(sig,N1,N2,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2;

for(int i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(int i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
        ic=ic1+ic2;

	NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi])*(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi]);
	errf+=((NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr])*(NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_sigg+=errf;	
   }
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]));
cout<<"err_sigg :"<<err_sigg<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c33=C2pC0[2][2];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c33=C1pC0[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
       }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(sig,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1]<<", "<<sig[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R coef1=1./G0;
R coef2=K0/(G0*(K0+G0));

R ccoef1=(K0+G0)/(4*K0*G0);
R ccoef2=(-K0+G0)/(4*K0*G0);
R ccoef3=1./(2*G0);

for(i2=0;i2<N2;i2++){

	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	R NYY=NY*NY;
		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;

    R NXX=NX*NX;
	R norm2=NXX+NYY;
	norm2=1./norm2;
	R norm4=norm2*norm2;
	
	        if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  

	
	        R divs1 = NX*sig[0][indr]+NY*sig[2][indr];
	        R divs2 = NX*sig[2][indr]+NY*sig[1][indr];
	        
	        R dd = NX*divs1+NY*divs2;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[2][indi];
	        R idivs2 = NX*sig[2][indi]+NY*sig[1][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = (NX*iu2+NY*iu1)/2.;	
	        
	        /*

			R gg0=4*coef1*NX*NX/norm2-coef2*NX*NX*NX*NX/norm4;
			R gg1=-coef2*NX*NX*NY*NY/norm4;
			R gg2=4*coef1*NY*NY/norm2-coef2*NY*NY*NY*NY/norm4;
			R gg3=coef1*(NX*NX+NY*NY)/norm2-coef2*NX*NX*NY*NY/norm4;
			R gg4=2*coef1*(NX*NY)/norm2-coef2*NX*NX*NX*NY/norm4;
			R gg5=2*coef1*(NX*NY)/norm2-coef2*NX*NY*NY*NY/norm4;

		  eps_co[0][indr]=gg0*sig[0][indr]+gg1*sig[1][indr]+2*gg4*sig[2][indr];
		  eps_co[1][indr]=gg1*sig[0][indr]+gg2*sig[1][indr]+2*gg5*sig[2][indr];
		  eps_co[2][indr]=gg4*sig[0][indr]+gg5*sig[1][indr]+2*gg3*sig[2][indr]; 
		  eps_co[0][indi]=gg0*sig[0][indi]+gg1*sig[1][indi]+2*gg4*sig[2][indi];
		  eps_co[1][indi]=gg1*sig[0][indi]+gg2*sig[1][indi]+2*gg5*sig[2][indi];
		  eps_co[2][indi]=gg4*sig[0][indi]+gg5*sig[1][indi]+2*gg3*sig[2][indi]; */


		}
		else if((i1!=0)||(i2!=0))
		{
		  eps_co[0][indr]=ccoef1*sig[0][indr]+ccoef2*sig[1][indr];
		  eps_co[1][indr]=ccoef2*sig[0][indr]+ccoef1*sig[1][indr];
		  eps_co[2][indr]=ccoef3*sig[2][indr];
		  eps_co[0][indi]=ccoef1*sig[0][indi]+ccoef2*sig[1][indi];
		  eps_co[1][indi]=ccoef2*sig[0][indi]+ccoef1*sig[1][indi];
		  eps_co[2][indi]=ccoef3*sig[2][indi];
		}
		else
		{
   		  eps_co[0][indr]=di[0];
		  eps_co[0][indi]=0.;
		  eps_co[1][indr]=di[1];
		  eps_co[1][indi]=0.	;
		  eps_co[2][indr]=di[2];
		  eps_co[2][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][2*i+j*ND1]<<", "<<eps_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(eps_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
R ic11,ic12,ic33;
R eps0,eps1,eps2;

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic33=iC2_C0[2][2];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic33=iC1_C0[2][2];		
			}
			
               eps0=eps_co[0][it]-eps[0][it];
               eps1=eps_co[1][it]-eps[1][it];
               eps2=eps_co[2][it]-eps[2][it];
               
                eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*eps1;
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*eps0;
                
                eps[2][it]-=ic33*eps2;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co+=err_co1*err_co1;              
               
                
       }
}


err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];

cmin =2e11;
cmax =-2e11;

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c33=C2[2][2];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c33=C1[2][2];		
			}
			
            if(dir2==1){LIST_S[itc]=c11*eps[0][it]+c12*eps[1][it]; }
            else if(dir2==2){LIST_S[itc]=c11*eps[1][it]+c12*eps[0][it];}
            else if(dir2==3){LIST_S[itc]=c33*eps[2][it];}
            
           cmin=(cmin<LIST_S[itc])?cmin:LIST_S[itc];
           cmax=(cmax>LIST_S[itc])?cmax:LIST_S[itc];
			
       }
}

fft2(eps_co,N1,N2,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}

void reso_fft2_elas_dpct(bool * LIST_N, R * LIST_U, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma, int dir, R & cmin, R & cmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;

int NDIM=3;
int N1,N2,ND1,IP1;
int i,j,i1,i2,it,itc;
N1=H1;
N2=H2;
R NX,NY;
R c11,c12,c33;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Deformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;
if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];

R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2;
	eps_co[j]=eps_co[0]+j*ND1*N2;
	sig[j]=sig[0]+j*ND1*N2;
	}		
	for(i=0;i<ND1*N2;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

// Initialization
fft2(sig,N1,N2,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-18){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c33=C2[2][2];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c33=C1[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];	

			err_sig1=max(sig[0][it],err_sig1);
			err_sig1=max(sig[1][it],err_sig1);
			err_sig1=max(sig[2][it],err_sig1);	
			
       }
}

err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft2(sig,N1,N2,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2;

for(int i2=0;i2<N2;i2++){
	ic2=0;
	if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
	
	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
		
   for(int i1=0;i1<(N1/2)+1;i1++){
	ic1=0;
	if((i1==N1/2)&&(N1%2==0)) ic1=1;	
        ic=ic1+ic2;

	NX=R(i1);
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;
	
	R errf=(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi])*(NX*sig[0][indi]+NY*sig[1][indi]+(NX+NY)*sig[2][indi]);
	errf+=((NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr])*(NX*sig[0][indr]+NY*sig[1][indr]+(NX+NY)*sig[2][indr]));
	
	if(ic==2){errf=errf*0.25;}
	else if(ic==1){errf=errf*0.5;}	
	if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
	
   err_sigg+=errf;	
   }
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]));
cout<<"err_sigg :"<<err_sigg<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c33=C2pC0[2][2];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c33=C1pC0[2][2];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*eps[1][it];
			sig[1][it]=c11*eps[1][it]+c12*eps[0][it];
			sig[2][it]=c33*eps[2][it];
       }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft2(sig,N1,N2,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1]<<", "<<sig[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R coef1=1./G0;
R coef2=K0/(G0*(K0+G0));

R ccoef1=(K0+G0)/(4*K0*G0);
R ccoef2=(-K0+G0)/(4*K0*G0);
R ccoef3=1./(2*G0);

for(i2=0;i2<N2;i2++){

	NY=R(i2);
	if(NY>N2/2){NY=NY-N2;}
	R NYY=NY*NY;
		
   for(i1=0;i1<(N1/2)+1;i1++){
	  	   
	NX=R(i1);
	
	int indr=2*i1+i2*ND1;
	int indi=2*i1+1+i2*ND1;

    R NXX=NX*NX;
	R norm2=NXX+NYY;
	norm2=1./norm2;
	R norm4=norm2*norm2;
	
	        if((i1!=N1/2)&&(i2!=N2/2)&&((i1!=0)||(i2!=0))){  

	
	        R divs1 = NX*sig[0][indr]+NY*sig[2][indr];
	        R divs2 = NX*sig[2][indr]+NY*sig[1][indr];
	        
	        R dd = NX*divs1+NY*divs2;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[2][indi];
	        R idivs2 = NX*sig[2][indi]+NY*sig[1][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = (NX*iu2+NY*iu1)/2.;	
	        
	        /*

			R gg0=4*coef1*NX*NX/norm2-coef2*NX*NX*NX*NX/norm4;
			R gg1=-coef2*NX*NX*NY*NY/norm4;
			R gg2=4*coef1*NY*NY/norm2-coef2*NY*NY*NY*NY/norm4;
			R gg3=coef1*(NX*NX+NY*NY)/norm2-coef2*NX*NX*NY*NY/norm4;
			R gg4=2*coef1*(NX*NY)/norm2-coef2*NX*NX*NX*NY/norm4;
			R gg5=2*coef1*(NX*NY)/norm2-coef2*NX*NY*NY*NY/norm4;

		  eps_co[0][indr]=gg0*sig[0][indr]+gg1*sig[1][indr]+2*gg4*sig[2][indr];
		  eps_co[1][indr]=gg1*sig[0][indr]+gg2*sig[1][indr]+2*gg5*sig[2][indr];
		  eps_co[2][indr]=gg4*sig[0][indr]+gg5*sig[1][indr]+2*gg3*sig[2][indr]; 
		  eps_co[0][indi]=gg0*sig[0][indi]+gg1*sig[1][indi]+2*gg4*sig[2][indi];
		  eps_co[1][indi]=gg1*sig[0][indi]+gg2*sig[1][indi]+2*gg5*sig[2][indi];
		  eps_co[2][indi]=gg4*sig[0][indi]+gg5*sig[1][indi]+2*gg3*sig[2][indi]; */


		}
		else if((i1!=0)||(i2!=0))
		{
		  eps_co[0][indr]=ccoef1*sig[0][indr]+ccoef2*sig[1][indr];
		  eps_co[1][indr]=ccoef2*sig[0][indr]+ccoef1*sig[1][indr];
		  eps_co[2][indr]=ccoef3*sig[2][indr];
		  eps_co[0][indi]=ccoef1*sig[0][indi]+ccoef2*sig[1][indi];
		  eps_co[1][indi]=ccoef2*sig[0][indi]+ccoef1*sig[1][indi];
		  eps_co[2][indi]=ccoef3*sig[2][indi];
		}
		else
		{
   		  eps_co[0][indr]=di[0];
		  eps_co[0][indi]=0.;
		  eps_co[1][indr]=di[1];
		  eps_co[1][indi]=0.	;
		  eps_co[2][indr]=di[2];
		  eps_co[2][indi]=0.	;
		} 	

   }
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][2*i+j*ND1]<<", "<<eps_co[jt][2*i+1+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft2(eps_co,N1,N2,NDIM,1); // re initialization ???

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
R ic11,ic12,ic33;
R eps0,eps1,eps2;

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic33=iC2_C0[2][2];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic33=iC1_C0[2][2];		
			}
			
               eps0=eps_co[0][it]-eps[0][it];
               eps1=eps_co[1][it]-eps[1][it];
               eps2=eps_co[2][it]-eps[2][it];
               
                eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*eps1;
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*eps0;
                
                eps[2][it]-=ic33*eps2;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co+=err_co1*err_co1;              
               
                
       }
}


err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1]<<", eps_co :"<<eps_co[jt][i+j*ND1]<<endl;
    }
  }
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];

cmin =2e11;
cmax =-2e11;

for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1;  
        it=i1+i2*ND1;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c33=C2[2][2];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c33=C1[2][2];		
			}
			
            if(dir==1){LIST_U[itc]=eps[0][it]*(i1+0.5);}
            else if(dir==2){LIST_U[itc]=eps[1][it]*(i2+0.5);}
            
           cmin=(cmin<LIST_U[itc])?cmin:LIST_U[itc];
           cmax=(cmax>LIST_U[itc])?cmax:LIST_U[itc];
			
       }
}

fft2(eps_co,N1,N2,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}

void reso_fft3_ther(bool * LIST_N, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=3;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R KK,NX,NY,NZ;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2*N3];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2*N3];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2*N3;
	grad_co[j]=grad_co[0]+j*ND1*N2*N3;
	flux[j]=flux[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(flux,N1,N2,N3,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-15){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];

			mf[0]+=flux[0][it];
			mf[1]+=flux[1][it];
			mf[2]+=flux[2][it];

			err_flux1=max(flux[0][it],err_flux1);
			err_flux1=max(flux[1][it],err_flux1);
			err_flux1=max(flux[2][it],err_flux1);	
			
       }
    }
}	

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(flux,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi])*(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi]);
		errf+=((NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr])*(NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_fluxx+=errf;	
	   }
	}
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]+flux[2][0]*flux[2][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];
			
       }
    }
}	

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", flux :"<<flux[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(flux,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1+k*ND1*N2]<<", "<<flux[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R cx,cy,cz,sx,sy,sz,cxx,cxy,cxz,sxy,sxz,cyy,cyz,syz,czz,norm2;
int indi,indr;

for(i3=0;i3<N3;i3++){
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	cz=cos(2*Pi*NZ/N3);
	sz=sin(2*Pi*NZ/N3);		
	czz=2*(1-cz);
	
	for(i2=0;i2<N2;i2++){

		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
		cy=cos(2*Pi*NY/N2);
		sy=sin(2*Pi*NY/N2);	
		cyy=2*(1-cy);	
	    cyz=((cz-1)*(cy-1)+sy*sz); 	
	    syz=(sy*(cz-1)-sz*(cy-1));   
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
			cx=cos(2*Pi*NX/N1);
			sx=sin(2*Pi*NX/N1);	
			cxx=2*(1-cx);
							    
			norm2=2*(1-cx)+2*(1-cy)+2*(1-cz);  
			norm2*=K0;           
            cxy=((cx-1)*(cy-1)+sx*sy);
            cxz=((cx-1)*(cz-1)+sx*sz);              
            sxy=(sx*(cy-1)-sy*(cx-1));
            sxz=(sx*(cz-1)-sz*(cx-1));                               

	  	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
			grad_co[0][indr]+=cxz*flux[2][indr]-flux[2][indi]*sxz;
		
			grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;			
            grad_co[0][indi]+=cxz*flux[2][indi]+flux[2][indr]*sxz;

            grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;  
            grad_co[1][indr]+=cyz*flux[2][indr]-flux[2][indi]*syz;

            grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
            grad_co[1][indi]+=cyz*flux[2][indi]+flux[2][indr]*syz;

            grad_co[2][indr]=czz*flux[2][indr]+cxz*flux[0][indr]+flux[0][indi]*sxz;  
            grad_co[2][indr]+=cyz*flux[1][indr]+flux[1][indi]*syz;

            grad_co[2][indi]=czz*flux[2][indi]+cxz*flux[0][indi]-flux[0][indr]*sxz;  
            grad_co[2][indi]+=cyz*flux[1][indi]-flux[1][indr]*syz;

			
			grad_co[0][indr]=grad_co[0][indr]/norm2;
			grad_co[0][indi]=grad_co[0][indi]/norm2;
			grad_co[1][indr]=grad_co[1][indr]/norm2;
			grad_co[1][indi]=grad_co[1][indi]/norm2;    
			grad_co[2][indr]=grad_co[2][indr]/norm2;
			grad_co[2][indi]=grad_co[2][indi]/norm2;   			
			
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
			grad_co[0][indr]=flux[0][indr]/K0;
			grad_co[0][indi]=flux[0][indi]/K0;
			grad_co[1][indr]=flux[1][indr]/K0;
			grad_co[1][indi]=flux[1][indi]/K0;
			grad_co[2][indr]=flux[2][indr]/K0;
			grad_co[2][indi]=flux[2][indi]/K0;			
			}
			else
			{
			grad_co[0][indr]=gti[0];
			grad_co[0][indi]=0.;
			grad_co[1][indr]=gti[1];
			grad_co[1][indi]=0.	;
			grad_co[2][indr]=gti[2];
			grad_co[2][indi]=0.	;			
			} 	

	   }
	}

}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<grad_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(grad_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			grad[2][it]-=KK*(-grad[2][it]+grad_co[2][it]);
			
			err_co+=pow(max(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),abs(grad[2][it]-grad_co[2][it])),2);
			
       }
    }
}	

err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]+gti[2]*gti[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];
flu[2]=mf[2];

fft3(grad_co,N1,N2,N3,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}


void reso_fft3_ther_temp(bool * LIST_N, R * LIST_T, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir,R & tmin, R & tmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=3;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R KK,NX,NY,NZ;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2*N3];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2*N3];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2*N3;
	grad_co[j]=grad_co[0]+j*ND1*N2*N3;
	flux[j]=flux[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(flux,N1,N2,N3,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-15){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];

			mf[0]+=flux[0][it];
			mf[1]+=flux[1][it];
			mf[2]+=flux[2][it];

			err_flux1=max(flux[0][it],err_flux1);
			err_flux1=max(flux[1][it],err_flux1);
			err_flux1=max(flux[2][it],err_flux1);	
			
       }
    }
}	

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(flux,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi])*(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi]);
		errf+=((NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr])*(NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_fluxx+=errf;	
	   }
	}
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]+flux[2][0]*flux[2][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];
			
       }
    }
}	

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", flux :"<<flux[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(flux,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1+k*ND1*N2]<<", "<<flux[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R cx,cy,cz,sx,sy,sz,cxx,cxy,cxz,sxy,sxz,cyy,cyz,syz,czz,norm2;
int indi,indr;

for(i3=0;i3<N3;i3++){
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	cz=cos(2*Pi*NZ/N3);
	sz=sin(2*Pi*NZ/N3);		
	czz=2*(1-cz);
	
	for(i2=0;i2<N2;i2++){

		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
		cy=cos(2*Pi*NY/N2);
		sy=sin(2*Pi*NY/N2);	
		cyy=2*(1-cy);	
	    cyz=((cz-1)*(cy-1)+sy*sz); 	
	    syz=(sy*(cz-1)-sz*(cy-1));   
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
			cx=cos(2*Pi*NX/N1);
			sx=sin(2*Pi*NX/N1);	
			cxx=2*(1-cx);
							    
			norm2=2*(1-cx)+2*(1-cy)+2*(1-cz);  
			norm2*=K0;           
            cxy=((cx-1)*(cy-1)+sx*sy);
            cxz=((cx-1)*(cz-1)+sx*sz);              
            sxy=(sx*(cy-1)-sy*(cx-1));
            sxz=(sx*(cz-1)-sz*(cx-1));                               

	  	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
			grad_co[0][indr]+=cxz*flux[2][indr]-flux[2][indi]*sxz;
		
			grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;			
            grad_co[0][indi]+=cxz*flux[2][indi]+flux[2][indr]*sxz;

            grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;  
            grad_co[1][indr]+=cyz*flux[2][indr]-flux[2][indi]*syz;

            grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
            grad_co[1][indi]+=cyz*flux[2][indi]+flux[2][indr]*syz;

            grad_co[2][indr]=czz*flux[2][indr]+cxz*flux[0][indr]+flux[0][indi]*sxz;  
            grad_co[2][indr]+=cyz*flux[1][indr]+flux[1][indi]*syz;

            grad_co[2][indi]=czz*flux[2][indi]+cxz*flux[0][indi]-flux[0][indr]*sxz;  
            grad_co[2][indi]+=cyz*flux[1][indi]-flux[1][indr]*syz;

			
			grad_co[0][indr]=grad_co[0][indr]/norm2;
			grad_co[0][indi]=grad_co[0][indi]/norm2;
			grad_co[1][indr]=grad_co[1][indr]/norm2;
			grad_co[1][indi]=grad_co[1][indi]/norm2;    
			grad_co[2][indr]=grad_co[2][indr]/norm2;
			grad_co[2][indi]=grad_co[2][indi]/norm2;   			
			
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
			grad_co[0][indr]=flux[0][indr]/K0;
			grad_co[0][indi]=flux[0][indi]/K0;
			grad_co[1][indr]=flux[1][indr]/K0;
			grad_co[1][indi]=flux[1][indi]/K0;
			grad_co[2][indr]=flux[2][indr]/K0;
			grad_co[2][indi]=flux[2][indi]/K0;			
			}
			else
			{
			grad_co[0][indr]=gti[0];
			grad_co[0][indi]=0.;
			grad_co[1][indr]=gti[1];
			grad_co[1][indi]=0.	;
			grad_co[2][indr]=gti[2];
			grad_co[2][indi]=0.	;			
			} 	

	   }
	}

}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<grad_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(grad_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			grad[2][it]-=KK*(-grad[2][it]+grad_co[2][it]);
			
			err_co+=pow(max(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),abs(grad[2][it]-grad_co[2][it])),2);
			
       }
    }
}	

err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]+gti[2]*gti[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];
flu[2]=mf[2];

tmin=0.;
tmax=0.;


for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}

            if(dir==1){LIST_T[itc]=grad[0][it]*(i1+0.5);}
            else if(dir==2){LIST_T[itc]=grad[1][it]*(i2+0.5);}
            else if(dir==3){LIST_T[itc]=grad[2][it]*(i3+0.5);}
		
		   tmin=(tmin<LIST_T[itc])?tmin:LIST_T[itc];
           tmax=(tmax>LIST_T[itc])?tmax:LIST_T[itc];
           	
       }
    }
}	

fft3(grad_co,N1,N2,N3,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}


void reso_fft3_ther_flux(bool * LIST_N, R * LIST_S, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir, int dir2, R & fmin, R & fmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=3;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R KK,NX,NY,NZ;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//gradient de temp imposé
R gti[NDIM];
for(int j=0;j<NDIM;j++){
gti[j]=0.;	
}
gti[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R mf[NDIM];

R ** grad  = new R * [NDIM];
grad[0] = new R [NDIM*ND1*N2*N3];
R ** grad_co  = new R * [NDIM];
grad_co[0] = new R [NDIM*ND1*N2*N3];
R ** flux  = new R * [NDIM];
flux[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	grad[j]=grad[0]+j*ND1*N2*N3;
	grad_co[j]=grad_co[0]+j*ND1*N2*N3;
	flux[j]=flux[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		grad[j][i]=gti[j];
	    }
	    else
	    {
		grad[j][i]=0.;	
     	}
		grad_co[j][i]=0.;
		flux[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(flux,N1,N2,N3,NDIM,0); //grad ? grad_co ?	

R err_eq=1.;
R err_co=1.;
R err_flux0=1.;
R err_flux1;
int n_it=-1;

while((max(err_eq,err_co)>1e-4)&&(n_it<2000)){
n_it++;

if(err_co<1e-15){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_flux1=0.;
for(j=0;j<NDIM;j++){
mf[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];

			mf[0]+=flux[0][it];
			mf[1]+=flux[1][it];
			mf[2]+=flux[2][it];

			err_flux1=max(flux[0][it],err_flux1);
			err_flux1=max(flux[1][it],err_flux1);
			err_flux1=max(flux[2][it],err_flux1);	
			
       }
    }
}	

err_eq=abs(err_flux1-err_flux0)/err_flux0;
err_flux0=err_flux1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(flux,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/flux/gradient de température

R err_fluxx=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi])*(NX*flux[0][indi]+NY*flux[1][indi]+NZ*flux[2][indi]);
		errf+=((NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr])*(NX*flux[0][indr]+NY*flux[1][indr]+NZ*flux[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_fluxx+=errf;	
	   }
	}
}
err_fluxx=sqrt(err_fluxx/(flux[0][0]*flux[0][0]+flux[1][0]*flux[1][0]+flux[2][0]*flux[2][0]));
cout<<"err_fluxx :"<<err_fluxx<<endl;*/

}

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2pK0;
			}else{
			KK=K1pK0;		
			}
			
			flux[0][it]=KK*grad[0][it];
			flux[1][it]=KK*grad[1][it];
			flux[2][it]=KK*grad[2][it];
			
       }
    }
}	

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", flux :"<<flux[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(flux,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", flux :"<<flux[jt][2*i+j*ND1+k*ND1*N2]<<", "<<flux[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

R cx,cy,cz,sx,sy,sz,cxx,cxy,cxz,sxy,sxz,cyy,cyz,syz,czz,norm2;
int indi,indr;

for(i3=0;i3<N3;i3++){
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	cz=cos(2*Pi*NZ/N3);
	sz=sin(2*Pi*NZ/N3);		
	czz=2*(1-cz);
	
	for(i2=0;i2<N2;i2++){

		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
		cy=cos(2*Pi*NY/N2);
		sy=sin(2*Pi*NY/N2);	
		cyy=2*(1-cy);	
	    cyz=((cz-1)*(cy-1)+sy*sz); 	
	    syz=(sy*(cz-1)-sz*(cy-1));   
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
			cx=cos(2*Pi*NX/N1);
			sx=sin(2*Pi*NX/N1);	
			cxx=2*(1-cx);
							    
			norm2=2*(1-cx)+2*(1-cy)+2*(1-cz);  
			norm2*=K0;           
            cxy=((cx-1)*(cy-1)+sx*sy);
            cxz=((cx-1)*(cz-1)+sx*sz);              
            sxy=(sx*(cy-1)-sy*(cx-1));
            sxz=(sx*(cz-1)-sz*(cx-1));                               

	  	    grad_co[0][indr]=cxx*flux[0][indr]+cxy*flux[1][indr]-flux[1][indi]*sxy;
			grad_co[0][indr]+=cxz*flux[2][indr]-flux[2][indi]*sxz;
		
			grad_co[0][indi]=cxx*flux[0][indi]+cxy*flux[1][indi]+flux[1][indr]*sxy;			
            grad_co[0][indi]+=cxz*flux[2][indi]+flux[2][indr]*sxz;

            grad_co[1][indr]=cyy*flux[1][indr]+cxy*flux[0][indr]+flux[0][indi]*sxy;  
            grad_co[1][indr]+=cyz*flux[2][indr]-flux[2][indi]*syz;

            grad_co[1][indi]=cyy*flux[1][indi]+cxy*flux[0][indi]-flux[0][indr]*sxy;  
            grad_co[1][indi]+=cyz*flux[2][indi]+flux[2][indr]*syz;

            grad_co[2][indr]=czz*flux[2][indr]+cxz*flux[0][indr]+flux[0][indi]*sxz;  
            grad_co[2][indr]+=cyz*flux[1][indr]+flux[1][indi]*syz;

            grad_co[2][indi]=czz*flux[2][indi]+cxz*flux[0][indi]-flux[0][indr]*sxz;  
            grad_co[2][indi]+=cyz*flux[1][indi]-flux[1][indr]*syz;

			
			grad_co[0][indr]=grad_co[0][indr]/norm2;
			grad_co[0][indi]=grad_co[0][indi]/norm2;
			grad_co[1][indr]=grad_co[1][indr]/norm2;
			grad_co[1][indi]=grad_co[1][indi]/norm2;    
			grad_co[2][indr]=grad_co[2][indr]/norm2;
			grad_co[2][indi]=grad_co[2][indi]/norm2;   			
			
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
			grad_co[0][indr]=flux[0][indr]/K0;
			grad_co[0][indi]=flux[0][indi]/K0;
			grad_co[1][indr]=flux[1][indr]/K0;
			grad_co[1][indi]=flux[1][indi]/K0;
			grad_co[2][indr]=flux[2][indr]/K0;
			grad_co[2][indi]=flux[2][indi]/K0;			
			}
			else
			{
			grad_co[0][indr]=gti[0];
			grad_co[0][indi]=0.;
			grad_co[1][indr]=gti[1];
			grad_co[1][indi]=0.	;
			grad_co[2][indr]=gti[2];
			grad_co[2][indi]=0.	;			
			} 	

	   }
	}

}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<grad_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(grad_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=iK2;
			}else{
			KK=iK1;		
			}
			
            grad[0][it]-=KK*(-grad[0][it]+grad_co[0][it]);
			grad[1][it]-=KK*(-grad[1][it]+grad_co[1][it]);
			grad[2][it]-=KK*(-grad[2][it]+grad_co[2][it]);
			
			err_co+=pow(max(max(abs(grad[0][it]-grad_co[0][it]),abs(grad[1][it]-grad_co[1][it])),abs(grad[2][it]-grad_co[2][it])),2);
			
       }
    }
}	

err_co/=N_TOT;
err_co=err_co/sqrt(gti[0]*gti[0]+gti[1]*gti[1]+gti[2]*gti[2]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", grad :"<<grad[jt][i+j*ND1+k*ND1*N2]<<", grad_co :"<<grad_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

flu[0]=mf[0];
flu[1]=mf[1];
flu[2]=mf[2];

fmin=0.;
fmax=0.;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			KK=K2;
			}else{
			KK=K1;		
			}

            if(dir2==1){LIST_S[itc]=KK*grad[0][it];}
            else if(dir2==2){LIST_S[itc]=KK*grad[1][it];}
            else if(dir2==3){LIST_S[itc]=KK*grad[2][it];}
 
           fmin=(fmin<LIST_S[itc])?fmin:LIST_S[itc];
           fmax=(fmax>LIST_S[itc])?fmax:LIST_S[itc];
                      
			
       }
    }
}	

fft3(grad_co,N1,N2,N3,NDIM,9);
delete [] grad[0];
delete [] grad_co[0];
delete [] flux[0];
delete [] grad;
delete [] grad_co;
delete [] flux;
}

return;
}

void reso_fft3_elas(bool * LIST_N, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma,int dir)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=6;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R NX,NY,NZ;
//R CC[6][6];
R c11,c12,c44;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Déformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];
			
R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2*N3];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2*N3];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2*N3;
	eps_co[j]=eps_co[0]+j*ND1*N2*N3;
	sig[j]=sig[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(sig,N1,N2,N3,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-3)&&(n_it<2000)){
n_it++;

if(err_co<1e-4){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c44=C2[3][3];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c44=C1[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c11*eps[1][it]+c12*(eps[2][it]+eps[0][it]);
			sig[2][it]=c11*eps[2][it]+c12*(eps[0][it]+eps[1][it]);
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];	
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];
			ms[3]+=sig[3][it];
			ms[4]+=sig[4][it];
			ms[5]+=sig[5][it];		

			err_sig1=max(abs(sig[0][it]),err_sig1);
			err_sig1=max(abs(sig[1][it]),err_sig1);
			err_sig1=max(abs(sig[2][it]),err_sig1);
			err_sig1=max(abs(sig[3][it]),err_sig1);
			err_sig1=max(abs(sig[4][it]),err_sig1);
			err_sig1=max(abs(sig[5][it]),err_sig1);		
			
       }
    }
}		

cout<<"err_sig0 :"<<err_sig0<<", "<<"err_sig1 :"<<err_sig1<<endl;
err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(sig,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi])*(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi]);
		errf+=((NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr])*(NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr]));
		
		errf+=(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi])*(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi]);
		errf+=((NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr])*(NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr]));
		
		errf+=(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi])*(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi]);
		errf+=((NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr])*(NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_sigg+=errf;	
	   }
	}
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]+sig[3][0]*sig[3][0]+sig[4][0]*sig[4][0]+sig[5][0]*sig[5][0]));
cout<<"err_sigg :"<<err_sigg<<endl;

*/

}

//R t4i=CPUtime();

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c44=C2pC0[3][3];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c44=C1pC0[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c12*(eps[0][it]+eps[2][it])+c11*eps[1][it];
			sig[2][it]=c12*(eps[0][it]+eps[1][it])+c11*eps[2][it];
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];		
			
       }
    }
}		

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", sig :"<<sig[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 4 :"<<(CPUtime()-t4i)<<endl;
//R t5i=CPUtime();

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(sig,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1+k*ND1*N2]<<", "<<sig[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 5 :"<<(CPUtime()-t5i)<<endl;
//R t6i=CPUtime();

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

	
int indi,indr;
R coef3=(3.*K0+G0)/(9.*G0*K0);
R coef4=(-3.*K0+2.*G0)/(12.*G0*K0);
R coef5=1./(2.*G0);		
/*		
for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    R NYZ=NY*NZ;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R NXY=NX*NY;
		R NXZ=NX*NZ;
		R norm2=NXX+NYY+NZZ;
		R norm4=norm2*norm2;
		R coef1=1./(4*G0*norm2);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0)*norm4);

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
			//0 -> G11			//1 -> G12			//2 -> G13
			//3 -> G14			//4 -> G15			//5 -> G16
			//6 -> G22			//7 -> G23			//8 -> G24
			//9 -> G25			//10-> G26			//11-> G33
			//12-> G34			//13-> G35			//14-> G36
			//15-> G44			//16-> G45			//17-> G46
			//18-> G55			//19-> G56			//20-> G66
	
			R gg0=(4.*coef1-coef2*NXX)*NXX;
			R gg1=-coef2*NXX*NYY;
			R gg2=-coef2*NXX*NZZ;
			R gg3=-coef2*NXX*NYZ;
			R gg4=(2.*coef1-coef2*NXX)*NXZ;
			R gg5=(2.*coef1-coef2*NXX)*NXY;
			R gg6=(4.*coef1-coef2*NYY)*NYY;
			R gg7=-coef2*NYY*NZZ;
            R gg8=(2.*coef1-coef2*NYY)*NYZ;
			R gg9=-coef2*NYY*NXZ;
            R gg10=(2.*coef1-coef2*NYY)*NXY;
            R gg11=(4.*coef1-coef2*NZZ)*NZZ;
            R gg12=(2.*coef1-coef2*NZZ)*NYZ;
			R gg13=(2.*coef1-coef2*NZZ)*NXZ;
		    R gg14=-coef2*NXY*NZZ;
            R gg15=coef1*(NYY+NZZ)-coef2*NYY*NZZ;
            R gg16=(coef1-coef2*NZZ)*NXY;
            R gg17=(coef1-coef2*NYY)*NXZ;
			R gg18=coef1*(NXX+NZZ)-coef2*NXX*NZZ;
			R gg19=(coef1-coef2*NXX)*NYZ;
			R gg20=coef1*(NXX+NYY)-coef2*NXX*NYY;

              eps_co[0][indr] =gg0*sig[0][indr]+gg1*sig[1][indr]+gg2*sig[2][indr];
			  eps_co[0][indr]+=2*(gg3*sig[3][indr]+gg4*sig[4][indr]+gg5*sig[5][indr]);
			  
			  eps_co[1][indr] =gg1*sig[0][indr]+gg6*sig[1][indr]+gg7*sig[2][indr];
			  eps_co[1][indr]+=2*(gg8*sig[3][indr]+gg9*sig[4][indr]+gg10*sig[5][indr]);
			  
			  eps_co[2][indr] =gg2*sig[0][indr]+gg7*sig[1][indr]+gg11*sig[2][indr];
			  eps_co[2][indr]+=2*(gg12*sig[3][indr]+gg13*sig[4][indr]+gg14*sig[5][indr]);
			  
			  eps_co[3][indr] =gg3*sig[0][indr]+gg8*sig[1][indr]+gg12*sig[2][indr];
			  eps_co[3][indr]+=2*(gg15*sig[3][indr]+gg16*sig[4][indr]+gg17*sig[5][indr]);
			  
			  eps_co[4][indr] =gg4*sig[0][indr]+gg9*sig[1][indr]+gg13*sig[2][indr];
			  eps_co[4][indr]+=2*(gg16*sig[3][indr]+gg18*sig[4][indr]+gg19*sig[5][indr]);
			   
			  eps_co[5][indr] =gg5*sig[0][indr]+gg10*sig[1][indr]+gg14*sig[2][indr];
			  eps_co[5][indr]+=2*(gg17*sig[3][indr]+gg19*sig[4][indr]+gg20*sig[5][indr]);                             
 
			  eps_co[0][indi] =gg0*sig[0][indi]+gg1*sig[1][indi]+gg2*sig[2][indi];
			  eps_co[0][indi]+=2*(gg3*sig[3][indi]+gg4*sig[4][indi]+gg5*sig[5][indi]);
			  
			  eps_co[1][indi] =gg1*sig[0][indi]+gg6*sig[1][indi]+gg7*sig[2][indi];
			  eps_co[1][indi]+=2*(gg8*sig[3][indi]+gg9*sig[4][indi]+gg10*sig[5][indi]);
			  
			  eps_co[2][indi] =gg2*sig[0][indi]+gg7*sig[1][indi]+gg11*sig[2][indi];
			  eps_co[2][indi]+=2*(gg12*sig[3][indi]+gg13*sig[4][indi]+gg14*sig[5][indi]);
			  
			  eps_co[3][indi] =gg3*sig[0][indi]+gg8*sig[1][indi]+gg12*sig[2][indi];
			  eps_co[3][indi]+=2*(gg15*sig[3][indi]+gg16*sig[4][indi]+gg17*sig[5][indi]);
			  
			  eps_co[4][indi] =gg4*sig[0][indi]+gg9*sig[1][indi]+gg13*sig[2][indi];
			  eps_co[4][indi]+=2*(gg16*sig[3][indi]+gg18*sig[4][indi]+gg19*sig[5][indi]);
			   
			  eps_co[5][indi] =gg5*sig[0][indi]+gg10*sig[1][indi]+gg14*sig[2][indi];
			  eps_co[5][indi]+=2*(gg17*sig[3][indi]+gg19*sig[4][indi]+gg20*sig[5][indi]);   
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

}
*/
 

for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R norm2=NXX+NYY+NZZ;
		norm2=1./norm2;
		R norm4=norm2*norm2;
		
		R coef1=1./(G0);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0));

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
	        R divs1 = NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr];
	        R divs2 = NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr];
	        R divs3 = NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr];	
	        
	        R dd = NX*divs1+NY*divs2+NZ*divs3;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        R u3 = coef1*norm2*divs3-coef2*norm4*NZ*dd;	
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = NZ*u3;
	        eps_co[3][indr] = (NY*u3+NZ*u2)/2.;
	        eps_co[4][indr] = (NX*u3+NZ*u1)/2.;
	        eps_co[5][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi];
	        R idivs2 = NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi];
	        R idivs3 = NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2+NZ*idivs3;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        R iu3 = coef1*norm2*idivs3-coef2*norm4*NZ*idd;	
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = NZ*iu3;
	        eps_co[3][indi] = (NY*iu3+NZ*iu2)/2.;
	        eps_co[4][indi] = (NX*iu3+NZ*iu1)/2.;
	        eps_co[5][indi] = (NX*iu2+NY*iu1)/2.;	        	
	        	
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

} 


/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<eps_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 6 :"<<(CPUtime()-t6i)<<endl;
//R t7i=CPUtime();

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(eps_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 7 :"<<(CPUtime()-t7i)<<endl;
//R t8i=CPUtime();

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
//R iCX_C0[6][6];
R ic11,ic12,ic44;
R eps0,eps1,eps2,eps3,eps4,eps5;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic44=iC2_C0[3][3];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic44=iC1_C0[3][3];		
			}
			
                eps0=eps_co[0][it]-eps[0][it];
                eps1=eps_co[1][it]-eps[1][it];
                eps2=eps_co[2][it]-eps[2][it];
                eps3=eps_co[3][it]-eps[3][it];
                eps4=eps_co[4][it]-eps[4][it];
                eps5=eps_co[5][it]-eps[5][it];  
                                
  	        	eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*(eps1+eps2);
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*(eps2+eps0);
                
	        	eps[2][it]-=ic12*(eps0+eps1);
                eps[2][it]-=ic11*eps2;

 	        	eps[3][it]-=ic44*eps3;  
                eps[4][it]-=ic44*eps4;  
                eps[5][it]-=ic44*eps5;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co1=max(err_co1,abs(eps[3][it]-eps_co[3][it]));
                err_co1=max(err_co1,abs(eps[4][it]-eps_co[4][it]));
                err_co1=max(err_co1,abs(eps[5][it]-eps_co[5][it]));
                err_co+=err_co1*err_co1;
			    
		//err_co+=pow(max(max(max(max(max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it])),abs(eps[2][it]-eps_co[2][it])),abs(eps[3][it]-eps_co[3][it])),abs(eps[4][it]-eps_co[4][it])),abs(eps[5][it]-eps_co[5][it])),2);
				
       }
    }
}		

err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]+di[3]*di[3]+di[4]*di[4]+di[5]*di[5]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
//cout<<"temps 8 :"<<(CPUtime()-t8i)<<endl;  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];
sigma[3]=ms[3];
sigma[4]=ms[4];
sigma[5]=ms[5];

fft3(eps_co,N1,N2,N3,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}


void reso_fft3_elas_ctr(bool * LIST_N, R  * LIST_S, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma, int dir, int dir2, R & cmin, R & cmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=6;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R NX,NY,NZ;
//R CC[6][6];
R c11,c12,c44;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Déformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];
			
R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2*N3];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2*N3];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2*N3;
	eps_co[j]=eps_co[0]+j*ND1*N2*N3;
	sig[j]=sig[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(sig,N1,N2,N3,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-3)&&(n_it<2000)){
n_it++;

if(err_co<1e-4){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c44=C2[3][3];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c44=C1[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c11*eps[1][it]+c12*(eps[2][it]+eps[0][it]);
			sig[2][it]=c11*eps[2][it]+c12*(eps[0][it]+eps[1][it]);
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];	
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];
			ms[3]+=sig[3][it];
			ms[4]+=sig[4][it];
			ms[5]+=sig[5][it];		

			err_sig1=max(abs(sig[0][it]),err_sig1);
			err_sig1=max(abs(sig[1][it]),err_sig1);
			err_sig1=max(abs(sig[2][it]),err_sig1);
			err_sig1=max(abs(sig[3][it]),err_sig1);
			err_sig1=max(abs(sig[4][it]),err_sig1);
			err_sig1=max(abs(sig[5][it]),err_sig1);		
			
       }
    }
}		

cout<<"err_sig0 :"<<err_sig0<<", "<<"err_sig1 :"<<err_sig1<<endl;
err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(sig,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi])*(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi]);
		errf+=((NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr])*(NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr]));
		
		errf+=(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi])*(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi]);
		errf+=((NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr])*(NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr]));
		
		errf+=(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi])*(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi]);
		errf+=((NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr])*(NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_sigg+=errf;	
	   }
	}
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]+sig[3][0]*sig[3][0]+sig[4][0]*sig[4][0]+sig[5][0]*sig[5][0]));
cout<<"err_sigg :"<<err_sigg<<endl;

*/

}

//R t4i=CPUtime();

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c44=C2pC0[3][3];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c44=C1pC0[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c12*(eps[0][it]+eps[2][it])+c11*eps[1][it];
			sig[2][it]=c12*(eps[0][it]+eps[1][it])+c11*eps[2][it];
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];		
			
       }
    }
}		

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", sig :"<<sig[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 4 :"<<(CPUtime()-t4i)<<endl;
//R t5i=CPUtime();

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(sig,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1+k*ND1*N2]<<", "<<sig[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 5 :"<<(CPUtime()-t5i)<<endl;
//R t6i=CPUtime();

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

	
int indi,indr;
R coef3=(3.*K0+G0)/(9.*G0*K0);
R coef4=(-3.*K0+2.*G0)/(12.*G0*K0);
R coef5=1./(2.*G0);		
/*		
for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    R NYZ=NY*NZ;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R NXY=NX*NY;
		R NXZ=NX*NZ;
		R norm2=NXX+NYY+NZZ;
		R norm4=norm2*norm2;
		R coef1=1./(4*G0*norm2);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0)*norm4);

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
			//0 -> G11			//1 -> G12			//2 -> G13
			//3 -> G14			//4 -> G15			//5 -> G16
			//6 -> G22			//7 -> G23			//8 -> G24
			//9 -> G25			//10-> G26			//11-> G33
			//12-> G34			//13-> G35			//14-> G36
			//15-> G44			//16-> G45			//17-> G46
			//18-> G55			//19-> G56			//20-> G66
	
			R gg0=(4.*coef1-coef2*NXX)*NXX;
			R gg1=-coef2*NXX*NYY;
			R gg2=-coef2*NXX*NZZ;
			R gg3=-coef2*NXX*NYZ;
			R gg4=(2.*coef1-coef2*NXX)*NXZ;
			R gg5=(2.*coef1-coef2*NXX)*NXY;
			R gg6=(4.*coef1-coef2*NYY)*NYY;
			R gg7=-coef2*NYY*NZZ;
            R gg8=(2.*coef1-coef2*NYY)*NYZ;
			R gg9=-coef2*NYY*NXZ;
            R gg10=(2.*coef1-coef2*NYY)*NXY;
            R gg11=(4.*coef1-coef2*NZZ)*NZZ;
            R gg12=(2.*coef1-coef2*NZZ)*NYZ;
			R gg13=(2.*coef1-coef2*NZZ)*NXZ;
		    R gg14=-coef2*NXY*NZZ;
            R gg15=coef1*(NYY+NZZ)-coef2*NYY*NZZ;
            R gg16=(coef1-coef2*NZZ)*NXY;
            R gg17=(coef1-coef2*NYY)*NXZ;
			R gg18=coef1*(NXX+NZZ)-coef2*NXX*NZZ;
			R gg19=(coef1-coef2*NXX)*NYZ;
			R gg20=coef1*(NXX+NYY)-coef2*NXX*NYY;

              eps_co[0][indr] =gg0*sig[0][indr]+gg1*sig[1][indr]+gg2*sig[2][indr];
			  eps_co[0][indr]+=2*(gg3*sig[3][indr]+gg4*sig[4][indr]+gg5*sig[5][indr]);
			  
			  eps_co[1][indr] =gg1*sig[0][indr]+gg6*sig[1][indr]+gg7*sig[2][indr];
			  eps_co[1][indr]+=2*(gg8*sig[3][indr]+gg9*sig[4][indr]+gg10*sig[5][indr]);
			  
			  eps_co[2][indr] =gg2*sig[0][indr]+gg7*sig[1][indr]+gg11*sig[2][indr];
			  eps_co[2][indr]+=2*(gg12*sig[3][indr]+gg13*sig[4][indr]+gg14*sig[5][indr]);
			  
			  eps_co[3][indr] =gg3*sig[0][indr]+gg8*sig[1][indr]+gg12*sig[2][indr];
			  eps_co[3][indr]+=2*(gg15*sig[3][indr]+gg16*sig[4][indr]+gg17*sig[5][indr]);
			  
			  eps_co[4][indr] =gg4*sig[0][indr]+gg9*sig[1][indr]+gg13*sig[2][indr];
			  eps_co[4][indr]+=2*(gg16*sig[3][indr]+gg18*sig[4][indr]+gg19*sig[5][indr]);
			   
			  eps_co[5][indr] =gg5*sig[0][indr]+gg10*sig[1][indr]+gg14*sig[2][indr];
			  eps_co[5][indr]+=2*(gg17*sig[3][indr]+gg19*sig[4][indr]+gg20*sig[5][indr]);                             
 
			  eps_co[0][indi] =gg0*sig[0][indi]+gg1*sig[1][indi]+gg2*sig[2][indi];
			  eps_co[0][indi]+=2*(gg3*sig[3][indi]+gg4*sig[4][indi]+gg5*sig[5][indi]);
			  
			  eps_co[1][indi] =gg1*sig[0][indi]+gg6*sig[1][indi]+gg7*sig[2][indi];
			  eps_co[1][indi]+=2*(gg8*sig[3][indi]+gg9*sig[4][indi]+gg10*sig[5][indi]);
			  
			  eps_co[2][indi] =gg2*sig[0][indi]+gg7*sig[1][indi]+gg11*sig[2][indi];
			  eps_co[2][indi]+=2*(gg12*sig[3][indi]+gg13*sig[4][indi]+gg14*sig[5][indi]);
			  
			  eps_co[3][indi] =gg3*sig[0][indi]+gg8*sig[1][indi]+gg12*sig[2][indi];
			  eps_co[3][indi]+=2*(gg15*sig[3][indi]+gg16*sig[4][indi]+gg17*sig[5][indi]);
			  
			  eps_co[4][indi] =gg4*sig[0][indi]+gg9*sig[1][indi]+gg13*sig[2][indi];
			  eps_co[4][indi]+=2*(gg16*sig[3][indi]+gg18*sig[4][indi]+gg19*sig[5][indi]);
			   
			  eps_co[5][indi] =gg5*sig[0][indi]+gg10*sig[1][indi]+gg14*sig[2][indi];
			  eps_co[5][indi]+=2*(gg17*sig[3][indi]+gg19*sig[4][indi]+gg20*sig[5][indi]);   
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

}
*/
 

for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R norm2=NXX+NYY+NZZ;
		norm2=1./norm2;
		R norm4=norm2*norm2;
		
		R coef1=1./(G0);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0));

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
	        R divs1 = NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr];
	        R divs2 = NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr];
	        R divs3 = NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr];	
	        
	        R dd = NX*divs1+NY*divs2+NZ*divs3;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        R u3 = coef1*norm2*divs3-coef2*norm4*NZ*dd;	
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = NZ*u3;
	        eps_co[3][indr] = (NY*u3+NZ*u2)/2.;
	        eps_co[4][indr] = (NX*u3+NZ*u1)/2.;
	        eps_co[5][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi];
	        R idivs2 = NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi];
	        R idivs3 = NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2+NZ*idivs3;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        R iu3 = coef1*norm2*idivs3-coef2*norm4*NZ*idd;	
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = NZ*iu3;
	        eps_co[3][indi] = (NY*iu3+NZ*iu2)/2.;
	        eps_co[4][indi] = (NX*iu3+NZ*iu1)/2.;
	        eps_co[5][indi] = (NX*iu2+NY*iu1)/2.;	        	
	        	
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

} 


/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<eps_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 6 :"<<(CPUtime()-t6i)<<endl;
//R t7i=CPUtime();

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(eps_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 7 :"<<(CPUtime()-t7i)<<endl;
//R t8i=CPUtime();

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
//R iCX_C0[6][6];
R ic11,ic12,ic44;
R eps0,eps1,eps2,eps3,eps4,eps5;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic44=iC2_C0[3][3];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic44=iC1_C0[3][3];		
			}
			
                eps0=eps_co[0][it]-eps[0][it];
                eps1=eps_co[1][it]-eps[1][it];
                eps2=eps_co[2][it]-eps[2][it];
                eps3=eps_co[3][it]-eps[3][it];
                eps4=eps_co[4][it]-eps[4][it];
                eps5=eps_co[5][it]-eps[5][it];  
                                
  	        	eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*(eps1+eps2);
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*(eps2+eps0);
                
	        	eps[2][it]-=ic12*(eps0+eps1);
                eps[2][it]-=ic11*eps2;

 	        	eps[3][it]-=ic44*eps3;  
                eps[4][it]-=ic44*eps4;  
                eps[5][it]-=ic44*eps5;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co1=max(err_co1,abs(eps[3][it]-eps_co[3][it]));
                err_co1=max(err_co1,abs(eps[4][it]-eps_co[4][it]));
                err_co1=max(err_co1,abs(eps[5][it]-eps_co[5][it]));
                err_co+=err_co1*err_co1;
			    
		//err_co+=pow(max(max(max(max(max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it])),abs(eps[2][it]-eps_co[2][it])),abs(eps[3][it]-eps_co[3][it])),abs(eps[4][it]-eps_co[4][it])),abs(eps[5][it]-eps_co[5][it])),2);
				
       }
    }
}		

err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]+di[3]*di[3]+di[4]*di[4]+di[5]*di[5]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
//cout<<"temps 8 :"<<(CPUtime()-t8i)<<endl;  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];
sigma[3]=ms[3];
sigma[4]=ms[4];
sigma[5]=ms[5];

cmin =2e11;
cmax =-2e11;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c44=C2[3][3];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c44=C1[3][3];		
			}
				
           if(dir2==1){LIST_S[itc]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);}
           else if(dir2==2){LIST_S[itc]=c11*eps[1][it]+c12*(eps[2][it]+eps[0][it]);}
           else if(dir2==3){LIST_S[itc]=c11*eps[2][it]+c12*(eps[0][it]+eps[1][it]);}
           else if(dir2==4){LIST_S[itc]=c44*eps[3][it];}
           else if(dir2==5){LIST_S[itc]=c44*eps[4][it];} 
           else if(dir2==6){LIST_S[itc]=c44*eps[5][it];} 

           cmin=(cmin<LIST_S[itc])?cmin:LIST_S[itc];
           cmax=(cmax>LIST_S[itc])?cmax:LIST_S[itc];
			
       }
    }
}

fft3(eps_co,N1,N2,N3,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}


void reso_fft3_elas_dpct(bool * LIST_N, R  * LIST_U, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma, int dir, R & cmin, R & cmax)
{
cout<<"Resolution fft - schema Eyre and Milton - dir : "<<dir<<endl;

bool boolin=1;
int NDIM=6;
int ND1,IP1;
int i,j,k,itc,it,i1,i2,i3;
R NX,NY,NZ;
//R CC[6][6];
R c11,c12,c44;

if(N1%2==0) {IP1=2;} else {IP1=1;}
ND1=N1+IP1;

//Déformation imposée
R di[NDIM];
for(int j=0;j<NDIM;j++){
di[j]=0.;	
}
di[dir-1]=1.;

if((dir<1)||(dir>NDIM)){
boolin=0;  
cout<<"Mauvaise direction"<<endl;  
}

if(boolin==1){

//valeurs moyennes à évaluer
R ms[NDIM];
			
R ** eps  = new R * [NDIM];
eps[0] = new R [NDIM*ND1*N2*N3];
R ** eps_co  = new R * [NDIM];
eps_co[0] = new R [NDIM*ND1*N2*N3];
R ** sig  = new R * [NDIM];
sig[0] = new R [NDIM*ND1*N2*N3];

for(j=0;j<NDIM;j++){
	if(j>0){
	eps[j]=eps[0]+j*ND1*N2*N3;
	eps_co[j]=eps_co[0]+j*ND1*N2*N3;
	sig[j]=sig[0]+j*ND1*N2*N3;
	}		
	for(i=0;i<ND1*N2*N3;i++){
		if((i%ND1)<N1){
		eps[j][i]=di[j];
	    }
	    else
	    {
		eps[j][i]=0.;	
     	}
		eps_co[j][i]=0.;
		sig[j][i]=0.;
	}
}

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
	 for(k=0;k<N3;k++){		
	  for(j=0;j<N2;j++){	
		for(i=0;i<N1;i++){
		cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
		}
	  }
	 } 
}
}*/

// Initialization
fft3(sig,N1,N2,N3,NDIM,0); // ?	

R err_eq=1.;
R err_co=1.;
R err_sig0=1.;
R err_sig1;
int n_it=-1;

while((max(err_eq,err_co)>1e-3)&&(n_it<2000)){
n_it++;

if(err_co<1e-4){
//if (DEBUG) cout<<"step1"<<endl;
//sig=C*eps quelque soit le pixel 

err_sig1=0.;
for(j=0;j<NDIM;j++){
ms[j]=0.;	
}

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c44=C2[3][3];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c44=C1[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c11*eps[1][it]+c12*(eps[2][it]+eps[0][it]);
			sig[2][it]=c11*eps[2][it]+c12*(eps[0][it]+eps[1][it]);
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];	
					
			ms[0]+=sig[0][it];
			ms[1]+=sig[1][it];
			ms[2]+=sig[2][it];
			ms[3]+=sig[3][it];
			ms[4]+=sig[4][it];
			ms[5]+=sig[5][it];		

			err_sig1=max(abs(sig[0][it]),err_sig1);
			err_sig1=max(abs(sig[1][it]),err_sig1);
			err_sig1=max(abs(sig[2][it]),err_sig1);
			err_sig1=max(abs(sig[3][it]),err_sig1);
			err_sig1=max(abs(sig[4][it]),err_sig1);
			err_sig1=max(abs(sig[5][it]),err_sig1);		
			
       }
    }
}		

cout<<"err_sig0 :"<<err_sig0<<", "<<"err_sig1 :"<<err_sig1<<endl;
err_eq=abs(err_sig1-err_sig0)/err_sig0;
err_sig0=err_sig1;

/*
//if(DEBUG) cout<<"step2"<<endl;
// Direct fft
fft3(sig,N1,N2,N3,NDIM,-1);

//if(DEBUG) cout<<"step3"<<endl;
//calcul de l'erreur/contraintes/deformation

R err_sigg=0.;
int ic,ic1,ic2,ic3;
for(i3=0;i3<N3;i3++){
	ic3=0;
	if((i3==N3/2)&&(N3%2==0)) {ic3=1;}		
	
		NZ=R(i3);
		if(NZ>N3/2){NZ=NZ-N3;}	
		
	for(i2=0;i2<N2;i2++){
		ic2=0;
		if((i2==N2/2)&&(N2%2==0)) {ic2=1;}	
		
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
			
	   for(i1=0;i1<(N1/2)+1;i1++){
		ic1=0;
		if((i1==N1/2)&&(N1%2==0)) ic1=1;	
		ic=ic1+ic2+ic3;

		NX=R(i1);
		int indr=2*i1+i2*ND1+i3*ND1*N2;
		int indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R errf=(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi])*(NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi]);
		errf+=((NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr])*(NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr]));
		
		errf+=(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi])*(NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi]);
		errf+=((NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr])*(NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr]));
		
		errf+=(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi])*(NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi]);
		errf+=((NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr])*(NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr]));
		
		if(ic==2){errf=errf*0.125;}
		else if(ic==2){errf=errf*0.25;}
		else if(ic==1){errf=errf*0.5;}	
		if((i1!=0)&&((i1!=N1/2)||(N1%2!=0))){errf=2.*errf;}
		
	   err_sigg+=errf;	
	   }
	}
}
err_sigg=sqrt(err_sigg/(sig[0][0]*sig[0][0]+sig[1][0]*sig[1][0]+sig[2][0]*sig[2][0]+sig[3][0]*sig[3][0]+sig[4][0]*sig[4][0]+sig[5][0]*sig[5][0]));
cout<<"err_sigg :"<<err_sigg<<endl;

*/

}

//R t4i=CPUtime();

//if(DEBUG) cout<<"step4"<<endl;
//calcul de tau
//tau=(C+C0)*eps quelquesoit le pixel 

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2pC0[0][0]; 
			c12=C2pC0[0][1]; 
			c44=C2pC0[3][3];
			}else{
			c11=C1pC0[0][0]; 
			c12=C1pC0[0][1]; 
			c44=C1pC0[3][3];		
			}
			
			sig[0][it]=c11*eps[0][it]+c12*(eps[1][it]+eps[2][it]);
			sig[1][it]=c12*(eps[0][it]+eps[2][it])+c11*eps[1][it];
			sig[2][it]=c12*(eps[0][it]+eps[1][it])+c11*eps[2][it];
			sig[3][it]=c44*eps[3][it];
			sig[4][it]=c44*eps[4][it];
			sig[5][it]=c44*eps[5][it];		
			
       }
    }
}		

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", sig :"<<sig[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 4 :"<<(CPUtime()-t4i)<<endl;
//R t5i=CPUtime();

//if(DEBUG) cout<<"step5"<<endl;
// Direct fft - tau
fft3(sig,N1,N2,N3,NDIM,-1);

/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<(N1/2)<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", sig :"<<sig[jt][2*i+j*ND1+k*ND1*N2]<<", "<<sig[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 5 :"<<(CPUtime()-t5i)<<endl;
//R t6i=CPUtime();

//if(DEBUG) cout<<"step6"<<endl;
//eps_co

	
int indi,indr;
R coef3=(3.*K0+G0)/(9.*G0*K0);
R coef4=(-3.*K0+2.*G0)/(12.*G0*K0);
R coef5=1./(2.*G0);		
/*		
for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    R NYZ=NY*NZ;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){
			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R NXY=NX*NY;
		R NXZ=NX*NZ;
		R norm2=NXX+NYY+NZZ;
		R norm4=norm2*norm2;
		R coef1=1./(4*G0*norm2);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0)*norm4);

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
			//0 -> G11			//1 -> G12			//2 -> G13
			//3 -> G14			//4 -> G15			//5 -> G16
			//6 -> G22			//7 -> G23			//8 -> G24
			//9 -> G25			//10-> G26			//11-> G33
			//12-> G34			//13-> G35			//14-> G36
			//15-> G44			//16-> G45			//17-> G46
			//18-> G55			//19-> G56			//20-> G66
	
			R gg0=(4.*coef1-coef2*NXX)*NXX;
			R gg1=-coef2*NXX*NYY;
			R gg2=-coef2*NXX*NZZ;
			R gg3=-coef2*NXX*NYZ;
			R gg4=(2.*coef1-coef2*NXX)*NXZ;
			R gg5=(2.*coef1-coef2*NXX)*NXY;
			R gg6=(4.*coef1-coef2*NYY)*NYY;
			R gg7=-coef2*NYY*NZZ;
            R gg8=(2.*coef1-coef2*NYY)*NYZ;
			R gg9=-coef2*NYY*NXZ;
            R gg10=(2.*coef1-coef2*NYY)*NXY;
            R gg11=(4.*coef1-coef2*NZZ)*NZZ;
            R gg12=(2.*coef1-coef2*NZZ)*NYZ;
			R gg13=(2.*coef1-coef2*NZZ)*NXZ;
		    R gg14=-coef2*NXY*NZZ;
            R gg15=coef1*(NYY+NZZ)-coef2*NYY*NZZ;
            R gg16=(coef1-coef2*NZZ)*NXY;
            R gg17=(coef1-coef2*NYY)*NXZ;
			R gg18=coef1*(NXX+NZZ)-coef2*NXX*NZZ;
			R gg19=(coef1-coef2*NXX)*NYZ;
			R gg20=coef1*(NXX+NYY)-coef2*NXX*NYY;

              eps_co[0][indr] =gg0*sig[0][indr]+gg1*sig[1][indr]+gg2*sig[2][indr];
			  eps_co[0][indr]+=2*(gg3*sig[3][indr]+gg4*sig[4][indr]+gg5*sig[5][indr]);
			  
			  eps_co[1][indr] =gg1*sig[0][indr]+gg6*sig[1][indr]+gg7*sig[2][indr];
			  eps_co[1][indr]+=2*(gg8*sig[3][indr]+gg9*sig[4][indr]+gg10*sig[5][indr]);
			  
			  eps_co[2][indr] =gg2*sig[0][indr]+gg7*sig[1][indr]+gg11*sig[2][indr];
			  eps_co[2][indr]+=2*(gg12*sig[3][indr]+gg13*sig[4][indr]+gg14*sig[5][indr]);
			  
			  eps_co[3][indr] =gg3*sig[0][indr]+gg8*sig[1][indr]+gg12*sig[2][indr];
			  eps_co[3][indr]+=2*(gg15*sig[3][indr]+gg16*sig[4][indr]+gg17*sig[5][indr]);
			  
			  eps_co[4][indr] =gg4*sig[0][indr]+gg9*sig[1][indr]+gg13*sig[2][indr];
			  eps_co[4][indr]+=2*(gg16*sig[3][indr]+gg18*sig[4][indr]+gg19*sig[5][indr]);
			   
			  eps_co[5][indr] =gg5*sig[0][indr]+gg10*sig[1][indr]+gg14*sig[2][indr];
			  eps_co[5][indr]+=2*(gg17*sig[3][indr]+gg19*sig[4][indr]+gg20*sig[5][indr]);                             
 
			  eps_co[0][indi] =gg0*sig[0][indi]+gg1*sig[1][indi]+gg2*sig[2][indi];
			  eps_co[0][indi]+=2*(gg3*sig[3][indi]+gg4*sig[4][indi]+gg5*sig[5][indi]);
			  
			  eps_co[1][indi] =gg1*sig[0][indi]+gg6*sig[1][indi]+gg7*sig[2][indi];
			  eps_co[1][indi]+=2*(gg8*sig[3][indi]+gg9*sig[4][indi]+gg10*sig[5][indi]);
			  
			  eps_co[2][indi] =gg2*sig[0][indi]+gg7*sig[1][indi]+gg11*sig[2][indi];
			  eps_co[2][indi]+=2*(gg12*sig[3][indi]+gg13*sig[4][indi]+gg14*sig[5][indi]);
			  
			  eps_co[3][indi] =gg3*sig[0][indi]+gg8*sig[1][indi]+gg12*sig[2][indi];
			  eps_co[3][indi]+=2*(gg15*sig[3][indi]+gg16*sig[4][indi]+gg17*sig[5][indi]);
			  
			  eps_co[4][indi] =gg4*sig[0][indi]+gg9*sig[1][indi]+gg13*sig[2][indi];
			  eps_co[4][indi]+=2*(gg16*sig[3][indi]+gg18*sig[4][indi]+gg19*sig[5][indi]);
			   
			  eps_co[5][indi] =gg5*sig[0][indi]+gg10*sig[1][indi]+gg14*sig[2][indi];
			  eps_co[5][indi]+=2*(gg17*sig[3][indi]+gg19*sig[4][indi]+gg20*sig[5][indi]);   
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

}
*/
 

for(i3=0;i3<N3;i3++){
	
    NZ=R(i3);
	if(NZ>N3/2){NZ=NZ-N3;}	
	R NZZ=NZ*NZ;
	
	for(i2=0;i2<N2;i2++){
		NY=R(i2);
		if(NY>N2/2){NY=NY-N2;}
	    R NYY=NY*NY;
	    
	   for(i1=0;i1<(N1/2)+1;i1++){			   
		NX=R(i1);
		
		indr=2*i1+i2*ND1+i3*ND1*N2;
		indi=2*i1+1+i2*ND1+i3*ND1*N2;
		
		R NXX=NX*NX;
		R norm2=NXX+NYY+NZZ;
		norm2=1./norm2;
		R norm4=norm2*norm2;
		
		R coef1=1./(G0);
		R coef2=(3*K0+G0)/(G0*(3*K0+4*G0));

		
			if((i1!=N1/2)&&(i2!=N2/2)&&(i3!=N3/2)&&((i1!=0)||(i2!=0)||(i3!=0))){  
	
	        R divs1 = NX*sig[0][indr]+NY*sig[5][indr]+NZ*sig[4][indr];
	        R divs2 = NX*sig[5][indr]+NY*sig[1][indr]+NZ*sig[3][indr];
	        R divs3 = NX*sig[4][indr]+NY*sig[3][indr]+NZ*sig[2][indr];	
	        
	        R dd = NX*divs1+NY*divs2+NZ*divs3;
	        
	        R u1 = coef1*norm2*divs1-coef2*norm4*NX*dd;
	        R u2 = coef1*norm2*divs2-coef2*norm4*NY*dd;
	        R u3 = coef1*norm2*divs3-coef2*norm4*NZ*dd;	
	        
	        eps_co[0][indr] = NX*u1;
	        eps_co[1][indr] = NY*u2;
	        eps_co[2][indr] = NZ*u3;
	        eps_co[3][indr] = (NY*u3+NZ*u2)/2.;
	        eps_co[4][indr] = (NX*u3+NZ*u1)/2.;
	        eps_co[5][indr] = (NX*u2+NY*u1)/2.;
	        	
	        R idivs1 = NX*sig[0][indi]+NY*sig[5][indi]+NZ*sig[4][indi];
	        R idivs2 = NX*sig[5][indi]+NY*sig[1][indi]+NZ*sig[3][indi];
	        R idivs3 = NX*sig[4][indi]+NY*sig[3][indi]+NZ*sig[2][indi];	
	        
	        R idd = NX*idivs1+NY*idivs2+NZ*idivs3;
	        
	        R iu1 = coef1*norm2*idivs1-coef2*norm4*NX*idd;
	        R iu2 = coef1*norm2*idivs2-coef2*norm4*NY*idd;
	        R iu3 = coef1*norm2*idivs3-coef2*norm4*NZ*idd;	
	        
	        eps_co[0][indi] = NX*iu1;
	        eps_co[1][indi] = NY*iu2;
	        eps_co[2][indi] = NZ*iu3;
	        eps_co[3][indi] = (NY*iu3+NZ*iu2)/2.;
	        eps_co[4][indi] = (NX*iu3+NZ*iu1)/2.;
	        eps_co[5][indi] = (NX*iu2+NY*iu1)/2.;	        	
	        	
			  
			}
			else if((i1!=0)||(i2!=0)||(i3!=0))
			{
				
			  eps_co[0][indr]=(coef3*sig[0][indr]+coef4*(sig[1][indr]+sig[2][indr]));  
			  eps_co[1][indr]=(coef4*(sig[0][indr]+sig[2][indr])+coef3*sig[1][indr]);
			  eps_co[2][indr]=(coef4*(sig[0][indr]+sig[1][indr])+coef3*sig[2][indr]);
			  eps_co[3][indr]=coef5*sig[3][indr];
			  eps_co[4][indr]=coef5*sig[4][indr];  
			  eps_co[5][indr]=coef5*sig[5][indr];   				

			  eps_co[0][indi]=(coef3*sig[0][indi]+coef4*(sig[1][indi]+sig[2][indi]));  
			  eps_co[1][indi]=(coef4*(sig[0][indi]+sig[2][indi])+coef3*sig[1][indi]);
			  eps_co[2][indi]=(coef4*(sig[0][indi]+sig[1][indi])+coef3*sig[2][indi]);
			  eps_co[3][indi]=coef5*sig[3][indi];
			  eps_co[4][indi]=coef5*sig[4][indi];  
			  eps_co[5][indi]=coef5*sig[5][indi];   	
		
			}
			else
			{
			eps_co[0][indr]=di[0];
			eps_co[0][indi]=0.;
			eps_co[1][indr]=di[1];
			eps_co[1][indi]=0.	;
			eps_co[2][indr]=di[2];
			eps_co[2][indi]=0.	;		
			eps_co[3][indr]=di[3];
			eps_co[3][indi]=0.	;
			eps_co[4][indr]=di[4];
			eps_co[4][indi]=0.	;
			eps_co[5][indr]=di[5];
			eps_co[5][indi]=0.	;											
				
			} 	

	   }
	}

} 


/*if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<(N1/2);i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][2*i+j*ND1+k*ND1*N2]<<", "<<eps_co[jt][2*i+1+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 6 :"<<(CPUtime()-t6i)<<endl;
//R t7i=CPUtime();

//if(DEBUG) cout<<"step7"<<endl;
// FFT reverse 
fft3(eps_co,N1,N2,N3,NDIM,1); // re initialization ???

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/

//cout<<"temps 7 :"<<(CPUtime()-t7i)<<endl;
//R t8i=CPUtime();

//if(DEBUG) cout<<"step8"<<endl;
//calcul de l'erreur comp//itération i+1
err_co=0.;
//R iCX_C0[6][6];
R ic11,ic12,ic44;
R eps0,eps1,eps2,eps3,eps4,eps5;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			ic11=iC2_C0[0][0]; 
			ic12=iC2_C0[0][1]; 
			ic44=iC2_C0[3][3];
			}else{
			ic11=iC1_C0[0][0]; 
			ic12=iC1_C0[0][1]; 
			ic44=iC1_C0[3][3];		
			}
			
                eps0=eps_co[0][it]-eps[0][it];
                eps1=eps_co[1][it]-eps[1][it];
                eps2=eps_co[2][it]-eps[2][it];
                eps3=eps_co[3][it]-eps[3][it];
                eps4=eps_co[4][it]-eps[4][it];
                eps5=eps_co[5][it]-eps[5][it];  
                                
  	        	eps[0][it]-=ic11*eps0;
                eps[0][it]-=ic12*(eps1+eps2);
                
                eps[1][it]-=ic11*eps1;
                eps[1][it]-=ic12*(eps2+eps0);
                
	        	eps[2][it]-=ic12*(eps0+eps1);
                eps[2][it]-=ic11*eps2;

 	        	eps[3][it]-=ic44*eps3;  
                eps[4][it]-=ic44*eps4;  
                eps[5][it]-=ic44*eps5;  
                
                R err_co1=max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it]));
                err_co1=max(err_co1,abs(eps[2][it]-eps_co[2][it]));
                err_co1=max(err_co1,abs(eps[3][it]-eps_co[3][it]));
                err_co1=max(err_co1,abs(eps[4][it]-eps_co[4][it]));
                err_co1=max(err_co1,abs(eps[5][it]-eps_co[5][it]));
                err_co+=err_co1*err_co1;
			    
		//err_co+=pow(max(max(max(max(max(abs(eps[0][it]-eps_co[0][it]),abs(eps[1][it]-eps_co[1][it])),abs(eps[2][it]-eps_co[2][it])),abs(eps[3][it]-eps_co[3][it])),abs(eps[4][it]-eps_co[4][it])),abs(eps[5][it]-eps_co[5][it])),2);
				
       }
    }
}		

err_co/=N_TOT;
err_co=err_co/sqrt(di[0]*di[0]+di[1]*di[1]+di[2]*di[2]+di[3]*di[3]+di[4]*di[4]+di[5]*di[5]);

/*
if(DEBUG){
for(int jt=0;jt<NDIM;jt++){
 for(k=0;k<N3;k++){		
  for(j=0;j<N2;j++){	
	for(i=0;i<N1;i++){
    cout<<"jt :"<<jt<<", i :"<<i<<", j :"<<j<<", k :"<<k<<", eps :"<<eps[jt][i+j*ND1+k*ND1*N2]<<", eps_co :"<<eps_co[jt][i+j*ND1+k*ND1*N2]<<endl;
    }
  }
 } 
}
}*/
  
//cout<<"temps 8 :"<<(CPUtime()-t8i)<<endl;  
cout<<"n_it : "<<n_it<<", err_eq : "<<err_eq<<", err_co :"<<err_co<<endl;


}

sigma[0]=ms[0];
sigma[1]=ms[1];
sigma[2]=ms[2];
sigma[3]=ms[3];
sigma[4]=ms[4];
sigma[5]=ms[5];

cmin =2e11;
cmax =-2e11;

for(i3=0;i3<N3;i3++){
	for(i2=0;i2<N2;i2++){
	   for(i1=0;i1<N1;i1++){
		itc=i1+i2*N1+i3*N1*N2;  
        it=i1+i2*ND1+i3*ND1*N2;
        
			if(LIST_N[itc]){
			c11=C2[0][0]; 
			c12=C2[0][1]; 
			c44=C2[3][3];
			}else{
			c11=C1[0][0]; 
			c12=C1[0][1]; 
			c44=C1[3][3];		
			}
	
            if(dir==1){LIST_U[itc]=eps[0][it]*(i1+0.5);}
            else if(dir==2){LIST_U[itc]=eps[1][it]*(i2+0.5);}
            else if(dir==3){LIST_U[itc]=eps[2][it]*(i3+0.5);}     
    
           cmin=(cmin<LIST_U[itc])?cmin:LIST_U[itc];
           cmax=(cmax>LIST_U[itc])?cmax:LIST_U[itc];
		
       }
    }
}

fft3(eps_co,N1,N2,N3,NDIM,9);
delete [] eps[0];
delete [] eps_co[0];
delete [] sig[0];
delete [] eps;
delete [] eps_co;
delete [] sig;
}

return;
}











