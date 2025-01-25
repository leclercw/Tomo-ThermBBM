#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "FFT.h"
#include <fftw3.h>

void fft2(double **buf, int N1,int N2,int NDIM,int isign)
{
// isign = -1 : direct FFT
// isign =  0 : initialization
// isign = +1 : reverse FFT	
 
  int i;	
  int i0, i1;	
  int ind;

  int ND_1;
  if((N1%2)==0){ND_1=N1+2;}else{ND_1=N1+1;}
    
  int ND_2=N2;
  
  int ND2_1=ND_1/2;
  int ND2_2=ND_2;   	  
	
  static fftw_plan pland[2];
  static fftw_plan plani[2];
  static double coef;
  
     if (isign==-1) {
      for(ind=0; ind<NDIM; ind++) {		 	  
		for(i1=0; i1<N2; i1++) {
		  fftw_execute_dft_r2c(
					   pland[0],
					   (buf[ind] + i1*ND_1), 
					   (fftw_complex *)(buf[ind] + i1*ND_1)
					   );		  
		}		
		
		for(i0=0; i0<ND2_1; i0++) {		  
		  fftw_execute_dft(
				   pland[1],
				   (fftw_complex *)(buf[ind] + 2*i0), 
				   (fftw_complex *)(buf[ind] + 2*i0)
				   );
		}	
	   }
	   	 
     /* for(i=0; i< ND_1*ND_2 ; i++) {
	   buf[i] *= coef;
      }	*/	
	  	
	 } 
	 
	 else if (isign==1) { 
      for(ind=0; ind<NDIM; ind++) {		 
		for(i0=0; i0<ND2_1; i0++) {
		  fftw_execute_dft(
				   plani[1],
				   (fftw_complex *)(buf[ind] + 2*i0), 
				   (fftw_complex *)(buf[ind] + 2*i0)
				   );
		}		
		for(i1=0; i1<N2; i1++) {		  
		  fftw_execute_dft_c2r(
					   plani[0],
					   (fftw_complex *)(buf[ind] + i1*ND_1), 
									   (buf[ind] + i1*ND_1)
					   );
		}	
      }
      
      for(ind=0; ind<NDIM; ind++) {
		  for(i=0; i< ND_1*ND_2 ; i++) {
		   buf[ind][i] *= coef;
		  }		
      } 		
      		
     }
    
	  else if (isign==0){

      coef = 1. / (N1*N2); 
	   
	  unsigned flag = FFTW_MEASURE;

 
		/* plans for direct transforms */
		pland[0] = fftw_plan_many_dft_r2c(
						  1,                  /* rank= 1 */
						  &N1,                 /* size of dimension */
						  1,                  /* howmany=1 */
						  
						  buf[0],                /* input array */
						  &ND_1,                 /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND_1,                 /* idist */
						  
						  (fftw_complex *)buf[0],                /* output array */
						  &ND2_1,                /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND2_1,                /* idist=1 */
						  flag
						  );
		
		pland[1] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N2,                 /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],                /* input array */
					  &ND2_2,                 /* size of dimension as it has been allocated */
					  ND2_1,                 /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],                /* output array */
					  &ND2_2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  FFTW_FORWARD,
					  flag
					  );  

	   plani[0] = fftw_plan_many_dft_c2r(
						  1,                  /* rank= 1 */
						  &N1,                 /* size of dimension */
						  1,                  /* howmany=1 */
						  
						  (fftw_complex *)buf[0],/* output array */
						  &ND2_1,                /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND2_1,                /* idist=1 */

						  buf[0],                /* input array */
						  &ND_1,                 /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND_1,                 /* idist */

						  flag
						  );

		
		plani[1] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N2,                 /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],/* input array */
					  &ND2_2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],/* output array */
					  &ND2_2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  FFTW_BACKWARD,
					  flag
					  );  
		}	    
		else if (isign==9){  
    fftw_destroy_plan(pland[0]);
    fftw_destroy_plan(pland[1]);
    fftw_destroy_plan(plani[0]);
    fftw_destroy_plan(plani[1]);  		
	}
        
}

void fft3(R **buf, int N1,int N2,int N3,int NDIM,int isign)
{
// isign = -1 : direct FFT
// isign =  0 : initialization
// isign = +1 : reverse FFT	
 
  int i;	
  int i0, i1, i2;	
  int ind;

  int ND_1;
  if((N1%2)==0){ND_1=N1+2;}else{ND_1=N1+1;}  
  int ND2_1=ND_1/2;
	
  static fftw_plan pland[3];
  static fftw_plan plani[3];
  static R coef;
  
     if (isign==-1) {
      for(ind=0; ind<NDIM; ind++) {
		  
		for(i2=0; i2<N3; i2++) {  		 	  
			for(i1=0; i1<N2; i1++) {
			  fftw_execute_dft_r2c(
						   pland[0],
						   (buf[ind] + i1*ND_1 + i2*ND_1*N2), 
						   (fftw_complex *)(buf[ind] + i1*ND_1 + i2*ND_1*N2)
						   );		  
			}	
	    }	
		
		for(i2=0; i2<N3; i2++) {  	
			for(i0=0; i0<ND2_1; i0++) {		  
			  fftw_execute_dft(
					   pland[1],
					   (fftw_complex *)(buf[ind] + 2*i0 + i2*ND_1*N2), 
					   (fftw_complex *)(buf[ind] + 2*i0 + i2*ND_1*N2)
					   );
			}	
	    }
	    
		for(i1=0; i1<N2; i1++) {  	
			for(i0=0; i0<ND2_1; i0++) {		  
			  fftw_execute_dft(
					   pland[2],
					   (fftw_complex *)(buf[ind] + 2*i0 + i1*ND_1), 
					   (fftw_complex *)(buf[ind] + 2*i0 + i1*ND_1)
					   );
			}	
	    }	    
	    
		
	   }
	   	 
     /* for(i=0; i< ND_1*N2 ; i++) {
	   buf[i] *= coef;
      }	*/	
	  	
	 } 
	 
	 else if (isign==1) { 
      for(ind=0; ind<NDIM; ind++) {	
		  
	   for(i1=0; i1<N2; i1++) {    
		for(i0=0; i0<ND2_1; i0++) {
		  fftw_execute_dft(
				   plani[2],
				   (fftw_complex *)(buf[ind] + 2*i0 + i1*ND_1), 
				   (fftw_complex *)(buf[ind] + 2*i0 + i1*ND_1)
				   );
		}
	   }
	   for(i2=0; i2<N3; i2++) {  				 
		for(i0=0; i0<ND2_1; i0++) {
		  fftw_execute_dft(
				   plani[1],
				   (fftw_complex *)(buf[ind] + 2*i0 + i2*ND_1*N2), 
				   (fftw_complex *)(buf[ind] + 2*i0 + i2*ND_1*N2)
				   );
		}		
	   }
	   for(i2=0; i2<N3; i2++) {
		for(i1=0; i1<N2; i1++) {		  
		  fftw_execute_dft_c2r(
					   plani[0],
					   (fftw_complex *)(buf[ind] + i1*ND_1 + i2*ND_1*N2), 
									   (buf[ind] + i1*ND_1 + i2*ND_1*N2)
					   );
		}	
	   }
		
      }
      
      for(ind=0; ind<NDIM; ind++) {
		  for(i=0; i< ND_1*N2*N3 ; i++) {
		   buf[ind][i] *= coef;
		  }		
      } 		
      		
     }
    
	  else if (isign==0){

      coef = 1. / (N1*N2*N3); 
	   
	  unsigned flag = FFTW_MEASURE;

 
		/* plans for direct transforms */
		pland[0] = fftw_plan_many_dft_r2c(
						  1,                  /* rank= 1 */
						  &N1,                 /* size of dimension */
						  1,                  /* howmany=1 */
						  
						  buf[0],                /* input array */
						  &ND_1,                 /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND_1,                 /* idist */
						  
						  (fftw_complex *)buf[0],                /* output array */
						  &ND2_1,                /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND2_1,                /* idist=1 */
						  flag
						  );
		
		pland[1] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N2,                 /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],                /* input array */
					  &N2,                 /* size of dimension as it has been allocated */
					  ND2_1,                 /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],                /* output array */
					  &N2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  FFTW_FORWARD,
					  flag
					  );  

		pland[2] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N3,             /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],   /* input array */
					  &N3,            /* size of dimension as it has been allocated */
					  ND2_1*N2,      /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],   /* output array */
					  &N3,            /* size of dimension as it has been allocated */
					  ND2_1*N2,      /* istride */
					  1,                 /* idist */
					  
					  FFTW_FORWARD,
					  flag
					  
					  );


	   plani[0] = fftw_plan_many_dft_c2r(
						  1,                  /* rank= 1 */
						  &N1,                 /* size of dimension */
						  1,                  /* howmany=1 */
						  
						  (fftw_complex *)buf[0],/* output array */
						  &ND2_1,                /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND2_1,                /* idist=1 */

						  buf[0],                /* input array */
						  &ND_1,                 /* size of dimension as it has been allocated */
						  1,                  /* istride = 1 */
						  ND_1,                 /* idist */

						  flag
						  );

		
		plani[1] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N2,                 /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],/* input array */
					  &N2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],/* output array */
					  &N2,                /* size of dimension as it has been allocated */
					  ND2_1,                /* istride */
					  1,                  /* idist */
					  
					  FFTW_BACKWARD,
					  flag
					  );  
					  
		plani[2] = fftw_plan_many_dft(
					  1,                  /* rank= 1 */
					  &N3,             /* size of dimension */
					  1,                  /* howmany=1 */
					  
					  (fftw_complex *)buf[0],   /* input array */
					  &N3,            /* size of dimension as it has been allocated */
					  ND2_1*N2,      /* istride */
					  1,                  /* idist */
					  
					  (fftw_complex *)buf[0],   /* output array */
					  &N3,            /* size of dimension as it has been allocated */
					  ND2_1*N2,      /* istride */
					  1,                 /* idist */
					  
					  FFTW_BACKWARD,
					  flag
					  
					  );					  
					  
					  
		}	    
		else if (isign==9){  
    fftw_destroy_plan(pland[0]);
    fftw_destroy_plan(pland[1]);
    fftw_destroy_plan(pland[2]);    
    fftw_destroy_plan(plani[0]);
    fftw_destroy_plan(plani[1]);  	
    fftw_destroy_plan(plani[2]);  	    	
	}
        
}
