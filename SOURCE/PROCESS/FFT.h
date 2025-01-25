#ifndef ___FFT___
#define ___FFT___

typedef double R;

void fft2(double **buf, int N1,int N2,int NDIM,int isign);
void fft3(R **buf, int N1,int N2,int N3,int NDIM,int isign);

#endif
