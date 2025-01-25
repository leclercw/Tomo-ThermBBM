#ifndef __RESO_FFT__
#define __RESO_FFT__

typedef double R;

void reso_fft2_ther(bool * LIST_N, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir);
void reso_fft2_ther_temp(bool * LIST_N, R * LIST_T, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu, int dir, R & tmin, R & tmax);
void reso_fft2_ther_flux(bool * LIST_N, R * LIST_S, int H1, int H2, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu, int dir, int dir2, R & fmin, R & fmax);

void reso_fft2_elas(bool * LIST_N, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma,int dir);
void reso_fft2_elas_dpct(bool * LIST_N, R * LIST_U, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma, int dir, R & cmin, R & cmax);
void reso_fft2_elas_ctr(bool * LIST_N, R * LIST_S, int H1, int H2, int N_TOT,  R G0, R K0, R C0[][3], R C1[][3], R C2[][3], R C1pC0[][3], R C2pC0[][3], R iC1_C0[][3], R iC2_C0[][3], R * sigma, int dir, int dir2, R & cmin, R & cmax);

void reso_fft3_ther(bool * LIST_N, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir);
void reso_fft3_ther_temp(bool * LIST_N, R * LIST_T, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu, int dir, R & tmin, R & tmax);
void reso_fft3_ther_flux(bool * LIST_N, R * LIST_S, int N1, int N2, int N3, int N_TOT, R K0, R K1, R K2, R K1pK0, R K2pK0, R iK1, R iK2, R * flu,int dir, int dir2, R & fmin, R & fmax);

void reso_fft3_elas(bool * LIST_N, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma,int dir);
void reso_fft3_elas_ctr(bool * LIST_N, R * LIST_S, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma, int dir, int dir2, R & cmin, R & cmax);
void reso_fft3_elas_dpct(bool * LIST_N, R  * LIST_U, int N1, int N2, int N3, int N_TOT, R G0, R K0, R C1[][6], R C2[][6], R C1pC0[][6], R C2pC0[][6], R iC1_C0[][6], R iC2_C0[][6],R * sigma, int dir, R & cmin, R & cmax);

#endif
