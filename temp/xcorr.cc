// xcorr.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 5月 12 20:46:31 2017 (+0800)
// Last-Updated: 五 5月 12 22:35:30 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 47
// URL: http://wuhongyi.cn 

#include "xcorr.hh"

#include <fftw3.h>
#include <complex>
// #include <complex.h>
#include <iostream>
#include <cstring>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void xcorr_timedomain(void *signala, void *signalb, void *result, int N)
// {
//   for (int tau = 0; tau < 2*N-1; ++tau) {
//     complex acf = 0 + 0*I;
//     for (int i = 0; i < N; ++i) {
//       const int signala_idx = (i+tau)%(2*N-1);
//       const complex conjb = conj(((complex*)signalb)[i]);
//       const double factor = (signala_idx >= N) ?
// 	((complex*)signala)[signala_idx-N] : 1.0;
//       acf += factor * conjb;
//     }
//     ((complex*)result)[tau] = acf;
//   }
//   return;
// }

void xcorr_timedomain(std::complex<double> *signala, std::complex<double> *signalb, std::complex<double> *result, int N)
{
  for (int tau = 0; tau < (2*N-1); ++tau)
    {
      std::complex<double> acf(0,0);
      for (int i = 0; i < N; ++i)
      	{
      	  const int signala_idx = (i+tau)%(2*N-1);
      	  const std::complex<double> conjb = std::conj(((std::complex<double>*)signalb)[i]);
	  // std::cout<<conjb<<std::endl;
      	  const std::complex<double> factor = (signala_idx >= N) ? ((std::complex<double>*)signala)[signala_idx-N] : 1.0;
      	  acf += factor * conjb;
      	}
      result[tau] = acf;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void xcorr_fftw(std::complex<double> *signala, std::complex<double> *signalb, std::complex<double> *result, int N)
{
  fftw_complex *resulttemp = reinterpret_cast<fftw_complex*>(result);
  
  int N2 = 2 * N - 1;

  fftw_complex * signala_ext = (fftw_complex *) calloc(N2, sizeof(fftw_complex));
  fftw_complex * signalb_ext = (fftw_complex *) calloc(N2, sizeof(fftw_complex));

  // zero padding
  memcpy(signala_ext + (N - 1), signala, sizeof(fftw_complex) * N);
  memcpy(signalb_ext, signalb, sizeof(fftw_complex) * N);

  fftw_complex * outa = fftw_alloc_complex(N2);
  fftw_complex * outb = fftw_alloc_complex(N2);
  fftw_complex * out = fftw_alloc_complex(N2);

  fftw_plan pa = fftw_plan_dft_1d(N2, signala_ext, outa, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan pb = fftw_plan_dft_1d(N2, signalb_ext, outb, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan px = fftw_plan_dft_1d(N2, out, resulttemp, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(pa);
  fftw_execute(pb);

  std::complex<double> *outatemp = reinterpret_cast<std::complex<double>*>(outa);
  std::complex<double> *outbtemp = reinterpret_cast<std::complex<double>*>(outb);
  std::complex<double> *outtemp = reinterpret_cast<std::complex<double>*>(out);
  
  std::complex<double> scale = 1.0/(2 * N -1);
  for (int i = 0; i < N2; ++i)
    outtemp[i] = outatemp[i] * std::conj(outbtemp[i]) * scale;

  fftw_execute(px);

  fftw_destroy_plan(pa);
  fftw_destroy_plan(pb);
  fftw_destroy_plan(px);

  fftw_free(out);
  fftw_free(outa);
  fftw_free(outb);

  fftw_cleanup();

  free(signala_ext);
  free(signalb_ext);
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void xcorr_fftw_r2c(double *signala, double *signalb, std::complex<double> *result, int N)
{
  fftw_complex *resulttemp = reinterpret_cast<fftw_complex*>(result);
  
  int N2 = 2 * N - 1;

  double * signala_ext = (double *) calloc(N2, sizeof(double));
  double * signalb_ext = (double *) calloc(N2, sizeof(double));

  // zero padding
  memcpy(signala_ext + (N - 1), signala, sizeof(double) * N);
  memcpy(signalb_ext, signalb, sizeof(double) * N);

  fftw_complex * outa = fftw_alloc_complex(N2);
  fftw_complex * outb = fftw_alloc_complex(N2);
  fftw_complex * out = fftw_alloc_complex(N2);

  fftw_plan pa = fftw_plan_dft_r2c_1d(N2, signala_ext, outa, FFTW_ESTIMATE);
  fftw_plan pb = fftw_plan_dft_r2c_1d(N2, signalb_ext, outb, FFTW_ESTIMATE);
  fftw_plan px = fftw_plan_dft_1d(N2, out, resulttemp, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(pa);
  fftw_execute(pb);

  std::complex<double> *outatemp = reinterpret_cast<std::complex<double>*>(outa);
  std::complex<double> *outbtemp = reinterpret_cast<std::complex<double>*>(outb);
  std::complex<double> *outtemp = reinterpret_cast<std::complex<double>*>(out);
  
  std::complex<double> scale = 1.0/(2 * N -1);
  for (int i = 0; i < N2; ++i)
    outtemp[i] = outatemp[i] * std::conj(outbtemp[i]) * scale;

  fftw_execute(px);

  fftw_destroy_plan(pa);
  fftw_destroy_plan(pb);
  fftw_destroy_plan(px);

  fftw_free(out);
  fftw_free(outa);
  fftw_free(outb);

  fftw_cleanup();

  free(signala_ext);
  free(signalb_ext);
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void xcorr_fftw_r2r(void *signala, void *signalb, void *result, int N)
// {
//   int N2 = 2 * N - 1;
//   complex * result_cmplx = (complex *) calloc(N2, sizeof(complex));
//   xcorr_fftw_r2c(signala, signalb, result_cmplx, N);
//   for (int x = 0; x < 2*N-1; ++x) {
//     double real = creal(result_cmplx[x]);
//     double real_dest = ((double*) result)[x];
//     memcpy(&real_dest, &real, sizeof(double));
//   }
// }




// 
// xcorr.cc ends here
