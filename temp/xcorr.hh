// xcorr.hh --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 5月 12 20:40:16 2017 (+0800)
// Last-Updated: 五 5月 12 22:32:39 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 5
// URL: http://wuhongyi.cn 

#ifndef _XCORR_H_
#define _XCORR_H_

#include <complex>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void xcorr_timedomain(std::complex<double> *signala, std::complex<double> *signalb, std::complex<double> *result, int N);

void xcorr_fftw(std::complex<double> *signala, std::complex<double> *signalb, std::complex<double> *result, int N);

void xcorr_fftw_r2c(double *signala, double *signalb, std::complex<double> *result, int N);


#endif /* _XCORR_H_ */

// 
// xcorr.hh ends here
