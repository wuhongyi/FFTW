// pkuFFTW.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 5月 13 19:11:24 2017 (+0800)
// Last-Updated: 六 5月 27 11:14:34 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 43
// URL: http://wuhongyi.cn 

#include "pkuFFTW.hh"
#include <iostream>
#include <cstring>
#include <cmath>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

fftw1d::fftw1d(int n,int sign,unsigned flags)
{
  N = n;
  Sign = sign;
  Flags = flags;

  fftw_complex *in, *out;
  in = Malloc_fftw_complex(N);
  out = Malloc_fftw_complex(N);
 
  fftwplan = fftw_plan_dft_1d(N,in,out,Sign,Flags);
  if(fftwplan) haveplan = true;

  Free_fftw_complex(in);
  Free_fftw_complex(out);
}

fftw1d::~fftw1d()
{

}

void fftw1d::Execute(fftw_complex *in, fftw_complex *out)
{
  if(haveplan) fftw_execute_dft(fftwplan,in,out);
  else std::cout<<"You need "<<std::endl;
}

void fftw1d::ExecuteNormalized(fftw_complex *in, fftw_complex *out)
{
  if(haveplan) fftw_execute_dft(fftwplan,in,out);
  else std::cout<<"You need fftw_plan first."<<std::endl;

  if(Sign == -1)
    {
      for (int i = 0; i < N; ++i)
	{
	  out[i][0] = out[i][0]/N*2;
	  out[i][1] = out[i][1]/N*2;
	}
      // out[0][0]/=2.;
      // out[0][1]/=2.;
    }
  else
    {
      for (int i = 0; i < N; ++i)
	{
	  out[i][0] = out[i][0]/2;
	  out[i][1] = out[i][1]/2;
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

fftw1d_r2c::fftw1d_r2c(int n,unsigned flags)
{
  N = n;
  Flags = flags;

  fftw_complex *out;
  double *in;
  out = Malloc_fftw_complex(N);
  in = Malloc_fftw_real(N);
 
  fftwplan = fftw_plan_dft_r2c_1d(N,in,out,Flags);
  if(fftwplan) haveplan = true;

  Free_fftw_complex(out);
  Free_fftw_real(in);
}

fftw1d_r2c::~fftw1d_r2c()
{

}

void fftw1d_r2c::Execute(double *in, fftw_complex *out)
{
  if(haveplan) fftw_execute_dft_r2c(fftwplan,in,out);
  else std::cout<<"You need fftw_plan first."<<std::endl;
}

void fftw1d_r2c::ExecuteNormalized(double *in, fftw_complex *out)
{
  if(haveplan) fftw_execute_dft_r2c(fftwplan,in,out);
  else std::cout<<"You need fftw_plan first."<<std::endl;

  for (int i = 0; i < N; ++i)
    {
      out[i][0] = out[i][0]/N*2;
      out[i][1] = out[i][1]/N*2;
    }
  // out[0][0]/=2.;
  // out[0][1]/=2.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

fftw1d_c2r::fftw1d_c2r(int n,unsigned flags)
{
  N = n;
  Flags = flags;

  fftw_complex *in;
  double *out;
  in = Malloc_fftw_complex(N);
  out = Malloc_fftw_real(N);
 
  fftwplan = fftw_plan_dft_c2r_1d(N,in,out,Flags);
  if(fftwplan) haveplan = true;

  Free_fftw_complex(in);
  Free_fftw_real(out);
}

fftw1d_c2r::~fftw1d_c2r()
{
  
}

void fftw1d_c2r::Execute(fftw_complex *in, double *out)
{
  if(haveplan) fftw_execute_dft_c2r(fftwplan,in,out);
  else std::cout<<"You need fftw_plan first."<<std::endl;
}

void fftw1d_c2r::ExecuteNormalized(fftw_complex *in, double *out)
{
  if(haveplan) fftw_execute_dft_c2r(fftwplan,in,out);
  else std::cout<<"You need fftw_plan first."<<std::endl;
  
  for (int i = 0; i < N; ++i)
    {
      out[i] = out[i]/2;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

corr_fftw::corr_fftw(int n,bool biased)
{
  N = n;
  N2 = 2*N;
  Biased = biased;
  
  fft1d = new fftw1d(N2,-1);
  fft1dback = new fftw1d(N2,+1);
  
  signala_ext = Malloc_fftw_complex(N2);
  signalb_ext = Malloc_fftw_complex(N2);
  signal_result = Malloc_fftw_complex(N2);
  outa = Malloc_fftw_complex(N2);
  outb = Malloc_fftw_complex(N2);
  out = Malloc_fftw_complex(N2);

  outatemp = reinterpret_cast<std::complex<double>*>(outa);
  outbtemp = reinterpret_cast<std::complex<double>*>(outb);
  outtemp = reinterpret_cast<std::complex<double>*>(out);
  scale = 1.0/N2;
  
}

corr_fftw::~corr_fftw()
{
  delete fft1d;
  delete fft1dback;
  // delete fft1dc2r;

  Free_fftw_complex(signala_ext);
  Free_fftw_complex(signalb_ext);
  Free_fftw_complex(outa);
  Free_fftw_complex(outb);
  Free_fftw_complex(out);
}

void corr_fftw::Execute(fftw_complex *in1, fftw_complex *in2, double *result)
{
  memset(signala_ext, 0, sizeof(fftw_complex) * N2);
  memset(signalb_ext, 0, sizeof(fftw_complex) * N2);

  // zero padding
  memcpy(signala_ext + (N - 1), in1, sizeof(fftw_complex) * N);
  memcpy(signalb_ext, in2, sizeof(fftw_complex) * N);

  fft1d->Execute(signala_ext, outa);
  fft1d->Execute(signalb_ext, outb);

  for (int i = 0; i < N2; ++i)
    outtemp[i] = outatemp[i] * std::conj(outbtemp[i]) * scale * scale ;

  fft1dback->Execute(out, signal_result);

  if(Biased)
    {
      for (int i = 0; i < N; ++i)
	{
	  result[i] = signal_result[N-1-i][0]*2;
	}
    }
  else
    {
      for (int i = 0; i < N; ++i)
	{
	  result[i] = signal_result[N-1-i][0]*2*N/(N-i);
	}
    }
}

void corr_fftw::Execute(double *in1, double *in2, double *result)
{
  memset(signala_ext, 0, sizeof(fftw_complex) * N2);
  memset(signalb_ext, 0, sizeof(fftw_complex) * N2);

  for (int i = 0; i < N; ++i)
    {
      signala_ext[i+N-1][0] = in1[i];
      signalb_ext[i][0] = in2[i];
    }

  fft1d->Execute(signala_ext, outa);
  fft1d->Execute(signalb_ext, outb);

  for (int i = 0; i < N2; ++i)
    outtemp[i] = outatemp[i] * std::conj(outbtemp[i]) * scale * scale ;

  fft1dback->Execute(out, signal_result);

  if(Biased)
    {
      for (int i = 0; i < N; ++i)
	{
	  result[i] = signal_result[N-1-i][0]*2;
	}
    }
  else
    {
      for (int i = 0; i < N; ++i)
	{
	  result[i] = signal_result[N-1-i][0]*2*N/(N-i);
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

corr_timedomain::corr_timedomain(bool biased)
{
  Biased = biased;
}

corr_timedomain::~corr_timedomain()
{

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

conv_fftw::conv_fftw()
{

}

conv_fftw::~conv_fftw()
{

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

fftwbase::fftwbase()
{
  haveplan = false;
  fftwplan = NULL;

  N = -1;
}

fftwbase::~fftwbase()
{
  if(haveplan)
    {
      fftw_destroy_plan(fftwplan);
      fftwplan = NULL;
    }
}


// 
// pkuFFTW.cc ends here
