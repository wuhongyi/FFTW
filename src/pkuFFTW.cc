// pkuFFTW.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 5月 13 19:11:24 2017 (+0800)
// Last-Updated: 五 5月 26 22:07:01 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 21
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
      out[0][0]/=2.;
      out[0][1]/=2.;
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
  out[0][0]/=2.;
  out[0][1]/=2.;
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


  // memset(signalb_ext, 0, sizeof(fftw_complex) * N2);
  
}

corr_fftw::~corr_fftw()
{
  delete fft1d;
  delete fft1dc2r;

  
}

void corr_fftw::Execute(fftw_complex *in1, fftw_complex *in2, double *result)
{


}

void corr_fftw::Execute(double *in1, double *in2, double *result)
{


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

corr_timedomain::corr_timedomain(bool biased)
{
  Biased = biased;
}

corr_timedomain::~corr_timedomain()
{

}

template<typename T>
void corr_timedomain::corr_n_n(int n, T *in1,T *in2,double *out)
{
  if(Biased)
    {
      for (int i = 0; i < n; ++i)
	{
	  out[i] = 0;
	  for (int j = 0; j <= n-1-i; ++j)
	    {
	      out[i] += in1[j]*in2[j+i];
	    }
	  out[i] /= n;
	}
    }
  else
    {
      for (int i = 0; i < n; ++i)
	{
	  out[i] = 0;
	  for (int j = 0; j <= n-1-i; ++j)
	    {
	      out[i] += in1[j]*in2[j+i];
	    }
	  out[i] /= n-i;
	}
    }
}

template<typename T>
void corr_timedomain::corr_n_n2(int n, T *in1,T *in2,double *out)
{
  int n2 = 2*n-1;

  if(Biased)
    {
      for (int i = 0; i < n2; ++i)
	{
	  double sum = 0;
	  for (int j = 0; j < n-1-std::abs(i-(n-1)); ++j)
	    {
	      sum += in1[j]*in2[j+std::abs(i-(n-1))];
	    }
	  out[i] = sum/n;
	}
    }
  else
    {
      for (int i = 0; i < n2; ++i)
	{
	  double sum = 0;
	  for (int j = 0; j < n-1-std::abs(i-(n-1)); ++j)
	    {
	      sum += in1[j]*in2[j+std::abs(i-(n-1))];
	    }
	  out[i] = sum/(n-std::abs(i-(n-1)));
	}
    }
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
