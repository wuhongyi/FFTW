// pkuFFTW.hh --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 5月 13 19:11:04 2017 (+0800)
// Last-Updated: 五 5月 26 22:07:01 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 20
// URL: http://wuhongyi.cn 

#ifndef _PKUFFTW_H_
#define _PKUFFTW_H_

#include <complex>
#include <fftw3.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline fftw_complex* Malloc_fftw_complex(int n)
{
  return (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
}

inline double* Malloc_fftw_real(int n)
{
  return (double*)fftw_alloc_real(n);
}

inline void Free_fftw_complex(fftw_complex* inout)
{
  fftw_free(inout);
}

inline void Free_fftw_real(double *inout)
{
  fftw_free(inout);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class fftwbase
{
public:
  fftwbase();
  virtual ~fftwbase();
  
protected:
  fftw_plan fftwplan;
  bool haveplan;
  
  // 1d
  int N;
  int Sign;
  unsigned Flags;
  fftw_r2r_kind Kind;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class fftw1d : public fftwbase
{
public:
  // n 输入数据点数； sign -1为正变换，+1为逆变换； flags FFTW_MEASURE或FFTW_ESTIMATE
  fftw1d(int n,int sign,unsigned flags = FFTW_MEASURE);// fftw_complex *in, *out;
  virtual ~fftw1d();

  //正变换返回值out数据结构具有对称性,因此只需取前一半
  void Execute(fftw_complex *in, fftw_complex *out);
  void ExecuteNormalized(fftw_complex *in, fftw_complex *out);

  // TODO 添加函数直接得到幅值、相位等

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 实输入数据，复Hermitian输出，正变换。
// 由于实数据的DFT具有 Hermitian对称性，所以只需要计算n/2+1（向下取整）个输出就可以了。比如对于r2c，输入in有n个数据，输出out有floor(n /2)＋1个数据。
class fftw1d_r2c : public fftwbase
{
public:
  fftw1d_r2c(int n,unsigned flags = FFTW_MEASURE);
  virtual ~fftw1d_r2c();

  void Execute(double *in, fftw_complex *out);
  void ExecuteNormalized(double *in, fftw_complex *out);
  
};

// 复Hermitian输入数据，实输出数据，逆变换。
class fftw1d_c2r : public fftwbase
{
public:
  fftw1d_c2r(int n,unsigned flags = FFTW_MEASURE);
  virtual ~fftw1d_c2r();

  void Execute(fftw_complex *in, double *out);
  void ExecuteNormalized(fftw_complex *in, double *out);
  
};


// class fftw1d_r2r : public fftwbase
// {
// public:
//   fftw1d_r2r(int n,fftw_r2r_kind kind,unsigned flags = FFTW_MEASURE);
//   virtual ~fftw1d_r2r();

// };


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class corr_fftw
{
public:
  corr_fftw(int n,bool biased = true);
  virtual ~corr_fftw();

  
  void Execute(fftw_complex *in1, fftw_complex *in2, double *result);// TODO
  void Execute(double *in1, double *in2, double *result);// TODO
  
private:
  int N;
  int N2;
  bool Biased;

  fftw1d *fft1d;
  fftw1d_c2r *fft1dc2r;
  
};


// 对于数据点少的情况，采用时域方法按照公式计算较快！！！
class corr_timedomain
{
public:
  corr_timedomain(bool biased = true);
  virtual ~corr_timedomain();

  // 遍历所有点
  template<typename T>
  void corr_n_n(int n, T *in1,T *in2,double *out);//输出out为n个点
  
  template<typename T>
  void corr_n_n2(int n, T *in1,T *in2,double *out);//输出out为2n-1个点
  
  // 计算稀疏点

  
  // TODO

private:
  bool Biased;
  
};

// 卷积
class conv_fftw
{
public:
  conv_fftw();
  virtual ~conv_fftw();

  // TODO
};




#endif /* _PKUFFTW_H_ */
// 
// pkuFFTW.hh ends here
