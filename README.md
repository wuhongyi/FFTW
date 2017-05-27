<!-- README.md --- 
;; 
;; Description: 
;; Author: Hongyi Wu(吴鸿毅)
;; Email: wuhongyi@qq.com 
;; Created: 六 5月 27 10:33:41 2017 (+0800)
;; Last-Updated: 六 5月 27 20:18:36 2017 (+0800)
;;           By: Hongyi Wu(吴鸿毅)
;;     Update #: 3
;; URL: http://wuhongyi.cn -->

# README

```cpp
fftw_complex* Malloc_fftw_complex(int n);
double* Malloc_fftw_real(int n);

void Free_fftw_complex(fftw_complex* inout);
void Free_fftw_real(double *inout)

// 以上函数用来申请、释放所需要的数组的
// 例如

int N = 1024;
fftw_complex *in1 = Malloc_fftw_complex(N);//申请一个长度为1024的复数数组

Free_fftw_complex(in1);//不使用的时候通过该方式给释放掉
```

----

```cpp
// class fftw1d

// n 输入数据点数； sign -1为正变换，+1为逆变换； flags FFTW_MEASURE或FFTW_ESTIMATE
fftw1d(int n,int sign,unsigned flags = FFTW_MEASURE);// fftw_complex *in, *out;

//正变换返回值out数据结构具有对称性,因此只需取前一半
void Execute(fftw_complex *in, fftw_complex *out);//执行变换
void ExecuteNormalized(fftw_complex *in, fftw_complex *out);//执行变换并归一化获得真实幅值（直流分量没有除以2）

void ForwardGetAmplitude(fftw_complex *in,double *out);//正变换获得真实幅值


// 构造函数中，n为变换的点数；sign控制是正变换还是逆变换，其中-1为正变换，+1为逆变换；flags为生成策略，FFTW_MEASURE生成策略可能比较慢，但是执行速度最优，FFTW_ESTIMATE生成策略较快，但是执行速度相对较优。
// 函数 Execute(fftw_complex *in, fftw_complex *out); in为输入量，out为输出量。正变换返回值out数据结构具有对称性,因此只需取前一半

//例如
int L = 1024;
fftw_complex *in,*out;
in = Malloc_fftw_complex(L);
out = Malloc_fftw_complex(L);
for (int i = 0; i < L; ++i)
{
  in[i][0] = 5+7*std::cos(2*3.14159*15*(i*T)-30*3.14159/180);
  in[i][1] = 0;
}

fftw1d fft1d(L,-1);
fft1d.ExecuteNormalized(in,out);
```

----

```cpp
// class fftw1d_r2c
// class fftw1d_c2r
// 以上两个类的使用与fftw1d类似，fftw1d_r2c只能做正变换，fftw1d_c2r只能做逆变换。
```

----

```cpp
// class corr_fftw   采用fftw的相关计算

corr_fftw(int n,bool biased = true);//n为相关计算的点数，biased为true表示有偏，false表示无偏

void Execute(fftw_complex *in1, fftw_complex *in2, double *result);//result长度为n
void Execute(double *in1, double *in2, double *result);//result长度为n 

// 例如
int N = 1024;
fftw_complex *in1 = Malloc_fftw_complex(N);
fftw_complex *in2 = Malloc_fftw_complex(N);
fftw_complex *out = Malloc_fftw_complex(N2);

// in1 in2 赋值

corr_fftw corrfftw1(N,true);
corrfftw1.Execute(in1,in2,out);


double *x = Malloc_fftw_real(N);
double *y = Malloc_fftw_real(N);
double *z = Malloc_fftw_real(N);

// x y 赋值

corr_fftw corrfftw2(N,true);
corrfftw2.Execute(x,y,z);
```

----

```cpp
// class corr_timedomain

corr_timedomain(bool biased = true);//biased为true表示有偏，false表示无偏

template<typename T>
void corr_n_n(int n, T *in1,T *in2,double *out);//输出out为n个点

template<typename T>
void corr_n_n2(int n, T *in1,T *in2,double *out);//输出out为2n-1个点

// 计算稀疏点
void corr_n(std::vector<int> *in1,std::vector<int> *in2,int n,double *out);//in1 in2 为原始数据，为有计数的bin值，如果同一个bin内有多个事件，则该bin值有多个
void corr_n(int n1,int *in1,int n2,int *in2,int n,double *out);
// 例如1024个bin数据。其中 bin0=1 bin14=3 bin168=1 bin658=2 bin1011=1
// 那么in1 中的数据为 0 14 14 14 168 658 658 1011，n1=8



// 例如
int N = 1024;

double *x = Malloc_fftw_real(N);
double *y = Malloc_fftw_real(N);
double *z = Malloc_fftw_real(N);

// x y 赋值

corr_timedomain corrtime(true);
corrtime.corr_n_n(N,x,y,z);

```

----

```cpp
// class conv_fftw

// 未实现

```



<!-- README.md ends here -->
