// main.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 5月 12 20:48:58 2017 (+0800)
// Last-Updated: 五 5月 26 23:35:08 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 45
// URL: http://wuhongyi.cn 

#include "pkuFFTW.hh"

#include "TRandom.h"
#include "TRint.h"
#include "TGraph.h"
#include "TBenchmark.h"
#include "TCanvas.h"

#include <complex>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char *argv[])
{

  TRint *theApp = new TRint("Rint", &argc, argv);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // TCanvas *c1 = new TCanvas("c1","",600,400);
  // c1->Divide(2/*col*/,1/*raw*/);

  // gRandom->SetSeed(0);
  
  // int Fs = 128;
  // double T = 1.0/Fs;
  // int L = 256;

  // fftw_complex *in,*out;
  // in = Malloc_fftw_complex(L);
  // out = Malloc_fftw_complex(L);
  
  // double data[1024];
  // for (int i = 0; i < L; ++i)
  //   {
  //     data[i] = 5+7*std::cos(2*3.14159*15*(i*T)-30*3.14159/180)+3*std::cos(2*3.14159*40*(i*T)-90*3.14159/180)+gRandom->Uniform();
  //     in[i][0] = data[i];
  //     in[i][1] = 0;
  //   }

  // fftw1d fft1d(L,-1);
  // fft1d.ExecuteNormalized(in,out);
  

  // TGraph *g = new TGraph();
  // // for (int i = 0; i < L/2; ++i)
  // //   {
  // //     g->SetPoint(i,double(Fs)*i/L,std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]));
  // //     // std::cout<<i<<"  "<<std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1])<<std::endl;
  // //   }

  // for (int i = 0; i < L; ++i)
  //   {
  //     g->SetPoint(i,i,data[i]);
  //   }
  
  // c1->cd(1);
  // g->Draw();


  // fftw1d fft1df(L,1);
  // fft1df.ExecuteNormalized(out,in);


  // TGraph *gg = new TGraph();
  // for (int i = 0; i < L; ++i)
  //   {
  //     gg->SetPoint(i,i,in[i][0]);
  //   }
  // c1->cd(2);
  // gg->Draw();

  // c1->Update();

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  TCanvas *c1 = new TCanvas("c1","",600,400);
  c1->Divide(1/*col*/,5/*raw*/);

  int N = 1024;
  int N2 = N*2;

  double x[10240];
  double y[10240];
  double z[10240];
  double zz[10240];
  
  fftw_complex *in1 = Malloc_fftw_complex(N);
  fftw_complex *in2 = Malloc_fftw_complex(N);
  fftw_complex *out = Malloc_fftw_complex(N2);
  for (int i = 0; i < N; ++i)
    {
      x[i] = 1024-i;  //3*std::sin(0.1*i);//1024-i;
      y[i] = 0.5*i;  //std::cos(3*0.1*i);//0.5*i;

      in1[i][0] = x[i];
      in1[i][1] = 0;
      in2[i][0] = y[i];
      in2[i][1] = 0;
    }
  
  corr_timedomain corrtime(true);
  corr_fftw corrfftw(N,true);
  corrfftw.Execute(in1,in2,zz);
  
  TGraph *gx = new TGraph();
  TGraph *gy = new TGraph();
  TGraph *gg = new TGraph();
  TGraph *ggg = new TGraph();
  TGraph *gfft = new TGraph();

  
  corrtime.corr_n_n(N,x,y,z);
  
  for (int i = 0; i < N; ++i)
    {
      gx->SetPoint(i,i,x[i]);
      gy->SetPoint(i,i,y[i]);
      gg->SetPoint(i,i,z[i]);
      gfft->SetPoint(i,i,zz[i]);
      std::cout<<z[i]<<"  "<<zz[i]<<std::endl;
    }

  corrtime.corr_n_n2(N,x,y,z);
  for (int i = -(N-1); i < N; ++i)
    {
      ggg->SetPoint(i+N-1,i,z[i+N-1]);
      // gfft->SetPoint(i+N-1,i,out[i+N-1][0]);
    }

  c1->cd(1);
  gx->Draw();
  c1->cd(2);
  gy->Draw();
  
  c1->cd(3);
  gg->Draw();

  c1->cd(4);
  ggg->Draw();

  c1->cd(5);
  gfft->Draw();
  
  c1->Update();

  
  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  theApp->Run();

  delete theApp;
  
  return 0;
}

// 
// main.cc ends here
