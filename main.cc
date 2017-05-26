// main.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 5月 12 20:48:58 2017 (+0800)
// Last-Updated: 五 5月 26 22:09:14 2017 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 34
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
  
  TCanvas *c1 = new TCanvas("c1","",600,400);
  c1->Divide(2/*col*/,1/*raw*/);

  gRandom->SetSeed(0);
  
  int Fs = 128;
  double T = 1.0/Fs;
  int L = 256;

  fftw_complex *in,*out;
  in = Malloc_fftw_complex(L);
  out = Malloc_fftw_complex(L);
  
  double data[1024];
  for (int i = 0; i < L; ++i)
    {
      data[i] = 5+7*std::cos(2*3.14159*15*(i*T)-30*3.14159/180)+3*std::cos(2*3.14159*40*(i*T)-90*3.14159/180)+gRandom->Uniform();
      in[i][0] = data[i];
      in[i][1] = 0;
    }

  fftw1d fft1d(L,-1);
  fft1d.ExecuteNormalized(in,out);
  

  TGraph *g = new TGraph();
  // for (int i = 0; i < L/2; ++i)
  //   {
  //     g->SetPoint(i,double(Fs)*i/L,std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]));
  //     // std::cout<<i<<"  "<<std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1])<<std::endl;
  //   }

  for (int i = 0; i < L; ++i)
    {
      g->SetPoint(i,i,data[i]);
    }
  
  c1->cd(1);
  g->Draw();


  fftw1d fft1df(L,1);
  fft1df.ExecuteNormalized(out,in);


  TGraph *gg = new TGraph();
  for (int i = 0; i < L; ++i)
    {
      gg->SetPoint(i,i,in[i][0]);
    }
  c1->cd(2);
  gg->Draw();

  c1->Update();


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  theApp->Run();

  delete theApp;
  
  return 0;
}

// 
// main.cc ends here
