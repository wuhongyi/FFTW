// main.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 5月 12 20:48:58 2017 (+0800)
// Last-Updated: 四 12月 30 18:55:34 2021 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 97
// URL: http://wuhongyi.cn 

// g++ main.cc pkuFFTW.cc `root-config --cflags` -lfftw3 `root-config --glibs` -o 123

#include "pkuFFTW.hh"

#include "TRandom.h"
#include "TRint.h"
#include "TGraph.h"
#include "TBenchmark.h"
#include "TCanvas.h"

#include <complex>
#include <iostream>
#include <vector>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char *argv[])
{

  TRint *theApp = new TRint("Rint", &argc, argv);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  TCanvas *c1 = new TCanvas("c1","",600,400);
  c1->Divide(2/*col*/,3/*raw*/);

  gRandom->SetSeed(0);
  
  int Fs = 128;
  double T = 1.0/Fs;
  int L = 256;

  fftw_complex *in,*out;
  in = Malloc_fftw_complex(L);
  out = Malloc_fftw_complex(L);
  
  double data[1024];
  double aa,bb;
  TGraph *ga = new TGraph();
  TGraph *gb = new TGraph();
  TGraph *g = new TGraph();
  TGraph *gfft = new TGraph();

  
  for (int i = 0; i < L; ++i)
    {
      aa = 7*std::cos(2*3.14159*15*(i*T)-30*3.14159/180);
      bb= 3*std::cos(2*3.14159*40*(i*T)-90*3.14159/180);
      data[i] = 5+aa+bb+gRandom->Uniform();
      // data[i] = 5+7*std::cos(2*3.14159*15*(i*T)-30*3.14159/180);
      in[i][0] = data[i];
      in[i][1] = 0;

      ga->SetPoint(i,i,aa);
      gb->SetPoint(i,i,bb);
      g->SetPoint(i,i,data[i]);
    }

  fftw1d fft1d(L,-1);
  fft1d.ExecuteNormalized(in,out);
  


  
  for (int i = 0; i < L/2; ++i)
    {
      gfft->SetPoint(i,double(Fs)*i/L,std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]));
      // std::cout<<i<<"  "<<std::sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1])<<std::endl;
    }


  



  fftw1d fft1df(L,1);
  fft1df.ExecuteNormalized(out,in);


  TGraph *gg = new TGraph();
  for (int i = 0; i < L; ++i)
    {
      gg->SetPoint(i,i,in[i][0]);
    }

  // out[30][0] = 0;
  // out[30][1] = 0;
  // out[0][0] = 0;
  // out[0][1] = 0;
  // out[80][0] = 0;
  // out[80][1] = 0;

  for (int i = 0; i < L; ++i)
    {
      if(i < 2) continue;
      out[i][0] = 0;
      out[i][1] = 0; 
    }
  
  fftw1d fft1dff(L,1);
  fft1dff.ExecuteNormalized(out,in);  
  TGraph *ggg = new TGraph();
  for (int i = 0; i < L; ++i)
    {
      ggg->SetPoint(i,i,in[i][0]);
    }

  c1->cd(1);
  ga->Draw();
  gb->Draw("same");

  
  
  c1->cd(2);
  g->Draw();

  c1->cd(3);
  gfft->Draw();
  
  c1->cd(4);
  gg->Draw();
  ggg->Draw("same");
  
  c1->cd(5);
  ggg->Draw();
  
  c1->Update();

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
//   TCanvas *c1 = new TCanvas("c1","",600,400);
//   c1->Divide(1/*col*/,5/*raw*/);

//   int N = 1024;
//   int N2 = N*2;

//   double x[10240];
//   double y[10240];
//   double z[10240];
//   double Z[10240];
//   double zz[10240];
  
//   fftw_complex *in1 = Malloc_fftw_complex(N);
//   fftw_complex *in2 = Malloc_fftw_complex(N);
//   fftw_complex *out = Malloc_fftw_complex(N2);
  
//   for (int i = 0; i < N; ++i)
//     {
//       x[i] = 1024-i;  //3*std::sin(0.1*i);//1024-i;
//       y[i] = 0.5*i;  //std::cos(3*0.1*i);//0.5*i;

//       in1[i][0] = x[i];
//       in1[i][1] = 0;
//       in2[i][0] = y[i];
//       in2[i][1] = 0;
//     }
  
//   corr_timedomain corrtime(true);
//   corr_fftw corrfftw(N,true);
//   corrfftw.Execute(in1,in2,zz);
//   // corrfftw.Execute(x,y,zz);
  
//   TGraph *gx = new TGraph();
//   TGraph *gy = new TGraph();
//   TGraph *gg = new TGraph();
//   TGraph *ggg = new TGraph();
//   TGraph *gfft = new TGraph();

  
//   corrtime.corr_n_n(N,x,y,z);
  
//   for (int i = 0; i < N; ++i)
//     {
//       gx->SetPoint(i,i,x[i]);
//       gy->SetPoint(i,i,y[i]);
//       gg->SetPoint(i,i,z[i]);
//       gfft->SetPoint(i,i,zz[i]);
//       std::cout<<z[i]<<"  "<<zz[i]<<std::endl;
//     }

//   corrtime.corr_n_n2(N,x,y,Z);
//   for (int i = -(N-1); i < N; ++i)
//     {
//       ggg->SetPoint(i+N-1,i,Z[i+N-1]);
//       // gfft->SetPoint(i+N-1,i,out[i+N-1][0]);
//     }

// for (int i = 0; i < N; ++i)
//   {
//     std::cout<<i <<"  "<<z[i]<<"  "<<Z[i+N-1]<<std::endl;
//   }
  
//   c1->cd(1);
//   gx->Draw();
//   c1->cd(2);
//   gy->Draw();
  
//   c1->cd(3);
//   gg->Draw();

//   c1->cd(4);
//   ggg->Draw();

//   c1->cd(5);
//   gfft->Draw();

//   // c1->cd(1);
//   // gg->Draw();
//   // c1->cd(2);
//   // gfft->Draw();

  
//   c1->Update();

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // TCanvas *c1 = new TCanvas("c1","",600,400);
  // c1->Divide(1/*col*/,2/*raw*/);

  // int N = 1024;
  // double x[10240];
  // double y[10240];
  // double z[10240];

  // memset(x, 0, sizeof(double) * N);
  // memset(y, 0, sizeof(double) * N);
  // int xx[10240];
  // int yy[10240];
  // double zz[10240];

  // std::vector<int> xv(13);
  // std::vector<int> yv(16);
  // std::map<int,double> xm;
  // std::map<int,double> ym;
  
  // //13
  // x[0] = 1;
  // x[12] = 3;
  // x[99] = 2;
  // x[125] = 1;
  // x[368] = 1;
  // x[669] = 4;
  // x[951] = 1;
  // xm[0] = 1;
  // xm[12] = 3;
  // xm[99] = 2;
  // xm[125] = 1;
  // xm[368] = 1;
  // xm[669] = 4;
  // xm[951] = 1;
  // xx[0] = 0;
  // xx[1] = 12;
  // xx[2] = 12;
  // xx[3] = 12;
  // xx[4] = 99;
  // xx[5] = 99;
  // xx[6] = 125;
  // xx[7] = 368;
  // xx[8] = 669;
  // xx[9] = 669;
  // xx[10] = 669;
  // xx[11] = 669;
  // xx[12] = 951;
  // xv.at(0) = 0;
  // xv.at(1) = 12;
  // xv.at(2) = 12;
  // xv.at(3) = 12;
  // xv.at(4) = 99;
  // xv.at(5) = 99;
  // xv.at(6) = 125;
  // xv.at(7) = 368;
  // xv.at(8) = 669;
  // xv.at(9) = 669;
  // xv.at(10) = 669;
  // xv.at(11) = 669;
  // xv.at(12) = 951;
  
  // //16
  // y[25] = 2;
  // y[55] = 1;
  // y[126] = 1;
  // y[364] = 3;
  // y[444] = 2;
  // y[548] = 1;
  // y[987] = 4;
  // y[1011] = 2;
  // ym[25] = 2;
  // ym[55] = 1;
  // ym[126] = 1;
  // ym[364] = 3;
  // ym[444] = 2;
  // ym[548] = 1;
  // ym[987] = 4;
  // ym[1011] = 2;
  // yy[0] = 25;
  // yy[1] = 25;
  // yy[2] = 55;
  // yy[3] = 126;
  // yy[4] = 364;
  // yy[5] = 364;
  // yy[6] = 364;
  // yy[7] = 444;
  // yy[8] = 444;
  // yy[9] = 548;
  // yy[10] = 987;
  // yy[11] = 987;
  // yy[12] = 987;
  // yy[13] = 987;
  // yy[14] = 1011;
  // yy[15] = 1011;
  // yv.at(0) =25;
  // yv.at(1) =25;
  // yv.at(2) =55;
  // yv.at(3) =126;
  // yv.at(4) =364;
  // yv.at(5) =364;
  // yv.at(6) =364;
  // yv.at(7) =444;
  // yv.at(8) =444;
  // yv.at(9) =548;
  // yv.at(10) =987;
  // yv.at(11) =987;
  // yv.at(12) =987;
  // yv.at(13) =987;
  // yv.at(14) =1011;
  // yv.at(15) =1011;
 
  // corr_fftw corrfftw(N,false);

  // gBenchmark->Start("test1");//计时开始
  // for (int i = 0; i < 100000; ++i)
  //   {
  //     corrfftw.Execute(x,y,z);
  //   }
  // gBenchmark->Show("test1");//计时结束并输出时间

  // corr_timedomain corrtime(false);

  // gBenchmark->Start("test2");//计时开始
  // for (int i = 0; i < 100000; ++i)
  //   {
  //     // corrtime.corr_n_n(N,x,y,zz);
  //     // corrtime.corr_n(13,xx,16,yy,N,zz);
  //     // corrtime.corr_n(&xv,&yv,N,zz);
  //     corrtime.corr_n(&xm,&ym,N,zz);
  //   }
  // gBenchmark->Show("test2");//计时结束并输出时间


  // TGraph *gg = new TGraph();
  // TGraph *ggg = new TGraph();

  // for (int i = 0; i < N; ++i)
  //   {
  //     gg->SetPoint(i,i,z[i]);
  //     ggg->SetPoint(i,i,zz[i]);
  //     // std::cout<<z[i]<<"  "<<zz[i]<<std::endl;
  //   }

  // c1->cd(1);
  // gg->Draw();
  // c1->cd(2);
  // ggg->Draw();

  // c1->Update();
  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  theApp->Run();

  delete theApp;
  
  return 0;
}

// 
// main.cc ends here
