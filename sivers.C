#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TAxis.h"
#include "TCanvas.h"

using namespace std; 

//global variable & function

bool draw = true;
bool draw2 = false;

#define PI 3.14159265359
double dil_factor = /*0.18*/1.0;
double pol_target = /*0.85*/1.0;
double sivers_asym = 0.03;
double other_asym = 0.1;
double sys_pol = 0.0;
double sys_unpol = 0.0;
double N_mlm[4][200][200];
double N_hist[4] = {0.0};
double center_bin_x[4][100];
double center_bin_y[4][100];
double fc[4] = {4.0,2.0,4.0,4.0};

int N[4] = {54000, 45000, 57000, 55000};
int bin_x = 20;
int bin_y = 20;
int bin_mlm;

double gen_func(double *x, double *par);
double fit_func(double *x, double *par);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t par[], Int_t iflag);
double modulation(double *x, int par);

void sivers()
{
//Histogram for Asymmetry
char name[100];
TH1D *hist_asym[4][3];
double bin_edge[5] = {0.10,0.16,0.19,0.24,0.60};

for (int i=0; i<4; i++)
 {
  for (int j=0; j<3; j++)
   {
     sprintf(name,"Asym_Modulation_%i_Methode_%i",i,j);
     hist_asym[i][j] = new TH1D(name, name, 4, bin_edge);
     hist_asym[i][j]->Sumw2();
   }
 }

// put the Asymmetry
 double Asym[4][4];

 for (int i=0; i<4; i++)
 {
  for (int j=0; j<4; j++)
   {
     if (i==1)
      {
	Asym[i][j] = sivers_asym;
      }
     else
      {
	Asym[i][j] = other_asym;
      }
   }
 }
 

 //Generate pseudo data
 char name1[100];
 TH2D *hist_xt[4];

 for (int ii=0; ii<4; ii++) 
     {
       sprintf(name1,"xTarget_bin_%i",ii);
       hist_xt[ii] = new TH2D(name1,name1, bin_x, -1*PI, PI,bin_y, -1*PI, PI);
       hist_xt[ii]->Sumw2();

       int npar = 5;
       double params[5];

       for (int j=0; j<5; j++)
        {
	 if (j<4)
	  {
	    params[j] = Asym[j][ii];
	  }
	 else
	  {
	    params[j] = N[ii];
	  }
	}
  
        TF2 *f = new TF2("f",gen_func,-1*PI,PI,-1*PI,PI,npar);
        f->SetParameters(params);
	hist_xt[ii]->FillRandom("f", N[ii]);
	
	for (int k=0; k<bin_x; k++)
        {
	 TAxis *xaxis = hist_xt[ii]->GetXaxis();
	 center_bin_x[ii][k] = xaxis->GetBinCenter(k+1);	 
	}

        for (int l=0; l<bin_y; l++)
        {
	 
	 TAxis *yaxis = hist_xt[ii]->GetYaxis();
	 center_bin_y[ii][l] = yaxis->GetBinCenter(l+1);
	}
                 
     }

 // phi-s and phi distribution. In case we need it
 char name2[100];
 char name3[100];
 TH1D *hist_phi_s[4];
 TH1D *hist_phi[4];

 for (int i=0; i<4; i++)
  {
    sprintf(name2,"hist_phi_s_bin_%i",i);
    hist_phi_s[i] = new TH1D(name2,name2, bin_x, -1*PI, PI);
    hist_phi_s[i]->Sumw2();

    sprintf(name3,"hist_phi_bin_%i",i);
    hist_phi[i] = new TH1D(name3,name3, bin_y, -1*PI, PI);
    hist_phi[i]->Sumw2();

    for (int ii = 0; ii < bin_x; ii++)
         {
	   double sum_phi =0;
	   for (int j=0; j<bin_y; j++)
	    {
	      sum_phi = sum_phi + hist_xt[i]->GetBinContent(ii+1,j+1);
	    }
	   hist_phi_s[i]->SetBinContent(ii+1,sum_phi);
	   hist_phi_s[i]->SetBinError(ii+1,TMath::Sqrt(sum_phi));
	 }

    for (int ii = 0; ii < bin_x; ii++)
         {
	   double sum_phi_s =0;
	   for (int j=0; j<bin_y; j++)
	    {
	      sum_phi_s = sum_phi_s + hist_xt[i]->GetBinContent(j+1,ii+1);
	    }
	   hist_phi[i]->SetBinContent(ii+1,sum_phi_s);
	   hist_phi[i]->SetBinError(ii+1,TMath::Sqrt(sum_phi_s));
	 }
  }


 //Standard root-TF2 fit (least square minimization)
 for (int ii=0; ii<4; ii++)
     {

       int npar = 5;
       double params[5];
       for (int j=0; j<5; j++)
        {
	 if (j<4)
	  {
	    params[j] = Asym[j][ii];
	  }
	 else
	  {
	    params[j] = N[ii];
	  }
	}

       TF2 *fit_f = new TF2("fit_f",fit_func,-1*PI,PI,-1*PI,PI,npar);
       fit_f->SetParameters(params);
       fit_f->FixParameter(4,N[ii]);
       hist_xt[ii]->Fit("fit_f","0");
       
       for (int k=0; k<4; k++)
        {
          hist_asym[k][0]->SetBinContent(ii+1,fit_f->GetParameter(k));
	  hist_asym[k][0]->SetBinError(ii+1,fit_f->GetParError(k));
        }
     }

 


 //fourier projection method
 cerr<<"Fourier-Projection Method"<<endl;
 cerr<<"========================="<<endl;

 double An[4][4]={0};
 double error[4][4]={0};


 for (int i=0; i<4; i++)
  {
       
   double sum[4]={0};
   double sum_N=0;
   
   for (int l=0; l<bin_x; l++)
    {
     double phi_s = center_bin_x[i][l];
      for (int m=0; m<bin_y; m++)
       {
	 double phi = center_bin_y[i][m];
         double kin[2] = {phi_s, phi};
         double cont = hist_xt[i]->GetBinContent(l+1,m+1);
         sum_N = sum_N + cont;
         for (int nn=0; nn<4; nn++)
          {
	    sum[nn] = sum[nn] + cont*modulation(kin, nn);
          }
       }
    }

    for (int l=0; l<bin_x; l++)
     {
     double phi_s = center_bin_x[i][l];
      for (int m=0; m<bin_y; m++)
       {
	 double phi = center_bin_y[i][m];
	 double kin[2] = {phi_s, phi};
         double err = hist_xt[i]->GetBinError(l+1,m+1);

         for (int nn=0; nn<4; nn++)
          {
            double add_err = pow(err * (modulation(kin, nn) * sum_N - sum[nn]) / pow(sum_N, 2), 2);
            error[nn][i] = error[nn][i] + add_err;
	    
          }

       }
     }
    
    
    for (int p=0; p<4; p++)
     {
      	An[p][i] = fc[p]*sum[p]/(N[i]*dil_factor*pol_target);
	error[p][i] = fc[p]*TMath::Sqrt(error[p][i]);

	hist_asym[p][1]->SetBinContent(i+1,An[p][i]);
        hist_asym[p][1]->SetBinError(i+1,error[p][i]);

	cerr<<An[p][i]<<" "<<error[p][i]<<endl;
     }
   
   
   cerr<<"=================="<<endl;
 }



 
 //Binned maximum-likelihood method
 cerr<<"Binned-maximum likelihood method"<<endl;
 cerr<<"========================="<<endl;

 for (int i=0; i<4; i++)
  {
   for (int j=0; j<bin_x; j++)
    { 
     for (int k=0; k<bin_y; k++)
       {
         N_mlm[i][j][k] = hist_xt[i]->GetBinContent(j+1,k+1);
       }

    }
  }

 for (int i=0; i<4; i++)
  {
    cerr<<"Start MLM Fitting"<<endl;
    bin_mlm = i;

    int flag;
    gMinuit = new TMinuit(4); 
    gMinuit -> SetPrintLevel(0); 
    gMinuit->SetFCN(fcn);
  
    char name_par[100];
    for (int j=0; j<4; j++)
     {
       sprintf(name_par,"par_%i",j);
       gMinuit->mnparm(j,name_par,0.0025,0.001,-1.0,1.0,flag);
     }
    gMinuit->Migrad();
    
    for (int k=0; k<4; k++)
     {
      double par_value, par_error;
      gMinuit->GetParameter(k, par_value, par_error);
      hist_asym[k][2]->SetBinContent(i+1,par_value);
      hist_asym[k][2]->SetBinError(i+1, par_error);
     }
  }

cerr<<"========="<<endl; 
cerr<<"DONE"<<endl;
cerr<<"========="<<endl; 

//start drawing
 if (draw)
  {
    for (int i=0; i<4; i++)
     {
     for (int j=0; j<3; j++)
       {
	 if (i==0) hist_asym[i][j]->SetTitle("Asymmetry Acos(2#phi); xt; A_{cos(2#phi)}");
         if (i==1) hist_asym[i][j]->SetTitle("Asymmetry Asin(#phi_{s}); xt; A_{sin(#phi_{s})}");
	 if (i==2) hist_asym[i][j]->SetTitle("Asymmetry Acos(2#phi+#phi_{s}); xt; A_{cos(2#phi+#phi_{s})}");
         if (i==3) hist_asym[i][j]->SetTitle("Asymmetry Acos(2#phi-#phi_{s}); xt; A_{cos(2#phi-#phi_{s})}");

         hist_asym[i][j]->SetLineColor(j+1);
	 hist_asym[i][j]->SetMarkerColor(j+1);
	 hist_asym[i][j]->SetMarkerStyle(kFullCircle);
         if (i==1)
          {
           hist_asym[i][j]->SetMinimum(0.);
	   hist_asym[i][j]->SetMaximum(0.05);
          }
	 else
	  {
           hist_asym[i][j]->SetMinimum(0.);
	   hist_asym[i][j]->SetMaximum(0.2);
          }
	 hist_asym[i][j]->SetStats(0);
	 
       }
     }

 //draw all asymmetry superimposed for all methods
    TCanvas *c0[4];
    for (int ii=0; ii<4; ii++)
     {
       c0[ii] = new TCanvas();
       c0[ii]->SetGrid();
       c0[ii]->cd();

       TLine *line1 = new TLine(0.1,0.03,0.6,0.03);
       line1->SetLineStyle(9);
       TLine *line2 = new TLine(0.1,0.1,0.6,0.1);
       line2->SetLineStyle(9);
       
       for (int jj=0; jj<3; jj++)
       {
	  
	  
         if(jj==0) 
          {
           hist_asym[ii][jj]->Draw("E1");
          }
         else 
          { 
           hist_asym[ii][jj]->Draw("E1 SAME");
          } 
       }
        if(ii==1) line1->Draw();
        else line2->Draw();
     }

 //draw only sivers asymmetry individually for each method
    TCanvas *c1[3];
    for (int ii=0; ii<3; ii++)
     {
       c1[ii] = new TCanvas();
       c1[ii]->SetGrid();
       c1[ii]->cd();

       TLine *line1 = new TLine(0.1,0.03,0.6,0.03);
       line1->SetLineStyle(9);
       hist_asym[1][ii]->Draw("E1");
       line1->Draw();
     }

 //draw the phi_s and phi distributions
    TCanvas *c2[4];
    char dist_name[100];
    for (int ii=0; ii<4; ii++)
     {
       c2[ii] = new TCanvas();
       c2[ii]->SetGrid();
       c2[ii]->cd();
       
       sprintf(dist_name,"Event distribution bin-xt = %i; #phi_{s}; #phi",ii+1);
       hist_xt[ii]->SetTitle(dist_name);
       hist_xt[ii]->Draw("COLZ");
       
    }
  }

 if (draw2)
  {
    TCanvas *c3[4];
    char hist_name[100];

    TCanvas *c4[4];
    char hist_name2[100];
    for (int ii=0; ii<4; ii++)
     {
       c3[ii] = new TCanvas();
       c3[ii]->SetGrid();
       c3[ii]->cd();

       sprintf(hist_name,"#phi_{s} distribution bin-xt = %i; #phi_{s}; N",ii+1);
       hist_phi_s[ii]->SetTitle(hist_name);
       hist_phi_s[ii]->SetLineColor(kRed);
       hist_phi_s[ii]->SetMarkerColor(kRed);
       hist_phi_s[ii]->SetMarkerStyle(kFullCircle);
       hist_phi_s[ii]->SetMinimum(1500);
       hist_phi_s[ii]->SetMaximum(3500);
       hist_phi_s[ii]->SetStats(0);
       hist_phi_s[ii]->Draw("E1");

       c4[ii] = new TCanvas();
       c4[ii]->SetGrid();
       c4[ii]->cd();

       sprintf(hist_name2,"#phi distribution bin-xt = %i; #phi; N",ii+1);
       hist_phi[ii]->SetTitle(hist_name2);
       hist_phi[ii]->SetLineColor(kRed);
       hist_phi[ii]->SetMarkerColor(kRed);
       hist_phi[ii]->SetMarkerStyle(kFullCircle);
       hist_phi[ii]->SetMinimum(1500);
       hist_phi[ii]->SetMaximum(3500);
       hist_phi[ii]->SetStats(0);
       hist_phi[ii]->Draw("E1");
       
    }

  }  

}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t par[], Int_t iflag)
{
  
 double nll = 0.0;
 double prob = 0.0;
 double norm = 1.0;

 for (int i =0; i<bin_x; i++)
  {
    for (int j=0; j<bin_y; j++)
    {
      double phi_s = center_bin_x[bin_mlm][i];
      double phi = center_bin_y[bin_mlm][j];
      double N = N_mlm[bin_mlm][i][j];
    
      prob = norm*((1+0.5*par[0]*cos(2*phi))+pol_target*dil_factor*(par[1]*sin(phi_s)+0.5*par[2]*sin(2*phi+phi_s)+0.5*par[3]*sin(2*phi-phi_s)));
      if(prob>0) nll = nll + N*TMath::Log(prob);
    
    }
  }

f = -2*nll;
}

double gen_func(double *x, double *par)
{
 double phi_s = x[0];
 double phi = x[1];

 TRandom rndm;
 rndm.SetSeed(0);

 double sys_random_pol = 1.0 - rndm.Gaus(0,sys_pol);
 double sys_random_unpol = 1.0 - rndm.Gaus(0,sys_unpol);
 double f = par[4]*((1+0.5*sys_random_unpol*par[0]*cos(2*phi))+pol_target*dil_factor*sys_random_pol*(par[1]*sin(phi_s)+0.5*par[2]*sin(2*phi+phi_s)+0.5*par[3]*sin(2*phi-phi_s)));
 return f;
}

double fit_func(double *x, double *par)
{
 double phi_s = x[0];
 double phi = x[1];

 double f = par[4]*((1+0.5*par[0]*cos(2*phi))+pol_target*dil_factor*(par[1]*sin(phi_s)+0.5*par[2]*sin(2*phi+phi_s)+0.5*par[3]*sin(2*phi-phi_s)));
 return f;
}

double modulation(double *x, int par)
{
 double phi_s = x[0];
 double phi = x[1];
 int bin = par;
 
 double f;
 if (bin==0) f = cos(2*phi);
 if (bin==1) f = sin(phi_s);
 if (bin==2) f = sin(2*phi+phi_s);
 if (bin==3) f = sin(2*phi-phi_s);
 
 return f;
}



 

