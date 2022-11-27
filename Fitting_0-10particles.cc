#include <iostream>
#include "TF1.h"
#include "TObject.h"
#include "TMath.h"
#include "TF3.h"
#include "TH1.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;

Double_t IntegrandBG(const double * x, const double* p){
  // integrand for boltzman-gibbs blast wave

  double x0 = x[0]; 
  
  double mass = p[0];
  double pT   = p[1];
  double beta = p[2];
  double temp    = p[3];
  
  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta*x0);  
  double arg00 = pT*TMath::SinH(rho0)/temp;
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  return f0;
}

Double_t StaticBGdNdPt(Double_t* x, Double_t* p) {

  // implementation of BGBW (1/pt dNdpt)

  double pT = x[0];;
  

  double mass = p[0];
  double beta = p[1];
  double temp    = p[2];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 4);

  fIntBG->SetParameters(mass, pT, beta, temp);
  double result = fIntBG->Integral(0,1);
  return result*p[3];//*1e30;;

}

void Fitting_0_10particles () {
//TFile *file = new TFile("ppcent05.root","read");
//TGraphErrors* gr = (TGraphErrors*)file->FindObject("Graph;1");
  // BGBW 1/pt dNdpt

   TF1* fLastFunc1 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);
   TF1* fLastFunc2 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);
   TF1* fLastFunc3 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);
   TF1* fLastFunc4 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);
   TF1* fLastFunc5 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);
   TF1* fLastFunc6 = new TF1 ("name", StaticBGdNdPt, 0.55, 1.35, 4.0);

  //fLastFunc->SetParameters(0.139,0.35,0.10,2.0);
  fLastFunc1->SetParameters(1.0,0.472,0.119, 1.0);
  fLastFunc1->SetParNames("mass1", "#beta1", "T1", "norm1");
  fLastFunc1->FixParameter(0,0.13960);
  //fLastFunc->SetParLimits(1,0.2,0.7);
  //fLastFunc->FixParameter(2,0.100);
  fLastFunc1->SetParLimits(1,0.1,0.8);
  fLastFunc1->SetParLimits(2,0.080,0.180);
//fLastFunc1->Draw();
fLastFunc1->SetLineColor(1);
fLastFunc1->SetLineStyle(7);
fLastFunc1->SetLineWidth(1);

   //fLastFunc->SetParameters(0.139,0.35,0.10,2.0);
  fLastFunc2->SetParameters(1.0,0.472,0.119,1.0);
  fLastFunc2->SetParNames("mass2", "#beta2", "T2", "norm2");
  fLastFunc2->FixParameter(0,0.4936);
  //fLastFunc->FixParameter(1,0.50);
  //fLastFunc->FixParameter(3,1);
  fLastFunc2->SetParLimits(1,0.1,0.8);
  fLastFunc2->SetParLimits(2,0.100,0.180);

fLastFunc2->SetLineColor(1);
fLastFunc2->SetLineStyle(7);
fLastFunc2->SetLineWidth(1);
 
//fLastFunc->SetParameters(0.139,0.50,0.10,2.0);
  fLastFunc3->SetParameters(1.0,0.472,0.119,1.0);
  fLastFunc3->SetParNames("mass3", "#beta3", "T3", "norm3");
  fLastFunc3->FixParameter(0,0.938);
  //fLastFunc->FixParameter(1,0.50);
  //fLastFunc->FixParameter(3,1);
   fLastFunc3->SetParLimits(1,0.1,0.8);
   fLastFunc3->SetParLimits(2,0.100,0.200);
fLastFunc3->SetLineColor(1);
fLastFunc3->SetLineStyle(7);
fLastFunc3->SetLineWidth(1);

 //fLastFunc->SetParameters(0.139,0.35,0.10,2.0);      
  fLastFunc4->SetParameters(1.0,0.472, 0.119, 1.0);
  fLastFunc4->SetParNames("mass4", "#beta4", "T4", "norm4");
  fLastFunc4->FixParameter(0,0.13960);
  //fLastFunc->FixParameter(1,0.50);
  //fLastFunc->FixParameter(3,1);
   fLastFunc4->SetParLimits(2,0.080,0.180);
   fLastFunc5->SetParLimits(1,0.1,0.8);
fLastFunc4->SetLineColor(1);
fLastFunc4->SetLineStyle(7);
fLastFunc4->SetLineWidth(1);

  //fLastFunc->SetParameters(0.139,.35,0.10,2.0);
  fLastFunc5->SetParameters(1.0,0.472,0.119, 1.0);
  fLastFunc5->SetParNames("mass5", "#beta5", "T5", "norm5");
  fLastFunc5->FixParameter(0,0.4936);
  //fLastFunc->FixParameter(1,0.50);
  //fLastFunc->FixParameter(3,1);
   fLastFunc5->SetParLimits(2,0.100,0.200);
   fLastFunc5->SetParLimits(1,0.1,0.8);
fLastFunc5->SetLineColor(1);
fLastFunc5->SetLineStyle(7);
fLastFunc5->SetLineWidth(1);

  fLastFunc6->SetParameters(1.0,0.472,0.119, 1.0);
  fLastFunc6->SetParNames("mass6", "#beta6", "T6", "norm6");
  fLastFunc6->FixParameter(0,0.9386);
  //fLastFunc->FixParameter(1,0.50);
  //fLastFunc->FixParameter(3,1);
   fLastFunc6->SetParLimits(2,0.100,0.200);
   fLastFunc5->SetParLimits(1,0.1,0.8);
fLastFunc6->SetLineColor(1);
fLastFunc6->SetLineStyle(7);
fLastFunc6->SetLineWidth(1);


fstream file1;
file1.open("0-10.dat",ios::out);

TFile *f1 = new TFile("../../DATA/39_NTMAX2.root");


  TH1F *h1 = (TH1F*)f1->Get("hRefMult");
  Int_t evnts = h1->Integral();
    
  int colorwheel[9] = {kGreen,kRed,kBlue,kCyan,kOrange,6,kBlack,30,38};

    
    const int kCentBin = 5;
    Float_t  centbd[kCentBin+1]={-0.5,2.5,3.5,5.5,8.5,9.5};
    const char* CentName[kCentBin] = {"60-80%","40-60%","20-40%","10-20%","0-10%"};

    const int kPtBin = 20;
       //double ptbd[kPtBin+1]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0};
       //double ptbd[kPtBin+1]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4};
	double ptbd[kPtBin+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2};
       // double pt[kPtBin]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.3,2.6,3.0,3.4,3.9,4.4,5.0,6.0};
       //double dpt[kPtBin+1]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.4};
	double dpt[kPtBin+1]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
     Int_t npt = sizeof(ptbd)/sizeof(ptbd[0]) -1;
    
TH1D *hPiPlus_Rebin = new TH1D("hpiplud_Rebin","hpiplus_Rebin", npt, ptbd);
TH1D *hKPlus_Rebin = new TH1D("hKaon_Rebin","hK0s_Rebin", npt, ptbd);
TH1D *hProton_Rebin = new TH1D("hproton_Rebin","hK0s_Rebin", npt, ptbd);
TH1D *hPiMinus_Rebin = new TH1D("hpiminus_Rebin","hK0s_Rebin", npt, ptbd);
TH1D *hKMinus_Rebin = new TH1D("hakaon_Rebin","hK0s_Rebin", npt, ptbd);
TH1D *hAProton_Rebin = new TH1D("haproton_Rebin","hK0s_Rebin", npt, ptbd);
        
    
    Float_t ptbin[kPtBin],pterr[kPtBin];
      for(int i=0;i<kPtBin;i++){
       ptbin[i]=(ptbd[i]+ptbd[i+1])/2.;
       pterr[i]=0;
      }
   // pterr[]=0;
    
    const double dy=0.1;
    const float err = 1e-4;
  const Int_t nDim = 8;
    double x[kPtBin];
       int ptmin[kCentBin][kPtBin], ptmax[kCentBin][kPtBin];
        int ptmin1[kCentBin][kPtBin], ptmax1[kCentBin][kPtBin];
        int ptmin2[kCentBin][kPtBin], ptmax2[kCentBin][kPtBin];
        int ptmin3[kCentBin][kPtBin], ptmax3[kCentBin][kPtBin];
        int ptmin4[kCentBin][kPtBin], ptmax4[kCentBin][kPtBin];
        int ptmin5[kCentBin][kPtBin], ptmax5[kCentBin][kPtBin];



    Float_t int_yield[kCentBin][kPtBin];
    Float_t int_yield1[kCentBin][kPtBin];
    Float_t int_yield2[kCentBin][kPtBin];
    Float_t int_yield3[kCentBin][kPtBin];
    Float_t int_yield4[kCentBin][kPtBin];
    Float_t int_yield5[kCentBin][kPtBin];

     Float_t yield[kCentBin][kPtBin];
     Float_t yield1[kCentBin][kPtBin];
     Float_t yield2[kCentBin][kPtBin];
     Float_t yield3[kCentBin][kPtBin];
     Float_t yield4[kCentBin][kPtBin];
     Float_t yield5[kCentBin][kPtBin];
      
      Float_t staterr[kCentBin][kPtBin];
      Float_t staterr1[kCentBin][kPtBin];
      Float_t staterr2[kCentBin][kPtBin];
      Float_t staterr3[kCentBin][kPtBin];
      Float_t staterr4[kCentBin][kPtBin];
      Float_t staterr5[kCentBin][kPtBin];
    
    

      TString oname = "spectra_piplus_39.txt";
    TString oname1= "Spectra_kPlus_39.txt";
    TString oname2= "Spectra_Proton_39.txt";
    TString oname3= "Spectra_piminus_39.txt";
    TString oname4= "Spectra_kminus_39.txt";
    TString oname5= "Spectra_aproton_39.txt";

    ofstream oyield(oname);
    ofstream oyield1(oname1);
    ofstream oyield2(oname2);
    ofstream oyield3(oname3);
    ofstream oyield4(oname4);
    ofstream oyield5(oname5);

    TH1D *hPiPlus[nDim];
    TH1D *hKPlus[nDim];
    TH1D *hProton[nDim];
    TH1D *hPiMinus[nDim];
    TH1D *hKMinus[nDim];
    TH1D *hAProton[nDim];

   
    for(Int_t i=0;i<kCentBin;i++)
    {
        hPiPlus[i] = (TH1D*)f1->Get(Form("hPiPlus_Spectra_Cen%d",i));
        hKPlus[i] = (TH1D*)f1->Get(Form("hKPlus_Spectra_Cen%d",i));
        hProton[i] = (TH1D*)f1->Get(Form("hProton_Spectra_Cen%d",i));
        hPiMinus[i] = (TH1D*)f1->Get(Form("hPiMinus_Spectra_Cen%d",i));
        hKMinus[i] = (TH1D*)f1->Get(Form("hKMinus_Spectra_Cen%d",i));
        hAProton[i] = (TH1D*)f1->Get(Form("hAProton_Spectra_Cen%d",i));

	oyield<<"#---------"<<CentName[i]<<"---------"<<endl;
        oyield1<<"#---------"<<CentName[i]<<"---------"<<endl;
        oyield2<<"#---------"<<CentName[i]<<"---------"<<endl;
        oyield3<<"#---------"<<CentName[i]<<"---------"<<endl;
        oyield4<<"#---------"<<CentName[i]<<"---------"<<endl;
        oyield5<<"#---------"<<CentName[i]<<"---------"<<endl;


        for(int j=0;j<kPtBin;j++)
        {
            ptmin[i][j] = hPiPlus[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax[i][j] = hPiPlus[i]->GetXaxis()->FindBin(ptbd[j+1]-err);
	    ptmin1[i][j] = hKPlus[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax1[i][j] = hKPlus[i]->GetXaxis()->FindBin(ptbd[j+1]-err);
            ptmin2[i][j] = hProton[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax2[i][j] = hProton[i]->GetXaxis()->FindBin(ptbd[j+1]-err);
            ptmin3[i][j] = hPiMinus[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax3[i][j] = hPiMinus[i]->GetXaxis()->FindBin(ptbd[j+1]-err);
            ptmin4[i][j] = hKMinus[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax4[i][j] = hKMinus[i]->GetXaxis()->FindBin(ptbd[j+1]-err);
            ptmin5[i][j] = hAProton[i]->GetXaxis()->FindBin(ptbd[j]+err);
            ptmax5[i][j] = hAProton[i]->GetXaxis()->FindBin(ptbd[j+1]-err);


            x[j]=(ptbd[j]+ptbd[j+1])*0.5;
            //cout<<"Usman"<<x[j]<<endl;
            
            int_yield[i][j]=hPiPlus[i]->Integral(ptmin[i][j], ptmax[i][j]);
	    int_yield1[i][j]=hKPlus[i]->Integral(ptmin1[i][j], ptmax1[i][j]);
            int_yield2[i][j]=hProton[i]->Integral(ptmin2[i][j], ptmax2[i][j]);
            int_yield3[i][j]=hPiMinus[i]->Integral(ptmin3[i][j], ptmax3[i][j]);
            int_yield4[i][j]=hKMinus[i]->Integral(ptmin4[i][j], ptmax4[i][j]);
            int_yield5[i][j]=hAProton[i]->Integral(ptmin5[i][j], ptmax5[i][j]);
            if(j<14)
            {
                  //yield[i][j]=int_yield[i][j]/2;
                yield[i][j]=int_yield[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
		yield1[i][j]=int_yield1[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield2[i][j]=int_yield2[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield3[i][j]=int_yield3[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield4[i][j]=int_yield4[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield5[i][j]=int_yield5[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());

            }
            else if(j<17)
            {
                    //yield[i][j]=int_yield[i][j]/5;
                yield[i][j]=int_yield[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
		yield1[i][j]=int_yield1[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield2[i][j]=int_yield2[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield3[i][j]=int_yield3[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield4[i][j]=int_yield4[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield5[i][j]=int_yield5[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());

            }
            else
            {
                   // yield[i][j]=int_yield[i][j]/10;
                yield[i][j]=int_yield[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
		yield1[i][j]=int_yield1[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield2[i][j]=int_yield2[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield3[i][j]=int_yield3[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield4[i][j]=int_yield4[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());
                yield5[i][j]=int_yield5[i][j]/x[j]/dy/dpt[j]/evnts/(2*TMath::Pi());

            }
                
            staterr[i][j]= sqrt(int_yield[i][j]) / int_yield[i][j] * yield[i][j];
	    staterr1[i][j]= sqrt(int_yield1[i][j]) / int_yield1[i][j] * yield1[i][j];
            staterr2[i][j]= sqrt(int_yield2[i][j]) / int_yield2[i][j] * yield2[i][j];
            staterr3[i][j]= sqrt(int_yield3[i][j]) / int_yield3[i][j] * yield3[i][j];
            staterr4[i][j]= sqrt(int_yield4[i][j]) / int_yield4[i][j] * yield4[i][j];
            staterr5[i][j]= sqrt(int_yield5[i][j]) / int_yield5[i][j] * yield5[i][j];
            
            hPiPlus_Rebin->SetBinContent(j+1, yield[i][j]);
            hPiPlus_Rebin->SetBinError(j+1, staterr[i][j]);

             hKPlus_Rebin->SetBinContent(j+1, yield1[i][j]);
            hKPlus_Rebin->SetBinError(j+1, staterr1[i][j]);

	     hProton_Rebin->SetBinContent(j+1, yield2[i][j]);
            hProton_Rebin->SetBinError(j+1, staterr2[i][j]);

	     hPiMinus_Rebin->SetBinContent(j+1, yield3[i][j]);
            hPiMinus_Rebin->SetBinError(j+1, staterr3[i][j]);

	     hKMinus_Rebin->SetBinContent(j+1, yield4[i][j]);
            hKMinus_Rebin->SetBinError(j+1, staterr4[i][j]);

	     hAProton_Rebin->SetBinContent(j+1, yield5[i][j]);
            hAProton_Rebin->SetBinError(j+1, staterr5[i][j]);

	    cout<<yield[i][j]<<" "<<staterr[i][j]<<endl;
            cout<<yield1[i][j]<<" "<<staterr1[i][j]<<endl;
            cout<<yield2[i][j]<<" "<<staterr2[i][j]<<endl;
            cout<<yield3[i][j]<<" "<<staterr3[i][j]<<endl;
            cout<<yield4[i][j]<<" "<<staterr4[i][j]<<endl;
            cout<<yield5[i][j]<<" "<<staterr5[i][j]<<endl;

            if(j>=0&&j<kPtBin)oyield<<x[j]<<" "<<yield[i][j]<<" "<<staterr[i][j]<<endl;
            if(j>=0&&j<kPtBin)oyield1<<x[j]<<" "<<yield1[i][j]<<" "<<staterr1[i][j]<<endl;
            if(j>=0&&j<kPtBin)oyield2<<x[j]<<" "<<yield2[i][j]<<" "<<staterr2[i][j]<<endl;
            if(j>=0&&j<kPtBin)oyield3<<x[j]<<" "<<yield3[i][j]<<" "<<staterr3[i][j]<<endl;
            if(j>=0&&j<kPtBin)oyield4<<x[j]<<" "<<yield4[i][j]<<" "<<staterr4[i][j]<<endl;
            if(j>=0&&j<kPtBin)oyield5<<x[j]<<" "<<yield5[i][j]<<" "<<staterr5[i][j]<<endl;
        }
        
    }
    

   TGraph *gr[kCentBin];
   TGraph *gr1[kCentBin];
   TGraph *gr2[kCentBin];
   TGraph *gr3[kCentBin];
   TGraph *gr4[kCentBin];
   TGraph *gr5[kCentBin];

   for(Int_t i=0;i<kCentBin;i++){
	 gr[i] = new TGraphErrors(kPtBin,ptbin,yield[i],0,staterr[i]);
  	 gr1[i] = new TGraphErrors(kPtBin,ptbin,yield1[i],0,staterr1[i]);
        gr2[i] = new TGraphErrors(kPtBin,ptbin,yield2[i],0,staterr2[i]);
        gr3[i] = new TGraphErrors(kPtBin,ptbin,yield3[i],0,staterr3[i]);
        gr4[i] = new TGraphErrors(kPtBin,ptbin,yield4[i],0,staterr4[i]);
        gr5[i] = new TGraphErrors(kPtBin,ptbin,yield5[i],0,staterr5[i]);}
 





////////////////////////////////////
TCanvas *c = new TCanvas("c", "canvas", 900,6000 );
     gROOT->SetStyle("Modern");
   // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1.0);
    pad1->SetFillStyle( 4050 );
    pad1->SetFillColor( 0 );
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();   

gPad->SetLogy();
gPad->SetTickx();
gPad->SetTicky();
   //gStyle->SetCanvasDefH(400);
   //gStyle->SetCanvasDefW(400);
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1111);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.1);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.02);
   gStyle->SetLabelOffset(0.009,"Y");
   //gStyle->SetTitleOffset(1.2,"Y");
   gStyle->SetLabelOffset(0.009,"x");
   gStyle->SetTitleOffset(1.2,"x");
   gStyle->SetTickLength(0.03,"x");
   gStyle->SetTickLength(0.02,"y");
   //gStyle->SetTitleFillColor(11);
   //gStyle->SetTitleOffset(0,"Y");
   gStyle->SetLabelSize(0.04,"Y");
   gStyle->SetLabelSize(0.04,"X");
   gStyle->SetFuncWidth(0);
   gStyle->SetFillColor(11);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetNdivisions(510,"X");
   gStyle->SetNdivisions(507,"Y");
   //TString htitle("K_{s}^{0} Spectra39 AMPT NTMAX150");
   gr[4]->SetTitle("");
   gr[4]->GetXaxis()->SetLimits(0,2.25);
   gr[4]->SetMinimum(1e-6);
   gr[4]->SetMaximum(1e3);

   gr[4]->GetXaxis()->SetTitle("#bf{#it{p}_{T} (GeV/#it{c})}");
   gr[4]->GetYaxis()->SetTitle("#bf{1/(N_{ev}2#pi#it{p}_{T})dN^{2}/d#it{p}_{T}dy (GeV/#it{c})^{-2}}");
   gr[4]->GetYaxis()->SetTitleSize(0.04);
//////////////////////////////////////////////////////
   //hdata->SetTitle("");
   //hdata->GetXaxis()->SetLimits(0,2.25);
   //hdata->SetMinimum(1e-4);
   //hdata->SetMaximum(1e5);

  //hdata->GetXaxis()->SetTitle("#bf{#it{p}_{T} (GeV/#it{c})}");
  //hdata->GetYaxis()->SetTitle("#bf{1/(N_{ev}2#pi#it{p}_{T})dN^{2}/d#it{p}_{T}dy (GeV/#it{c})^{-2}}");

   for(int i=0;i<kCentBin;i++){ // gr[i]->SetMarkerColor(colorwheel[i]);
   if(i==4) gr[i]->SetMarkerColor(kRed);
   if(i==4) gr1[i]->SetMarkerColor(4);
   if(i==4) gr2[i]->SetMarkerColor(kCyan+2);
   if(i==4) gr3[i]->SetMarkerColor(kRed);
   if(i==4) gr4[i]->SetMarkerColor(4);
   if(i==4) gr5[i]->SetMarkerColor(kCyan+2);

/*
   if(i==1) gr[i]->SetMarkerColor(4);
   if(i==2) gr[i]->SetMarkerColor(kCyan+2);
   if(i==3) gr[i]->SetMarkerColor(6);
   if(i==4) gr[i]->SetMarkerColor(kBlack);*/
   }
   for(int i=0;i<kCentBin;i++){
        //if(i==4) gr[i]->SetMarkerSize(2);
        //if(i==3) gr[i]->SetMarkerSize(2);
        //if(i==2) gr[i]->SetMarkerSize(2);
        //if(i==1) gr[i]->SetMarkerSize(2);

        if(i==4) gr[i]->SetMarkerSize(1.5);
        if(i==4) gr1[i]->SetMarkerSize(2);
        if(i==4) gr2[i]->SetMarkerSize(2);
        if(i==4) gr3[i]->SetMarkerSize(1.5);
        if(i==4) gr4[i]->SetMarkerSize(2);
        if(i==4) gr5[i]->SetMarkerSize(2);



	if(i==4) gr[i]->SetMarkerStyle(20);
	if(i==4) gr1[i]->SetMarkerStyle(23);
	if(i==4) gr2[i]->SetMarkerStyle(33);
	if(i==4) gr3[i]->SetMarkerStyle(24);
	if(i==4) gr4[i]->SetMarkerStyle(32);
	if(i==4) gr5[i]->SetMarkerStyle(27);

	//if(i==1) gr[i]->SetMarkerStyle(23);
	//if(i==2) gr[i]->SetMarkerStyle(33);
	//if(i==3) gr[i]->SetMarkerStyle(22);
	//if(i==4) gr[i]->SetMarkerStyle(29);
        //gr[i]->SetMarkerStyle(23+i);
        if(i==4) gr[i]->Draw("AP");
	if(i==4) gr1[i]->Draw("PSAME");
	if(i==4) gr2[i]->Draw("PSAME");
	if(i==4) gr3[i]->Draw("PSAME");
	if(i==4) gr4[i]->Draw("PSAME");
	if(i==4) gr5[i]->Draw("PSAME");
//gr010->Draw("AP");
//gr010->Fit(fLastFunc1,"r");
//hdata->Draw("AEP");
//hdata1->Draw("EP,same");
//hdata2->Draw("EP,same");
//hdata->Fit(fLastFunc1,"r");
//hdata1->Fit(fLastFunc2,"r");
//hdata2->Fit(fLastFunc3,"r");
        //else gr[i]->Draw("PSAME");
        //if(i==4) gr[i]->Fit(fLastFunc1,"r");
	//if(i==3) gr[i]->Fit(fLastFunc2,"r");
	//if(i==2) gr[i]->Fit(fLastFunc3,"r");
	//if(i==1) gr[i]->Fit(fLastFunc4,"r");
	if(i==4) gr[i]->Fit(fLastFunc1,"r");
	if(i==4) gr1[i]->Fit(fLastFunc2,"r");
	if(i==4) gr2[i]->Fit(fLastFunc3,"r");
	if(i==4) gr3[i]->Fit(fLastFunc4,"r");
	if(i==4) gr4[i]->Fit(fLastFunc5,"r");
	if(i==4) gr5[i]->Fit(fLastFunc6,"r");

	//gr[i]->GetFunction("fLastFunc")->SetLineColor(6);

   }


           Double_t mass1 = fLastFunc1->GetParameter(0);
           Double_t B1 = fLastFunc1->GetParameter(1);
           Double_t T1 = fLastFunc1->GetParameter(2);
	   Double_t B1_err1 = fLastFunc1->GetParError(1);
	   Double_t T1_err1 = fLastFunc1->GetParError(2);
           Double_t chi2_1 = fLastFunc1->GetChisquare();
	   Double_t ndf1 =fLastFunc1->GetNDF();


	   Double_t mass2 = fLastFunc2->GetParameter(0);
           Double_t B2 = fLastFunc2->GetParameter(1);
           Double_t T2 = fLastFunc2->GetParameter(2);
           Double_t B2_err1 = fLastFunc2->GetParError(1);
           Double_t T2_err1 = fLastFunc2->GetParError(2);
           Double_t chi2_2 = fLastFunc2->GetChisquare();
           Double_t ndf2 =fLastFunc2->GetNDF();


	   Double_t mass3 = fLastFunc3->GetParameter(0);
           Double_t B3 = fLastFunc3->GetParameter(1);
           Double_t T3 = fLastFunc3->GetParameter(2);
           Double_t B3_err1 = fLastFunc3->GetParError(1);
           Double_t T3_err1 = fLastFunc3->GetParError(2);
           Double_t chi2_3 = fLastFunc3->GetChisquare();
           Double_t ndf3 =fLastFunc3->GetNDF();
	 
	   Double_t mass4 = fLastFunc4->GetParameter(0);
           Double_t B4 = fLastFunc4->GetParameter(1);
           Double_t T4 = fLastFunc4->GetParameter(2);
           Double_t B4_err1 = fLastFunc4->GetParError(1);
           Double_t T4_err1 = fLastFunc4->GetParError(2);
           Double_t chi2_4 = fLastFunc4->GetChisquare();
           Double_t ndf4 =fLastFunc4->GetNDF();
	   
	    Double_t mass5 = fLastFunc5->GetParameter(0);
           Double_t B5 = fLastFunc5->GetParameter(1);
           Double_t T5 = fLastFunc5->GetParameter(2);
           Double_t B5_err1 = fLastFunc5->GetParError(1);
           Double_t T5_err1 = fLastFunc5->GetParError(2);
           Double_t chi2_5 = fLastFunc5->GetChisquare();
           Double_t ndf5 =fLastFunc5->GetNDF();

	   Double_t mass6 = fLastFunc6->GetParameter(0);
           Double_t B6 = fLastFunc6->GetParameter(1);
           Double_t T6 = fLastFunc6->GetParameter(2);
           Double_t B6_err1 = fLastFunc6->GetParError(1);
           Double_t T6_err1 = fLastFunc6->GetParError(2);
           Double_t chi2_6 = fLastFunc6->GetChisquare();
           Double_t ndf6 =fLastFunc6->GetNDF();


	  // file1<<"cent0-10"<<endl;
	   file1<<mass1<<"   "<<B1<<"   "<<B1_err1<<"   "<<T1<<"   "<<T1_err1<<"  "<<chi2_1<<"   "<<ndf1<<endl;

	   //file1<<"cent10-20"<<endl;
           file1<<mass2<<"   "<<B2<<"   "<<B2_err1<<"   "<<T2<<"   "<<T2_err1<<"  "<<chi2_2<<"   "<<ndf2<<endl;

	  // file1<<"cent20-40"<<endl;
           file1<<mass3<<"   "<<B3<<"   "<<B3_err1<<"   "<<T3<<"   "<<T3_err1<<"  "<<chi2_3<<"   "<<ndf3<<endl;

	   //file1<<"cent40-60"<<endl;
           file1<<mass4<<"   "<<B4<<"   "<<B4_err1<<"   "<<T4<<"   "<<T4_err1<<"  "<<chi2_4<<"   "<<ndf4<<endl;

	  // file1<<"cent60-80"<<endl;
           file1<<mass5<<"   "<<B5<<"   "<<B5_err1<<"   "<<T5<<"   "<<T5_err1<<"  "<<chi2_5<<"   "<<ndf5<<endl;
           file1<<mass6<<"   "<<B6<<"   "<<B6_err1<<"   "<<T6<<"   "<<T6_err1<<"  "<<chi2_6<<"   "<<ndf6<<endl;


	   
//   gr1_1->Draw("EP");
//gr1_1->SetMarkerStyle(24);
//gr1_1->SetMarkerColor(15);
//gr1_1->SetLineColor(15);
////gr1_1->SetLineWidth(2);
//gr1_1->SetMarkerSize(1.5);
// gr1_1->Fit(fLastFunc,"r");



TLegend *leg1 = new TLegend(0.2,0.15,0.4,0.2,NULL,"brNDC");
leg1->SetHeader("");

   TLegendEntry *le1 = leg1->AddEntry(fLastFunc1,"Boltz-Gibbs. fit ","l");
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.03);
   leg1->Draw();


TLegend *leg = new TLegend(0.65,0.7,0.85,0.95,NULL,"brNDC");
leg->SetHeader("#bf{ AMPT NTMAX150}","c");

   for(int i=kCentBin-1;i>=0;i--){
        if(i-kCentBin+1==0) leg->AddEntry(gr[i],Form("#pi^{+} %s",CentName[i]),"P");
        if(i-kCentBin+1==0) leg->AddEntry(gr1[i],Form("K^{+} %s",CentName[i]),"P");
        if(i-kCentBin+1==0) leg->AddEntry(gr2[i],Form("p %s",CentName[i]),"P");
        if(i-kCentBin+1==0) leg->AddEntry(gr3[i],Form("#pi^{-} %s",CentName[i]),"P");
        if(i-kCentBin+1==0) leg->AddEntry(gr4[i],Form("k^{-} %s",CentName[i]),"P");
        if(i-kCentBin+1==0) leg->AddEntry(gr5[i],Form("#bar{p} %s",CentName[i]),"P");
   }
        //else leg->AddEntry(gr[i],Form("%s (x10^{%d})",CentName[i],i-kCentBin+1),"P");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   //leg->Draw();

   TLatex *text =new TLatex(0.35 ,0.91,"AuAu #sqrt{#it{s}_{NN}} = 39 GeV");
   //TLatex *text =new TLatex(0.2 ,0.83,"pp at #sqrt{S}= 0.9 TeV ");
   text->SetNDC();
   text->SetTextColor(900);
   text->SetTextSize(0.04);
   text->Draw();



gPad->Update();
gPad->Modified();

c->cd();

/*
auto *tt = new TLatex(0.08,0.3,"1/(N_{ev}2#piP_{T})dN^{2}/dP_{T}dy (GeV/c)^{-2}");
   tt->SetTextAlign(12); tt->SetTextSize(0.04);
   tt->SetTextAngle(90);
   tt->Draw();
auto *tt1 = new TLatex(0.4,0.035,"#it{p}_{T} (GeV/c)");
   tt1->SetTextAlign(12); tt->SetTextSize(0.04);
   tt1->SetTextAngle(0);
   tt1->Draw();


  // TLatex *tex = new TLatex(0.2884417,0.4206471,"Transverse momentum p_{   T} (GeV/c)");
   //tex->SetTextSize(0.8974576);
   //tex->SetLineWidth(2);
   //tex->SetTextFont(22);
   //tex->Draw();


    //TFile *file = new TFile("APMTks39.root","RECREATE");
    //file->cd();
    
    
    //for(int i=0;i<kCentBin;i++)
    //{
   //gr[i]->Write("");
    //}
   //TString ofigname = "auau"+Ecm+"_ks_raw_spectra_"+cuttag+".pdf";
*/
}

