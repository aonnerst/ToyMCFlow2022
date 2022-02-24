#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <vector>
#include "../include/toyflowinputs.h"
#include "../include/rootcommon.h"

using namespace std;

Int_t gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};
const Int_t NPhiHist = 12;
TF1 *fourier[NPhiHist][NC];
TH1D *hPhiEvent[NPhiHist][NC];
void LoadData(TString);
void DrawEbyE(int);


//Main Loop
void PhiEventByEvent(TString infile="../output/toymcflowao_1d0eb42_10k.root")
{
    LoadData(infile);
//for loop of centrality
    for (int ic=0; ic<NC; ic++){
        DrawEbyE(ic);
    }
}

void LoadData(TString inputname)
{
    TFile *fIn = TFile::Open(inputname,"read");
    fIn->Print();

    for (Int_t i=0; i<=(NPhiHist-1); i++){
        for (Int_t ic=0; ic<NC; ic++){
            cout<<Form("fourierC%02d_E%02d",ic,i+1)<<endl;
            fourier[i][ic] = (TF1*)fIn->Get(Form("fourierC%02d_E%02d",ic,i));
            fourier[i][ic] -> SetParameter(0,10);
            hPhiEvent[i][ic]=(TH1D*)fIn->Get(Form("hPhiEvent_C%02d_E%02d",ic,i+1));
        }
    }
}



void DrawEbyE(int ic=0){
    gStyle->SetOptStat(0);
    TCanvas *cPhi = new TCanvas("cPhi","cPhi");
    TLegend *legendPhi = new TLegend(0.5,0.6,0.8,0.85,"","brNDC");
    legendPhi->SetTextSize(0.04);legendPhi->SetBorderSize(0);legendPhi->SetFillStyle(0);//legend settings;

    double lowx = 0.,highx=2*TMath::Pi();
    double ly=hPhiEvent[0][ic]->GetMinimum()*0.1,hy=hPhiEvent[0][ic]->GetMaximum()*1.4;
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "#phi", "dN/d#phi",0.45,0.45, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();

    for(Int_t i=0; i<NPhiHist; i++){
        hPhiEvent[i][ic]->Draw("esame");
        hPhiEvent[i][ic]->SetLineColor(gColors[i]);
        hPhiEvent[i][ic]->SetLineWidth(1);
        fourier[i][ic]->SetLineColor(gColors[i]);
        fourier[i][ic]->Draw("same");
        legendPhi->AddEntry(hPhiEvent[i][ic],Form("Event %d",i+1));
    }
    legendPhi -> Draw("same");
    gPad->GetCanvas()->SaveAs(Form("figs/PhiEventByEventC%02d.pdf",ic));
}
    





