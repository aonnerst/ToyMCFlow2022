#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphPainter.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <vector>
#include "../include/toyflowinputs.h"
#include "../include/rootcommon.h"

using namespace std;

int gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
int gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
int gStyles[]={1,2,3,4,5,6,7,8,9,10};

const int NConfigs = 3;
const int NMethod = 3;
const int Nbg = 2;
TString ConfigNames[NMethod]={
	"No BG",
	//"10 Uniform BG",
	//"50 Uniform BG",
	"10 NUE BG",
	"50 NUE BG",
};
TString BGNames[Nbg]={
	"Uniform",
	"NUE"
};

TString gr_Names[NMethod]={"SP","TP","EP"};
TGraphErrors *gr_vn_cent[NConfigs][NMethod][NH];

TString rootfiles[]={
	//"../vnoutput_BF00_Uni_1000Evts.root",
	//"../vnoutput_BF10_Uni_1000Evts.root",
	//"../vnoutput_BF50_Uni_1000Evts.root",
	"../vnoutput_BF00_NUE_1000Evts.root",
	"../vnoutput_BF10_NUE_1000Evts.root",
	"../vnoutput_BF50_NUE_1000Evts.root",

};

void LoadGraphs();
void DrawVnCentConfig(int i=0, int ib=0);

void ConfigCompVn(){
	LoadGraphs();
	for(int im=0; im<NMethod; im++){
		//DrawVnCentConfig(im,0); //Uniform background
		DrawVnCentConfig(im,1); // uniform background with NUE
	}
}


void LoadGraphs(){
	for(int icon=0; icon<NConfigs; icon++){
		TFile *fIn = TFile::Open(rootfiles[icon],"read");
		cout<<"Opening config file: "<< icon << endl;
		for(int i=0; i<NMethod; i++){
			for (int ih=0; ih<2; ih++){
				gr_vn_cent[icon][i][ih]=(TGraphErrors*)fIn->Get(Form("gr_v%02d_%s_cent",ih+1, gr_Names[i].Data()));

			}
		}
	}
}

void DrawVnCentConfig(int i = 0, int ib=0) // i = method
{	
	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("C","canvas",1024,740);
	can->SetLeftMargin(0.15);
	can->SetBottomMargin(0.15);
	can->SetFillStyle(4000);
	TLegend *legend = new TLegend(0.5,0.6,0.8,0.85,"","brNDC");
	legend->SetTextSize(0.04);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
	double lowx = -0.5,highx=2.5;
	double ly=-0.01,hy=0.4;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "Centrality %", "v_{n}",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
	hfr->Draw();

	legend->AddEntry((TObjArray*)NULL,Form("Method %s",gr_Names[i].Data())," ");

	for(int icon=0; icon<NConfigs; icon++){
		for(int ih=1; ih<2; ih++){
			gr_vn_cent[icon][i][ih]->SetLineColor(gColors[icon]);
			gr_vn_cent[icon][i][ih]->SetMarkerStyle(gMarkers[icon]);
			gr_vn_cent[icon][i][ih]->SetMarkerColor(gColors[icon]);
			gr_vn_cent[icon][i][ih]->Draw("lpsame");
			legend->AddEntry(gr_vn_cent[icon][i][ih],Form("n=%d %s", ih+1,ConfigNames[icon].Data()));
		}
	}
	
	legend -> Draw("same");
	gPad->GetCanvas()->SaveAs(Form("../figs/CentDepVn_%s_%s.pdf",BGNames[ib].Data(),gr_Names[i].Data()));
}


