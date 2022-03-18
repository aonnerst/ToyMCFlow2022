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
#include "../include/Filipad.h"

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
	"NUE",
};

TString gr_Names[NMethod]={"SP","TP","EP"};
TGraphErrors *gr_vn_cent[NConfigs][NMethod][NH];
TGraphErrors *gr_vn_true[NH];

TString rootfiles[]={
	//"../results/vnOutput_BF00_Uni_362af09_10000Evts_2holes.root",
	//"../results/vnOutput_BF10_Uni_362af09_10000Evts_2holes.root",
	//"../results/vnOutput_BF50_Uni_362af09_10000Evts_2holes.root",
	"../results/vnOutput_BF00_NUE_362af09_100kEvts_2holes.root",
	"../results/vnOutput_BF10_NUE_362af09_100kEvts_2holes.root",
	"../results/vnOutput_BF50_NUE_362af09_100kEvts_2holes.root",
	//"../results/vnOutput_BF00_NUE_weightTest_100Evts3.root",
	//"../results/vnOutput_BF10_NUE_weightTest_100Evts3.root",
	//"../results/vnOutput_BF50_NUE_weightTest_100Evts3.root",

};

void LoadGraphs();
void DrawVnCentConfig(int i=0, int ib=0, int ih=0);

void ConfigCompVn(){
	LoadGraphs();
	for (int ih=0; ih<NH; ih++){
		for(int im=0; im<NMethod; im++){
			//DrawVnCentConfig(im,0,ih); //Uniform background
			DrawVnCentConfig(im,1,ih); // uniform background with NUE
		}
	}
}


void LoadGraphs(){
	for(int icon=0; icon<NConfigs; icon++){
		TFile *fIn = TFile::Open(rootfiles[icon],"read");
		cout<<"Opening config file: "<< icon << endl;
		for(int i=0; i<NMethod; i++){
			for (int ih=0; ih<NH; ih++){
				gr_vn_cent[icon][i][ih]=(TGraphErrors*)fIn->Get(Form("gr_v%02d_%s_cent",ih+1, gr_Names[i].Data()));
				cout<<gr_vn_cent[icon][i][ih]<< endl;
				gr_vn_true[ih]=(TGraphErrors*)fIn->Get(Form("gr_v%02d_true_cent",ih+1));

			}
		}
	}
}

void DrawVnCentConfig(int i = 0, int ib=0, int ih=0) // i = method, ib= what sort of background
{	
	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("C","canvas",1024,740);
	can->SetLeftMargin(0.15);
	can->SetBottomMargin(0.15);
	can->SetFillStyle(4000);
	Filipad *fpad= new Filipad(1,1.1,0.4,100,100,1.2,5);
	fpad->Draw();
	//====Upper pad
	TPad *p = fpad->GetPad(1); 
	p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
	TLegend *legend = new TLegend(0.2,0.6,0.6,0.85,"","brNDC");
	legend->SetTextSize(0.04);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
	double lowx = -0.5,highx=2.5;
	double ly=-1.2*TMath::MinElement(NC,gr_vn_true[ih]->GetY()),hy=4.5*TMath::MaxElement(NC,gr_vn_true[ih]->GetY());
	//double ly=-0.05,hy=0.04;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "Centrality %", "v_{n}",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
	hfr->Draw();
	

	legend->AddEntry((TObjArray*)NULL,Form("Method %s",gr_Names[i].Data())," ");
	//legend->AddEntry((TObjArray*)NULL,Form("True value of v_{2}: %s=%0.4f,% s=%0.4f, %s=%0.4f",strCentrality[0].Data(),inputVn[1][0],strCentrality[1].Data(),inputVn[1][1],strCentrality[2].Data(),inputVn[1][2])," ");


	gr_vn_true[ih]->SetLineColor(gColors[0]);
	gr_vn_true[ih]->SetMarkerStyle(gMarkers[0]);
	gr_vn_true[ih]->SetMarkerColor(gColors[0]);
	gr_vn_true[ih]->Draw("lpsame");
	legend->AddEntry(gr_vn_true[ih],Form("n=%d true value", ih+1));
	cout<<"creating canvas" << endl;

	for(int icon=0; icon<NConfigs; icon++){
		cout<<"getting graphs" << endl;
		gr_vn_cent[icon][i][ih]->SetLineColor(gColors[icon+1]);
		gr_vn_cent[icon][i][ih]->SetMarkerStyle(gMarkers[icon+1]);
		gr_vn_cent[icon][i][ih]->SetMarkerColor(gColors[icon+1]);
		gr_vn_cent[icon][i][ih]->Draw("lpsame");
		legend->AddEntry(gr_vn_cent[icon][i][ih],Form("n=%d %s", ih+1,ConfigNames[icon].Data()));
		
	}
	

	legend -> Draw("same");

	//====Lower pad, i.e. ratios
	TGraphErrors *grRatio[NConfigs];
	
		
	
	grRatio[0]=GetRatio(gr_vn_cent[0][i][ih],gr_vn_true[ih]);
	grRatio[1]=GetRatio(gr_vn_cent[1][i][ih],gr_vn_true[ih]);
	grRatio[2]=GetRatio(gr_vn_cent[2][i][ih],gr_vn_true[ih]);
	

	
	p = fpad->GetPad(2);
	p->SetTickx(); p->SetGridy(0); p->SetLogx(0), p->SetLogy(0); p->cd();
	TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.2*TMath::MinElement(NC,grRatio[0]->GetY()), 1.5*TMath::MaxElement(NC,grRatio[0]->GetY()));
	//TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.0, 2.0);
	hset( *hfr1, "Centrality %", "Calculations/True value",0.7,0.7, 0.07,0.07, 0.03,0.03, 0.08,0.08, 510,505);
	hfr1->Draw();
	for(int ir=0;ir<NConfigs;ir++){
		for(int icon=1; icon<NConfigs; icon++){
			
			grRatio[ir]->SetLineColor(gColors[ir+1]);
			grRatio[ir]->SetMarkerStyle(gMarkers[ir+1]);
			grRatio[ir]->SetMarkerColor(gColors[ir+1]);
			grRatio[ir]->Draw("lpsame");
		
		}
	}
	
	gPad->GetCanvas()->SaveAs(Form("../figs/CentDepVn_v%2d_%s_%s_2Holes.pdf",ih+1,BGNames[ib].Data(),gr_Names[i].Data()));
}



