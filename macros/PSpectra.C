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

Int_t gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};
const int NMethod=3;
const int Nconfigs=3;
TGraphErrors *gr_pvn[Nconfigs][NC][NMethod];
TString gr_Names[NMethod]={"SP","TP","EP"};
Double_t vnIn[NC][NH]={{0.}};
TGraphErrors *gr_vnin[NC];
string strHDummy[] = {"2","3","4","5","6","7","8","9","10","11","12"};//vn
TString ConfigNames[NMethod]={
	"No BG",
	//"10 Uniform BG",
	//"50 Uniform BG",
	"10 NUE BG",
	"50 NUE BG",
};
const int Nbg = 2;
TString BGNames[Nbg]={
	"Uniform",
	"NUE",
};

TString rootfiles[]={
	//"../results/vnOutput_BF00_Uni_9385abc_1000Evts.root ",
	//"../results/vnOutput_BF10_Uni_9385abc_1000Evts.root ",
	//"../results/vnOutput_BF50_Uni_9385abc_1000Evts.root ",
	"../results/vnOutput_BF00_NUE_362af09_100kEvts_2holes.root",
	"../results/vnOutput_BF10_NUE_362af09_100kEvts_2holes.root",
	"../results/vnOutput_BF50_NUE_362af09_100kEvts_2holes.root",
	//"../results/vnOutput_BF00_NUE_weightTest_100Evts3.root",
	//"../results/vnOutput_BF10_NUE_weightTest_100Evts3.root",
	//"../results/vnOutput_BF50_NUE_weightTest_100Evts3.root",

};

void LoadData(); //Loading TGraphs
void DrawPSpectra(int ic=0, int i=0);
//void SaveGraphs(int);

//---Main Function------
void PSpectra()
{
	LoadData();
	for(int ic=0; ic<NC; ic++) {
		for(int i=0; i<NMethod; i++){
			DrawPSpectra(ic,i);
			//SaveGraphs(ic);
		}
		
	}
}

//------Member Functions-------
void LoadData()
{
	for(int icon=0; icon<Nconfigs; icon++){
		TFile *fIn = TFile::Open(rootfiles[icon],"read");
		cout<< "config"<<icon<<endl;
		for(int i=0; i<NMethod; i++){
			cout<<"method"<< i <<endl;
			for (int ic=0; ic<NC; ic++){
				cout<< "centrality"<<ic <<endl;
				gr_pvn[icon][ic][i]=(TGraphErrors*)fIn->Get(Form("gr_pv%02d_%s",ic+1, gr_Names[i].Data()));
				gr_pvn[icon][ic][i]->Print();
				gr_pvn[icon][ic][i]->RemovePoint(0);
			}
		}
	}// end of config loop
	//for loop for drawing input values
	for(int ic=0; ic<NC; ic++){
	  	for (int ih=0; ih<NH; ih++){
	  		vnIn[ic][ih]=inputVn[ih][ic];
	  	}
	}
  	double px[NH] = {0.};
	double pxe[NH] = {0.};
	double vnInError[NC][NH] = {{0.}};
	for (int ih=0; ih<NH; ih++){px[ih]=ih;pxe[ih]=0.;vnInError[0][ih]=vnInError[1][ih]=vnInError[2][ih]=0.;}
	// loop over NMethod
	for(int ic=0; ic<NC; ic++){
		gr_vnin[ic] = new TGraphErrors(NH,px,vnIn[ic],pxe,vnInError[ic]);
		gr_vnin[ic]->Print();
		gr_vnin[ic]->RemovePoint(0);
	} 
	//------end of drawing input values 
}

void DrawPSpectra(int ic=0, int i=0)
{
	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas(Form("C%02d",ic),"canvas",1024,740);
	can->SetFillStyle(4000);
	can->SetLeftMargin(0.15);
   	can->SetBottomMargin(0.15);
   	can->SetLogy(1);
	TLegend *legend = new TLegend(0.5,0.7,0.7,0.9,"","brNDC");
    legend->SetTextSize(0.04);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
	double lowx = 0.5,highx=6.5;
  	double ly=1.2e-5,hy=4e-1;
  	TH2F *hfr = new TH2F("hfr"," ", 7,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  	hset( *hfr, "n+1", "v_{n}",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
  	//Changelabel(hfr,gr_pvn[0][ic],strHDummy);
  	hfr->Draw();
  	legend->AddEntry((TObjArray*)NULL,Form("Centrality %s - %s method",strCentrality[ic].Data(), gr_Names[i].Data())," ");

	gr_vnin[ic]->SetLineColor(kSpring-6);
	gr_vnin[ic]->SetLineWidth(3);
	gr_vnin[ic]->Draw("lsame");
	legend->AddEntry(gr_vnin[ic],"Input","l");

  	for(int icon=0; icon<Nconfigs; icon++){
  	
		gr_pvn[icon][ic][i]->SetLineColor(gColors[icon]);
		gr_pvn[icon][ic][i]->SetMarkerStyle(gMarkers[icon]);
		gr_pvn[icon][ic][i]->SetMarkerColor(gColors[icon]);
		gr_pvn[icon][ic][i]->Draw("plsame");
		legend->AddEntry(gr_pvn[icon][ic][i],Form("%s", ConfigNames[icon].Data()),"pl");
  		
  	}
  	legend -> Draw("same");
   	gPad->GetCanvas()->SaveAs(Form("../figs/PSpectraC%02d%s.pdf",ic,gr_Names[i].Data()));
}
/*
void SaveGraphs(int ic = 0)
{
	TFile *output = new TFile("out_GraphsForRatio.root","recreate");
	for(int im=0; im<NMethod; im++) gr_pvn[ic][im]->Write(Form("gr_pv%02d_%s",ic+1, gr_Names[im].Data()));
	gr_vnin[ic]->Write(Form("gr_vnin_%02d",ic));
	output->Write();
	output->Close();
}
*/
