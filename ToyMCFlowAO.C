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
#include "include/toyflowinputs.h"
#include <TSystem.h> 

using namespace std;
enum{kK0, kK1, kK2, nKL}; // order
TFile *output;
TComplex QvectorQC[NH][nKL];

//OPTIONS TO SET! 
bool calculateTwoP = 0; //1=use 2p method (note: this is a lot slower than the others)
bool b2holes = 0; //1=two holes, 0=1 hole
bool useSampleHisto = 1; //0=use fNUE function, 1=use phi histogram for evaluating the correction terms

double DeltaPhi(double phi1, double phi2); // relative angle
void CalculateQvectors(vector <double> phiarray, vector <double> phiweight);
TComplex Q(int n, int p);
TComplex Two(int n1, int n2 );
void NormalizeSample(TH1D *hist); //to normalize the inclusive phi histo to match fNUE


int main(int argc, char **argv)
{
	TROOT root("flow","run mc");
	if ( argc<4 ) {
		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"+  "<<argv[0]<<" <outputFile> <Nevt> <random seed> <bgfrac> <bNUE>"<<endl;
		cout << endl << endl;
		exit(1);
	}

    // CONSTANT
	char *outFile = argv[1];
	Int_t Nevt= atoi(argv[2]);
	Int_t random_seed = atoi(argv[3]);
	Double_t bgfrac = atof(argv[4]);
	Int_t bNUE = atoi(argv[5]); //if bNUE = false = 0, uses the uniform acceptance for sampling
	// Declare variables
	cout<< strCentrality[0]<<endl;

	Int_t NPhiHist = 12; // Number of events for event-by-event phi dist.

	output = TFile::Open(outFile,"recreate");
	output->cd();

	//Define uniform function for sampling centrality
	TF1 *centSamp = new TF1("centSamp", "[0]",0.0,0.9);
	centSamp->SetParameter(0,1.0);
	TF1 *uniform[NH]; // uniform distribution of psi for each harmonic, range 0 to 2*pi

//Defining histograms------------------------------------------------
	TH1D *hCentSample = new TH1D("hCentSample","hCentSample",3,-0.1,2.1); //sampling centrality
	TH1D *hSample = new TH1D("hSample","hSample",200, 0.0, 2.0*TMath::Pi()); //sampling fNUE histo
	//one for each harmonic
	TH1D *hPhiPsi[NH][NC]; //phi-psi_n (symmetry plane)
	TH1D *hPhiPsiQ[NH][NC]; // Q vector (event plane)
	TH1D *hEventPlane[NH][NC]; //single particle delta phi hist phi-psi_n (symmetry plane)
	TH1D *hEventPlaneEP[NH][NC]; // single particle delta phi hist Q vector (event plane)
	TH1D *h2PCumulantVn[NH][NC]; // cumulant method
	TH1D *hTPcosDeltaPhi[NH][NC]; // two particle delta phi ( phi_i-phi_j)
	TH1D *hResolution[NH][NC];
	TH1D *hResolutionDist[NH][NC];
	TH1D *hResolutionDistA[NH][NC];
	//event-by-event
	TH1D *hPhiEvent[NPhiHist][NC]; // Event-by-event phi 
	TH1D *hBGEvent[NPhiHist][NC]; // Event-by-event background 
	TH1D *hTruthEvent[NPhiHist][NC]; // Event-by-event truth 
	//phi-distibution histos
	TH1D *hBgPhi = new TH1D("hBgPhi","hBgPhi",200, 0.0, 2.0*TMath::Pi()); 
	TH1D *hSignalPhi = new TH1D("hSignalPhi","hSignalPhi",200, 0.0, 2.0*TMath::Pi());
	TH1D *hInclusivePhi = new TH1D("hInclusivePhi","hInclusivePhi",200, 0.0, 2.0*TMath::Pi());
	//TH2D *hPhicorrNUE = new TH2D("hPhicorrNUE","hPhicorrNUE",200, 0.0, 2.0*TMath::Pi(),200,0.0,2.0);

	TH1D *hDeltaPhiSum[NC];

//creation of histograms -------------------------------------------------

	for (Int_t ih=0; ih<NH; ih++){ //harmonic loop
		//------Symmetry planes random ------
		uniform[ih]= new TF1(Form("uniform%02d",ih+1),"[0]", -1*TMath::Pi()/(ih+1), 1.0*TMath::Pi()/(ih+1));
		uniform[ih]->SetParameter(0,1.);

		for (Int_t ic=0; ic<NC; ic++){
			//-----Histograms---------
			hEventPlane[ih][ic]     = new TH1D(Form("hEventPlaneC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hEventPlaneEP[ih][ic]   = new TH1D(Form("hEventPlaneEPC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hTPcosDeltaPhi[ih][ic]    = new TH1D(Form("hTPcosDeltaPhiC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			h2PCumulantVn[ih][ic]    = new TH1D(Form("h2PCumulantVnC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hPhiPsi[ih][ic]         = new TH1D(Form("hPhiPsiC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hPhiPsiQ[ih][ic]        = new TH1D(Form("hPhiPsiQC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hResolution[ih][ic]     = new TH1D(Form("hResolutionC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-100, 100);
			hResolutionDist[ih][ic] = new TH1D(Form("hResolutionDistC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-10, 10);
			hResolutionDistA[ih][ic]= new TH1D(Form("hResolutionDistAC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-10, 10);
		}
	}
	for (Int_t ic=0; ic<NC; ic++){
		hDeltaPhiSum[ic] = new TH1D(Form("hDeltaPhiSum_C%02d",ic),Form("%s",strCentrality[ic].Data()),200, 0.0, 2.0*TMath::Pi());
	}
	//event-by-event
	for (Int_t iPhiEvt=0; iPhiEvt<NPhiHist; iPhiEvt++){
		for (Int_t ic=0; ic<NC; ic++){
			hPhiEvent[iPhiEvt][ic] = new TH1D(Form("hPhiEvent_C%02d_E%02d",ic,(iPhiEvt+1)),Form("Event=%02d,%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
			hBGEvent[iPhiEvt][ic] 	 = new TH1D(Form("hBGEvent_C%02d_E%02d",ic,(iPhiEvt+1)),Form("Event=%02d,%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
			hTruthEvent[iPhiEvt][ic] = new TH1D(Form("hTruthEvent_C%02d_E%02d",ic,(iPhiEvt+1)),Form("Event=%02d,%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
		}
	}
//--------------------------End of histogram ---------------------------------------------------

//-----------------------------Generating pdfs--------------------------------------
	//signal
	TString strformula = "[0]*(1";
	for (Int_t ih=0; ih<NH; ih++){
		strformula += Form("+2*[%d]*TMath::Cos(%d*(x-[%d]))",ih+1,ih+1,NH+ih+1);
	}
	strformula+=")";
	cout<<strformula<<endl;

	TF1 *fourier = new TF1("Fourier", strformula, 0.0, 2.0*TMath::Pi());
	//background
	TF1 *bgUniform = new TF1("bgUniform","[0]",0.0, 2.0*TMath::Pi());
	bgUniform->SetParameter(0,1.0);

	//Making acceptance function fNUE
	TString NUEFormula = "[0]*(1-(x > 1.65)*(x < 2.2)*0.5";
	if(b2holes)NUEFormula+="-(x > 0.3)*(x < 0.4)*0.7";//background with two gaps in phi
	NUEFormula+=")";
	if(!bNUE) NUEFormula="[0]";
	TF1 *fNUE = new TF1("fNUE",NUEFormula,0.0,2.0*TMath::Pi());
	fNUE->SetParameter(0,1.0);

	//-------Random number needed for sampling
	TRandom3 *prng = new TRandom3(random_seed);
	gRandom->SetSeed(random_seed);
	
	//---------------------------End of generating pdfs-------------------------------------
	int ieout = Nevt/20;
	if (ieout<1) ieout=1;
	TStopwatch timer;
	timer.Start();
	Double_t Psi_n[NH]={0.0};// symmetry plane angles for each n

	//Eventloop to fill hSample
	for (Int_t iEvent=0; iEvent<Nevt; iEvent++)
	{
		if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%, filling hSample" << endl ;}
		//--------Sample randomly from centSamp---------
		Double_t dice = centSamp->GetRandom();
		Int_t Nch=0;
		Int_t ic=0;
		if(dice>= 0.0 && dice<0.3) ic=0;
		if(dice>=0.3 && dice<0.6) ic=1;
		if(dice>=0.6 && dice<=0.9) ic=2;
		//-----------End of random sampling of centrality-------------
		Nch=inputNch[ic];
		//Get Psi for different harmonics
		for (Int_t n=0; n<NH; n++) Psi_n[n]=uniform[n]->GetRandom();//harmonic loop
		// Setting parameter values of pdf
		fourier->SetParameter(0,Nch); 
		for (Int_t i=0; i<NH-1; i++){
			fourier->SetParameter(i+1,inputVn[i][ic]); //Setting the vn parameters
		}
		for (Int_t i=NH; i<2*NH; i++){
			fourier->SetParameter(i+1,Psi_n[i-NH]); //Setting the Psi parameters
		}

		vector <double> phiarray; //pharray is now vector
		double phi = -999.;
		for (Int_t t=0; t<Nch; t++){
			phi=fourier->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			phiarray.push_back(phi); // generating signal
		}
		Int_t N_bg = Nch*bgfrac;
		for (Int_t t=0; t<N_bg; t++){ //generating background
			phi = bgUniform->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			phiarray.push_back(phi); // if bNUE'is kTRUE
		}
		Int_t N_tot = phiarray.size();
		for (Int_t t=0; t<N_tot; t++){ 
			hSample->Fill(phiarray[t]); //fill hSample with phiarray following shape of fNUE
		}
	}// End of first event loop
	NormalizeSample(hSample); 

	ieout = Nevt/20;
	if (ieout<1) ieout=1;
	//Psi_n[NH]={0.0};// symmetry plane angle

	//Event loop for all other filling
	for (Int_t iEvent=0; iEvent<Nevt; iEvent++)
	{
		if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%" << endl ;}
		//--------Sample randomly from centSamp---------
		Double_t dice = centSamp->GetRandom();
		Int_t Nch=0;
		Int_t ic=0;
		if(dice>= 0.0 && dice<0.3) ic=0;
		if(dice>=0.3 && dice<0.6) ic=1;
		if(dice>=0.6 && dice<=0.9) ic=2;
		hCentSample->Fill(ic);
		//-----------End of random sampling of centrality-------------
		Nch=inputNch[ic];
		//Get Psi for different harmonics
		for (Int_t n=0; n<NH; n++) Psi_n[n]=uniform[n]->GetRandom();//harmonic loop
		// Setting parameter values of pdf
		fourier->SetParameter(0,Nch); 
		for (Int_t i=0; i<NH-1; i++){
			fourier->SetParameter(i+1,inputVn[i][ic]); //Setting the vn parameters
		}
		for (Int_t i=NH; i<2*NH; i++){
			fourier->SetParameter(i+1,Psi_n[i-NH]); //Setting the Psi parameters
		}
		//-----End of setting parameter values------------------
			if(iEvent<NPhiHist){
				fourier->Write(Form("fourierC%02d_E%02d",ic,iEvent));
			}
		//-----------Putting particle into vector----------------
		vector <double> phiarray; //phiarray is now vector
		vector <double> phiweight; //phiweight is now vector
		double phi = -999.;
		for (Int_t t=0; t<Nch; t++){
			phi=fourier->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//sampling randomly based on fNUE
				continue; //this doesn't happen if bNUE==false, since then fNUE->Eval() is always 1
			phiarray.push_back(phi); // generating signal
			hSignalPhi->Fill(phi); 
			if(iEvent<NPhiHist){
				hTruthEvent[iEvent][ic]->Fill(phi);
			}
		}
		Int_t N_bg = Nch*bgfrac;
		for (Int_t t=0; t<N_bg; t++){ //generating background
			phi = bgUniform->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//sampling randomly based on fNUE
				continue;
			phiarray.push_back(phi); 
			hBgPhi->Fill(phi); 
			if(iEvent<NPhiHist){
				hBGEvent[iEvent][ic]->Fill(phi);
			}
		} 
		//-----------end of placing particles in vectors------------------
		Int_t N_tot = phiarray.size(); // getting total amount of tracks
		//Initializing 
		Double_t Qn_x[NH] = {0.0};//Should be 0 since we sum the qvectors
		Double_t Qn_y[NH] = {0.0};//
		TComplex QvectorsEP[NH];
		Double_t Psi_n_EP[NH]={0.0};
		Double_t AngleDiff[NH]={0.0};
		for(int iH=0;iH<NH;iH++) QvectorsEP[iH] = TComplex(0,0);
		//Declareing correction NUE factor
		Double_t corrNUEi= 1.; Double_t corrNUEj= 1.;// i is used for single particle
		for (Int_t t=0; t<N_tot; t++)// start track loop 1
		{
			if (useSampleHisto){
				corrNUEi = 1./hSample->GetBinContent(phiarray[t]);
				if(corrNUEi > 100) corrNUEi = 1.;
				if(corrNUEi < 0.001) corrNUEi = 1.;
				//cout << "corrNUEi = " << corrNUEi<< "at phi = " << phiarray[t] << endl;
			}
			else corrNUEi = 1./fNUE->Eval(phiarray[t],0,0);
			phiweight.push_back(corrNUEi);
			if(iEvent<NPhiHist) {
				hPhiEvent[iEvent][ic]->Fill(phiarray[t]);
			}
			hInclusivePhi->Fill(phiarray[t]);
			//Harmonic loop
			for (Int_t n=0; n<NH; n++)
			{
				// Analytic Event plane method
				hPhiPsi[n][ic]->Fill(DeltaPhi(phiarray[t],Psi_n[n]));
				hEventPlane[n][ic]->Fill(corrNUEi*TMath::Cos((n+1)*(DeltaPhi(phiarray[t], Psi_n[n]))));
				
				// calculating eventplane with Q-vectors
				Qn_x[n] += corrNUEi*TMath::Cos((n+1)*phiarray[t]);
				Qn_y[n]+= corrNUEi*TMath::Sin((n+1)*phiarray[t]);
				QvectorsEP[n] += TComplex(corrNUEi*TMath::Cos((n+1)*phiarray[t]),corrNUEi*TMath::Sin((n+1)*phiarray[t]));
			}
		}//End of track loop 1
		//This only after the track loop, must sum over the tracks first
		for (Int_t n=0; n<NH; n++) {
			Psi_n_EP[n]=(1/double(n+1))*TMath::ATan2(Qn_y[n],Qn_x[n]);
		}
		
		// Start EP and two-particle correlation
		double correlation[NH][NC] = {0.0};
		double weightsProdSum[NH][NC] = {0.0};
		for (Int_t i=0; i<N_tot; i++){ //start track loop 2
			if (useSampleHisto){
				corrNUEi = 1./hSample->GetBinContent(phiarray[i]); 
			}
			else corrNUEi = 1./fNUE->Eval(phiarray[i],0,0);
			//Event plane method calculated vn
			for (Int_t n=0; n<NH; n++) {
				// Q-vector calculated Event plane method using Q-vectors
				hEventPlaneEP[n][ic]->Fill(TMath::Cos((n+1)*(DeltaPhi(phiarray[i], Psi_n_EP[n])))); 
				hPhiPsiQ[n][ic]->Fill(DeltaPhi(phiarray[i], Psi_n_EP[n]));
			}
			//2 particle correlation method, if option on
			if (!calculateTwoP){continue;}
			for (Int_t j=0; j<N_tot;j++){
				if(i==j) continue;
				if (useSampleHisto){
					corrNUEj = 1./hSample->GetBinContent(phiarray[j]);
				}
				else corrNUEj = 1./fNUE->Eval(phiarray[j],0,0);
				hDeltaPhiSum[ic]->Fill(DeltaPhi(phiarray[i], phiarray[j]));//For fitting
				for (Int_t n=0; n<NH; n++){
					correlation[n][ic] += corrNUEi*corrNUEj*TMath::Cos((n+1)*(DeltaPhi(phiarray[i], phiarray[j])));
					weightsProdSum[n][ic] += corrNUEi*corrNUEj;
					//hTPcosDeltaPhi[n][ic]->Fill(corrNUEi*corrNUEj*TMath::Cos((n+1)*(DeltaPhi(phiarray[i], phiarray[j]))));
				}
			}
		} // end of track loop 2
		for (Int_t n=0; n<NH; n++)
		{
			correlation[n][ic] /= weightsProdSum[n][ic];
			hTPcosDeltaPhi[n][ic]->Fill(correlation[n][ic]);
		}
		//Resolution for every event
		for (Int_t n=0; n<NH; n++)
		{
			AngleDiff[n] = TMath::Cos((n+1)*(DeltaPhi(Psi_n[n], Psi_n_EP[n]))); //Analystical resultion
			//if(n==1) cout <<  Form("n=%d,psi=%.3f, %.3f, %.3f",n,Psi_n[n], Psi_n_EP[n], Psi_n_EPQ[n]) << endl;
			hResolution[n][ic]->Fill(AngleDiff[n]);
			hResolutionDist[n][ic]->Fill(DeltaPhi(Psi_n[n], Psi_n_EP[n]));
			hResolutionDistA[n][ic]->Fill(Psi_n[n]-Psi_n_EP[n]);
		}
		//End EP

		//Start 2P cumulant
		CalculateQvectors(phiarray,phiweight);
		for( int n=1; n<NH+1; n++){
			TComplex sctwo = Two(n, -(n)) / Two(0,0).Re();
			//cout << "sctwo" << sctwo.Re() << endl;
			h2PCumulantVn[n-1][ic]->Fill(sctwo.Re());
		}
		
	}// End of event loop
	fNUE->Write("fNUE");
	output->Write();
	output->Close();
	timer.Print();
}
//---------Member function---------
double DeltaPhi(double phi1, double phi2) {
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	return res>0 ? res : 2.*TMath::Pi()+res ;
}

void CalculateQvectors(vector <double> phiarray, vector <double> phiweight){
	//initiate qvectors as complex zeroes
	for(int n=0; n<NH; n++){
		for(int ik=0; ik<nKL; ++ik){
			QvectorQC[n][ik] = TComplex(0,0);
		}
	} 

	Int_t N_tot = phiarray.size();
	for (int t=0; t<N_tot; t++){
		for(int n=0; n<NH+1; n++){
			Double_t tf = 1.0;
			TComplex q[nKL];
			for(int ik=0; ik<nKL; ik++){
				q[ik] = TComplex(tf*TMath::Cos((n)*phiarray[t]),tf*TMath::Sin((n)*phiarray[t]));
				QvectorQC[n][ik] += q[ik];
				tf *= phiweight[t];
			}
		}
	}
}

TComplex Q(int n, int p){
	// Return QvectorQC
	// Q{-n, p} = Q{n, p}*
	if(n >= 0)
		return QvectorQC[n][p];
	return TComplex::Conjugate(QvectorQC[-n][p]);
}

TComplex Two(int n1, int n2 ){
	// two-particle correlation <exp[i(n1*phi1 + n2*phi2)]>
	//	cout << "TWO FUNCTION " << Q(n1,1) << "*" << Q(n2,1) << " - " << Q(n1+n2 , 2) << endl;
	TComplex two = Q(n1, 1) * Q(n2, 1) - Q( n1+n2, 2);
	return two;
}


//depending on the shape of the histo, this may need to be changed
void NormalizeSample(TH1D *hist)
{
		Int_t binxmin = hist->GetXaxis()->FindBin(TMath::Pi());	
		Int_t binxmax = hist->GetNbinsX();
		double yav =  2.0*hist->Integral(binxmin,binxmax)/binxmax;
		hist->Scale(1/yav);
		return;
}