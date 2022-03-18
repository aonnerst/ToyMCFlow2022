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

using namespace std;
enum{kK0, kK1, kK2, nKL}; // order
TComplex QvectorQC[NH][nKL];

double DeltaPhi(double phi1, double phi2); // relative angle
void CalculateQvectors(vector <double> phiarray, vector <double> phiweight);
TComplex Q(int n, int p);
TComplex Two(int n1, int n2 );


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
	Int_t bNUE = atoi(argv[5]); //
	// Declare variables
	cout<< strCentrality[0]<<endl;

	Int_t NPhiHist = 12; // Number of events for event-by-event phi dist.


	TFile *output = new TFile(outFile,"recreate");
	output->cd();

	//Define uniform function for option B
	TF1 *centSamp = new TF1("centSamp", "[0]",0.0,0.9);
	centSamp->SetParameter(0,1.0);



	//Define histogram for option B ------------------------------------------------
	TH1D *hCentSample = new TH1D("hCentSample","hCentSample",3,-0.1,2.1);

	TH1D *hPhiPsi[NH][NC]; //phi-psi_n (symmetry plane)
	TH1D *hPhiPsiQ[NH][NC]; // Q vector (event plane)
	TH1D *hEventPlane[NH][NC]; //single particle delta phi hist phi-psi_n (symmetry plane)
	TH1D *hEventPlaneEP[NH][NC]; // single particle delta phi hist Q vector (event plane)
	TH1D *h2PCumulantVn[NH][NC]; // cumulant method
	TH1D *hPhiEvent[NPhiHist][NC]; // Event-by-event phi 
	TH1D *hResolution[NH][NC];
	TH1D *hResolutionDist[NH][NC];
	TH1D *hResolutionDistA[NH][NC];
	TH1D *hBgPhi = new TH1D("hBgPhi","hBgPhi",200, 0.0, 2.0*TMath::Pi()); 
	TH1D *hSignalPhi = new TH1D("hSignalPhi","hSignalPhi",200, 0.0, 2.0*TMath::Pi());
	TH1D *hInclusivePhi = new TH1D("hInclusivePhi","hInclusivePhi",200, 0.0, 2.0*TMath::Pi());
	//TH2D *hPhicorrNUE = new TH2D("hPhicorrNUE","hPhicorrNUE",200, 0.0, 2.0*TMath::Pi(),200,0.0,2.0);
	

	TF1 *uniform[NH]; // uniform distribution of psi for each harmonic
	//range 0 to 2*pi
	for (Int_t ih=0; ih<NH; ih++){

		//------Symmetry planes random ------
		uniform[ih]= new TF1(Form("uniform%02d",ih+1),"[0]", -1*TMath::Pi()/(ih+1), 1.0*TMath::Pi()/(ih+1));
		uniform[ih]->SetParameter(0,1.);

		for (Int_t ic=0; ic<NC; ic++){
			//-----Histograms---------
			hEventPlane[ih][ic]     = new TH1D(Form("hEventPlaneC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hEventPlaneEP[ih][ic]   = new TH1D(Form("hEventPlaneEPC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			h2PCumulantVn[ih][ic]    = new TH1D(Form("h2PCumulantVnC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hPhiPsi[ih][ic]         = new TH1D(Form("hPhiPsiC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hPhiPsiQ[ih][ic]        = new TH1D(Form("hPhiPsiQC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hResolution[ih][ic]     = new TH1D(Form("hResolutionC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-100, 100);
			hResolutionDist[ih][ic] = new TH1D(Form("hResolutionDistC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-10, 10);
			hResolutionDistA[ih][ic]= new TH1D(Form("hResolutionDistAC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-10, 10);
		}


	}

	TH1D *hDeltaPhiSum[NC];
	for (Int_t ic=0; ic<NC; ic++){
		hDeltaPhiSum[ic] = new TH1D(Form("hDeltaPhiSum_C%02d",ic),Form("%s",strCentrality[ic].Data()),200, 0.0, 2.0*TMath::Pi());
	}
	//event-by-event
	for (Int_t iPhiEvt=0; iPhiEvt<NPhiHist; iPhiEvt++){
		for (Int_t ic=0; ic<NC; ic++){
			hPhiEvent[iPhiEvt][ic] = new TH1D(Form("hPhiEvent_C%02d_E%02d",ic,(iPhiEvt+1)),Form("Event=%02d,%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
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

	//background with two gaps in phi
	TString NUEFormula = "[0]*(1-(x > 1.65)*(x < 2.2)*0.5-(x > 0.3)*(x < 0.4)*0.7)";
	if(!bNUE) NUEFormula="[0]";
	TF1 *fNUE = new TF1("fNUE",NUEFormula,0.0,2.0*TMath::Pi());
	fNUE->SetParameter(0,1.0);
	TF1 *fInvertedNUE = new TF1("fInvertedNUE","1.0/([0]*(1-(x > 1.65)*(x < 2.2)*0.5))",0.0,2.0*TMath::Pi());
	fInvertedNUE->SetParameter(0,1.0);
	TF1 *fAdd = new TF1("fAdd","fNUE*fInvertedNUE");
	//-------Random number needed for removal
	TRandom3 *prng = new TRandom3(random_seed);
	gRandom->SetSeed(random_seed);
	
	
	//---------------------------End of generating pdfs-------------------------------------
    
	int ieout = Nevt/20;
	if (ieout<1) ieout=1;
	TStopwatch timer;
	timer.Start();
	//initializing necessary variables
	Double_t Psi_n[NH]={0.0};// symmetry plane angle


	//Event loop
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
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			
			phiarray.push_back(phi); // generating signal
			hSignalPhi->Fill(phi); 
		}
		Int_t N_bg = Nch*bgfrac;
		for (Int_t t=0; t<N_bg; t++){ //generating background
			phi = bgUniform->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			phiarray.push_back(phi); // if bNUE'is kTRUE
			hBgPhi->Fill(phi); 

		} 
		//-----------end of removal of particles------------------
		Int_t N_tot = phiarray.size(); // getting total amount of tracks
		//Initializing 
		Double_t Qn_x[NH] = {0.0};//Should be 0 since we sum the qvectors
		Double_t Qn_y[NH] = {0.0};//
		TComplex QvectorsEP[NH];
		Double_t Psi_n_EP[NH]={0.0};
		Double_t AngleDiff[NH]={0.0};
		for(int iH=0;iH<NH;iH++) QvectorsEP[iH] = TComplex(0,0);
		//Declareing correction NUE factor
		Double_t corrNUEi= 1.; // i is used for single particle
		for (Int_t t=0; t<N_tot; t++)//track loop
		{
			corrNUEi = 1./fNUE->Eval(phiarray[t],0,0);
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
		}//End of track loop
		//Only after the track loop, must sum over the tracks first
		for (Int_t n=0; n<NH; n++) {
			Psi_n_EP[n]=(1/double(n+1))*TMath::ATan2(Qn_y[n],Qn_x[n]);
		}
		
		// Start EP
		for (Int_t i=0; i<N_tot; i++){ //track loop 2
			corrNUEi = 1./fNUE->Eval(phiarray[i],0,0);
			//Evenplane method calculated vn
			for (Int_t n=0; n<NH; n++) {
				// Q-vector calculated Event plane method using Q-vectors
				hEventPlaneEP[n][ic]->Fill(TMath::Cos((n+1)*(DeltaPhi(phiarray[i], Psi_n_EP[n])))); 
				hPhiPsiQ[n][ic]->Fill(DeltaPhi(phiarray[i], Psi_n_EP[n]));
			}
			
		} // end of track-loop 2
		
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
			h2PCumulantVn[n-1][ic]->Fill(sctwo.Re());
		}
		
	}// End of event loop
	fNUE->Write("fNUE");
	fInvertedNUE->Write("fInvertedNUE");
	fAdd->Write("fAdd");
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
		for(int n=0; n<NH; n++){
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


