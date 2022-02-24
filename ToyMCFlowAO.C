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

double DeltaPhi(double phi1, double phi2); // relative angle

int main(int argc, char **argv)
{

	TROOT root("flow","run mc");

	if ( argc<4 ) {
		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"+  "<<argv[0]<<" <outputFile> <Nevt> <random seed> <bgfrac>"<<endl;
		cout << endl << endl;
		exit(1);
	}

    // CONSTANT
	char *outFile = argv[1];
	Int_t Nevt= atoi(argv[2]);
	Int_t random_seed = atoi(argv[3]);
    Int_t bgfrac = atoi(argv[4]); //
    TRandom *myRandom = new TRandom(random_seed);
	// Declare variables
    cout<< strCentrality[0]<<endl;

    Int_t NPhiHist = 12;
    Double_t Psi_n[NH]={0.0};
    Double_t vn_psi[NH][NC]={{0.}};
    Double_t vn_obs_EP[NH][NC]={{0.}};
    Double_t vn_phi[NH][NC]={{0.}};
    Double_t weight = 1.0;

    TFile *output = new TFile(outFile,"recreate");
    output->cd();
    
	//Define uniform function for option B
    TF1 *centSamp = new TF1("centSamp", "[0]",0.0,0.9);
    centSamp->SetParameter(0,1.0);

    

	//Define histogram for option B
    TH1D *hCentSample = new TH1D("hCentSample","hCentSample",3,-0.1,2.1);

    TF1 *uniform[NH];
	//TF1 *fourier[NH]; 

    TH1D *hEventPlane[NH][NC];
    TH1D *hEventPlaneEP[NH][NC];
    TH1D *hTwoParticle[NH][NC];
    TH1D *hPhiPsi[NH][NC];
    TH1D *hPhiPsiQ[NH][NC];
    TH1D *hPhiEvent[NPhiHist][NC];
    TH1D *hResolution[NH][NC];
    TH1D *hResolutionDist[NH][NC];
    TH1D *hResolutionDistA[NH][NC];
    TString strformula = "[0]*(1";
    for (Int_t ih=0; ih<NH; ih++){
    	strformula += Form("+2*[%d]*TMath::Cos(%d*(x-[%d]))",ih+1,ih+1,NH+ih+1);
    }
    strformula+=")";
    cout<<strformula<<endl;
	//signal
    TF1 *fourier = new TF1("Fourier", strformula, 0.0, 2.0*TMath::Pi());
	//background
    TF1 *bgUniform = new TF1("bgUniform","[0]",0.0, 2.0*TMath::Pi());
    bgUniform->SetParameter(0,50.0);
	//range 0 to 2*pi
    for (Int_t ih=0; ih<NH; ih++){

		//------Symmetry planes random ------
    	uniform[ih]= new TF1(Form("uniform%02d",ih+1),"[0]", -1*TMath::Pi()/(ih+1), 1.0*TMath::Pi()/(ih+1));
    	uniform[ih]->SetParameter(0,1.);

    	for (Int_t ic=0; ic<NC; ic++){
			//-----Histograms---------
    		hEventPlane[ih][ic]     = new TH1D(Form("hEventPlaneC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
    		hEventPlaneEP[ih][ic]   = new TH1D(Form("hEventPlaneEPC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
    		hTwoParticle[ih][ic]    = new TH1D(Form("hTwoParticleC%02dH%02d",ic,ih+1),Form("n=%d,%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
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
    

    for (Int_t iPhiEvt=0; iPhiEvt<NPhiHist; iPhiEvt++){
    	for (Int_t ic=0; ic<NC; ic++){
    		hPhiEvent[iPhiEvt][ic] = new TH1D(Form("hPhiEvent_C%02d_E%02d",ic,(iPhiEvt+1)),Form("Event=%02d,%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
    	}
    }
    int ieout = Nevt/20;
    if (ieout<1) ieout=1;
    TStopwatch timer;
    timer.Start();
	//Event loop
    for (Int_t iEvent=0; iEvent<Nevt; iEvent++)
    {
    	if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%" << endl ;}
		//Sample randomly from centSamp
    	Double_t dice = centSamp->GetRandom();
    	Int_t Nch=0;
    	Int_t ic=0;
    	if(dice>= 0.0 && dice<0.3) ic=0;
    	if(dice>=0.3 && dice<0.6) ic=1;
    	if(dice>=0.6 && dice<=0.9) ic=2;
    	hCentSample->Fill(ic);
    	Nch=inputNch[ic];
		//Get Psi for different harmonics
		for (Int_t n=0; n<NH; n++) Psi_n[n]=uniform[n]->GetRandom();//harmonic loop
			fourier->SetParameter(0,Nch); 
		for (Int_t i=0; i<NH-1; i++)fourier->SetParameter(i+1,inputVn[i][ic]); //Setting the vn parameters
		for (Int_t i=NH; i<2*NH; i++)fourier->SetParameter(i+1,Psi_n[i-NH]); //Setting the Psi parameters

			if(iEvent<NPhiHist)fourier->Write(Form("fourierC%02d_E%02d",ic,iEvent));
		//Initializing 
		Double_t Qn_x[NH] = {0.0};//Should be 0 since we sum the qvectors
		Double_t Qn_y[NH] = {0.0};//
		TComplex QvectorsEP[NH];
		Double_t Psi_n_EP[NH]={0.0};
		Double_t Psi_n_EPQ[NH]={0.0};
		Double_t AngleDiff[NH]={0.0};
		for(int iH=0;iH<NH;iH++) QvectorsEP[iH] = TComplex(0,0);

    	vector <double> phiarray; //pharray is now vector
		for (Int_t t=0; t<Nch; t++) phiarray.push_back(fourier->GetRandom()); // generating signal
			Int_t N_bg = Nch*bgfrac;

		for (Int_t t=0; t<N_bg; t++) phiarray.push_back(bgUniform->GetRandom()); //generating background
			
			Int_t N_tot = phiarray.size();
		for (Int_t t=0; t<N_tot; t++)//track loop
		{
			if(iEvent<NPhiHist) {
				hPhiEvent[iEvent][ic]->Fill(phiarray[t]);
			}
			
			//Harmonic loop
			for (Int_t n=0; n<NH; n++)
			{
				hPhiPsi[n][ic]->Fill(DeltaPhi(phiarray[t],Psi_n[n]));
				vn_psi[n][ic] = TMath::Cos((n+1)*(DeltaPhi(phiarray[t], Psi_n[n]))); //Change to symmetry plane
				hEventPlane[n][ic]->Fill(vn_psi[n][ic]);
				
				// calculating eventplane with Q-vectors
				Qn_x[n] += weight*TMath::Cos((n+1)*phiarray[t]);
				Qn_y[n]+= weight*TMath::Sin((n+1)*phiarray[t]);
				QvectorsEP[n] += TComplex(TMath::Cos((n+1)*phiarray[t]),TMath::Sin((n+1)*phiarray[t]));
			}

		}//End of track loop

		//Only after the track loop, must sum over the tracks first
		for (Int_t n=0; n<NH; n++) Psi_n_EP[n]=(1/double(n+1))*TMath::ATan2(Qn_y[n],Qn_x[n]); 
			for (Int_t n=0; n<NH; n++) Psi_n_EPQ[n] = QvectorsEP[n].Theta()/double(n+1);

		//Two Particle correlation
				for (Int_t i=0; i<N_tot; i++){
			//Evenplane method calculated vn
					for (Int_t n=0; n<NH; n++) {
						hEventPlaneEP[n][ic]->Fill(TMath::Cos((n+1)*(DeltaPhi(phiarray[i], Psi_n_EP[n])))); 
						hPhiPsiQ[n][ic]->Fill(DeltaPhi(phiarray[i], Psi_n_EP[n]));
					}
					for (Int_t j=0; j<N_tot;j++){
						if(i==j) continue;
				hDeltaPhiSum[ic]->Fill(DeltaPhi(phiarray[i], phiarray[j]));//For fitting
				for (Int_t n=0; n<NH; n++){
					vn_phi[n][ic] = TMath::Cos((n+1)*(DeltaPhi(phiarray[i], phiarray[j])));//Analytic solution
					hTwoParticle[n][ic]->Fill(vn_phi[n][ic]);

				}
			}		
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
		
	}// End of event loop



	Double_t MeanArrayTwoParticle[NH][NC]={{0.}};
	Double_t MeanArrayEventPlane[NH][NC]={{0.}};
	Double_t MeanArrayEventPlaneQVec[NH][NC]={{0.}};
	Double_t MeanArrayEvtPlError[NH][NC]={{0.}};
	Double_t MeanArrayEvtPlErrorQvec[NH][NC]={{0.}};
	Double_t MeanArrayTwoPartError[NH][NC]={{0.}};
	Double_t MeanArrayResolution[NH][NC]={{0.}};
	Double_t MeanArrayResolutionError[NH][NC]={{0.}};
	Double_t vn_obs_ERROR[NH][NC]={{0.}};
	Double_t vn_TwoPart[NH][NC]={{0.}};
	Double_t vn_EvtPl[NH][NC]={{0.}};
	Double_t vn_EvtPlQvec[NH][NC]={{0.}};
	Double_t vn_TwoPartError[NH][NC]={{0.}};



	// Calculating the avarage over event
	for (Int_t n=0; n<NH; n++)
	{
		for (Int_t ic=0; ic<NC; ic++){
			MeanArrayTwoParticle[n][ic]=hTwoParticle[n][ic]->GetMean();
			MeanArrayTwoPartError[n][ic]=hTwoParticle[n][ic]->GetMeanError();
			MeanArrayEventPlane[n][ic]=hEventPlane[n][ic]->GetMean();
			MeanArrayEvtPlError[n][ic]=hEventPlane[n][ic]->GetMeanError();
			MeanArrayEventPlaneQVec[n][ic]=hEventPlaneEP[n][ic]->GetMean();
			MeanArrayEvtPlErrorQvec[n][ic]=hEventPlaneEP[n][ic]->GetMeanError();
			MeanArrayResolution[n][ic]=hResolution[n][ic]->GetMean();
			MeanArrayResolutionError[n][ic]=hResolution[n][ic]->GetMeanError();



			//for loop for swapping  vn arrays for vn[ic][n] !transformation!

			vn_TwoPart[n][ic]=TMath::Sqrt(TMath::Abs(MeanArrayTwoParticle[n][ic]));
			vn_TwoPartError[n][ic]=0.5*TMath::Power(MeanArrayTwoPartError[n][ic],-0.5);//Check error propagation in textbook
			vn_EvtPl[n][ic]=MeanArrayEventPlane[n][ic];
			vn_EvtPlQvec[n][ic]=MeanArrayEventPlaneQVec[n][ic]/MeanArrayResolution[n][ic];
			vn_obs_ERROR[n][ic]=TMath::Abs(vn_EvtPlQvec[n][ic])*TMath::Sqrt(TMath::Power(MeanArrayEvtPlErrorQvec[n][ic]/MeanArrayEventPlaneQVec[n][ic],2)+TMath::Power(MeanArrayResolutionError[n][ic]/MeanArrayResolution[n][ic],2));	
		}

	}
	// For Power spectra
	Double_t pvn_obs_ERROR[NC][NH]={{0.}};
	Double_t pvn_TwoPart[NC][NH]={{0.}};
	Double_t pvn_EvtPl[NC][NH]={{0.}};
	Double_t pvn_EvtPlQvec[NC][NH]={{0.}};
	Double_t pvn_TwoPartError[NC][NH]={{0.}};
	Double_t pMeanArrayEvtPlError[NC][NH]={{0.}};

	for (int ic=0; ic<NC; ic++){
		for (int ih=0; ih<NH; ih++){
			pvn_obs_ERROR[ic][ih]        = vn_obs_ERROR[ih][ic];
			pvn_TwoPart[ic][ih]          = vn_TwoPart[ih][ic];
			pvn_EvtPl[ic][ih]            = vn_EvtPl[ih][ic];
			pvn_EvtPlQvec[ic][ih]        = vn_EvtPlQvec[ih][ic];
			pvn_TwoPartError[ic][ih]     = vn_TwoPartError[ih][ic];
			pMeanArrayEvtPlError[ic][ih] = pMeanArrayEvtPlError[ih][ic];
		}
	}

	const int NMethod=3;
	TString gr_Names[NMethod]={"SP","TP","EP"};

	//Fill graphs I drew in the sktech
	double Cent[NC] = {0.,1.,2.};
	double eCent[NC] = {0.,0.,0.};
	TGraphErrors *gr_vn_cent[NMethod][NH];

	for (int ih=0; ih<NH; ih++) gr_vn_cent[0][ih]= new TGraphErrors(NC,Cent,vn_EvtPl[ih],eCent,MeanArrayEvtPlError[ih]); //add error here
	for (int ih=0; ih<NH; ih++) gr_vn_cent[1][ih]= new TGraphErrors(NC,Cent,vn_TwoPart[ih],eCent,vn_TwoPartError[ih]); //add error here
	for (int ih=0; ih<NH; ih++) gr_vn_cent[2][ih]= new TGraphErrors(NC,Cent,vn_EvtPlQvec[ih],eCent,vn_obs_ERROR[ih]); //add error here

	//Fill graphs I drew in the sktech
		double px[NH] = {0.};
	double pxe[NH] = {0.};
	for (int ih=0; ih<NH; ih++){px[ih]=ih;pxe[ih]=0.;}

		TGraphErrors *gr_pvn[NC][NMethod];

	for (int ic=0; ic<NC; ic++) gr_pvn[ic][0]= new TGraphErrors(NH,px,pvn_EvtPl[ic],pxe,pMeanArrayEvtPlError[ic]); //add error here
	for (int ic=0; ic<NC; ic++) gr_pvn[ic][1]= new TGraphErrors(NH,px,pvn_TwoPart[ic],pxe,pvn_TwoPartError[ic]); //add error here
	for (int ic=0; ic<NC; ic++) gr_pvn[ic][2]= new TGraphErrors(NH,px,pvn_EvtPlQvec[ic],pxe,pvn_obs_ERROR[ic]); //add error here
		
		for(int i=0; i<NMethod; i++){
			for (int ih=0; ih<NH; ih++){
				gr_vn_cent[i][ih]->SetTitle(Form("Centrality dependence %s Method", gr_Names[i].Data()));
				gr_vn_cent[i][ih]->Write(Form("gr_v%02d_%s_cent",ih+1, gr_Names[i].Data()));
			}
			for (int ic=0; ic<NC; ic++){
				gr_pvn[ic][i]->SetTitle(Form("n dependence %s Method", gr_Names[i].Data()));
				gr_pvn[ic][i]->Write(Form("gr_pv%02d_%s",ic+1, gr_Names[i].Data()));
			}
		}

		
		output->Write();
		output->Close();

		timer.Print();

	}


	double DeltaPhi(double phi1, double phi2) {
  // dphi
		double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
		return res>0 ? res : 2.*TMath::Pi()+res ;
	}