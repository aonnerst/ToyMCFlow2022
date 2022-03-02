TFile *_file0 = TFile::Open("flowOutput1000Evts10background.root")
new TBrowser
.q
TF1 *fNUE = new TF1("bgGap","[0]*(1-(x > 1.65)*(x < 2.2)*0.5-(x > 0.3)*(x < 0.4)*0.7)",0.0,2.0*TMath::Pi())
fNUE->SetParameter(0,1.0);
fNUE->Draw()
TFile *_file0 = TFile::Open("test.root")
.ls
hBg->Draw()
hSignal->Draw()
.q

