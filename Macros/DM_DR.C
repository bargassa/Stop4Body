void DM_DR(){

gStyle->SetOptStat(000000);
gStyle->SetPalette(1);

TChain DM("bdttree");
DM.Add("~cbeiraod/local-area/Stop4Body/nTuples_v2016-10-23/T2DegStop_deltaM30.root");


/////////////////////////////////////////////////////////////
// Define cuts

// Lepton selection
// Muon channel
TCut muon = "(nGoodMu==1)";
// Electron channel
TCut electron = "(nGoodEl==1)";
// Decide
//  TCut singlep = muon;
//  TCut singlep = electron;
TCut singlep = muon||electron;

// Jets
TCut ISRjet = "Jet1Pt > 110.";
TCut dphij1j2 = "DPhiJet1Jet2 < 2.5";

// MET
TCut met    = "Met>300.";

// Preselection
TCut presel = singlep+ISRjet+dphij1j2+met;

TCut s250n220 = "(genStopM==250.)";
TCut s300n270 = "(genStopM==300.)";
TCut s400n370 = "(genStopM==400.)";
TCut s500n470 = "(genStopM==500.)";

/////////////////////////////////////////////////////////////

const Int_t Xbin=50;
const Double_t Xmin=0.;
const Double_t Xmax=10.;

TH1F *dm1 = new TH1F("dm1","#Delta R(l,b)",Xbin,Xmin,Xmax);
DM.Draw("DrJetHBLep>>dm1",presel+s250n220);
double tn = dm1->GetEntries();
dm1->Scale(1./tn);
dm1->SetLineColor(1);
dm1->SetLineWidth(2);
TH1F *dm2 = new TH1F("dm2","#Delta R(l,b)",Xbin,Xmin,Xmax);
DM.Draw("DrJetHBLep>>dm2",presel+s300n270);
double sn = dm2->GetEntries();
dm2->Scale(1./sn);
dm2->SetLineColor(2);
dm2->SetLineWidth(2);
TH1F *dm3 = new TH1F("dm3","#Delta R(l,b)",Xbin,Xmin,Xmax);
DM.Draw("DrJetHBLep>>dm3",presel+s400n370);
double w1n = dm3->GetEntries();
dm3->Scale(1./w1n);
dm3->SetLineColor(3);
dm3->SetLineWidth(2);
TH1F *dm4 = new TH1F("dm4","#Delta R(l,b)",Xbin,Xmin,Xmax);
DM.Draw("DrJetHBLep>>dm4",presel+s500n470);
double w2n = dm4->GetEntries();
dm4->Scale(1./w2n);
dm4->SetLineColor(4);
dm4->SetLineWidth(2);

/////////////////////////////////////////////////////////////

TCanvas *c1 = new TCanvas("c1"," ",1400,1000);
c1->SetFillColor(10);
c1->SetBorderSize(2);

//c1->Divide(3,3);

//c1->cd(1);
c1->SetGridx();
c1->SetGridy();
dm1->Draw();
dm4->Draw("same");
dm3->Draw("same");
dm2->Draw("same");
 TLegend *lg = new TLegend(.65,.65,0.99,1.);
 lg->SetFillColor(0);
 lg->AddEntry(dm1,"Signal (250,220)","l");
 lg->AddEntry(dm2,"Signal (300,270)","l");
 lg->AddEntry(dm3,"Signal (400,370)","l");
 lg->AddEntry(dm4,"Signal (500,470)","l");
 lg->Draw();
 c1->SaveAs("DM_DR.png");

}
