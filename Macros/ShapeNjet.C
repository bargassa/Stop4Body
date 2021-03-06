void ShapeNjet(){

gStyle->SetOptStat(000000);
gStyle->SetPalette(1);

// TChain l0("bdttree");
// l0.Add("SET3110_DM80/data.root");
TChain TT("bdttree");
TT.Add("SET3110_DM80/TTJets_bdt.root");
TChain W1("bdttree");
W1.Add("SET3110_DM80/Wjets_100to200_bdt.root");
TChain W2("bdttree");
W2.Add("SET3110_DM80/Wjets_200to400_bdt.root");
TChain W3("bdttree");
W3.Add("SET3110_DM80/Wjets_400to600_bdt.root");
TChain W4("bdttree");
W4.Add("SET3110_DM80/Wjets_600to800_bdt.root");
TChain S1("bdttree");
S1.Add("SET3110_DM80/T2DegStop_350_270_bdt.root");


/////////////////////////////////////////////////////////////
// Define cuts


/////////////////////////////////////////////////////////////

const Int_t Xbin=10;
const Double_t Xmin=0.;
const Double_t Xmax=10.;

TH1F *ht = new TH1F("ht","Njet30",Xbin,Xmin,Xmax);
TT.Draw("Njet30>>ht");
double tn = ht->GetEntries();
ht->Scale(1./tn);
ht->SetLineColor(1);
ht->SetLineWidth(2);
TH1F *hs = new TH1F("hs","Njet30",Xbin,Xmin,Xmax);
S1.Draw("Njet30>>hs");
double sn = hs->GetEntries();
hs->Scale(1./sn);
hs->SetLineColor(2);
hs->SetLineWidth(2);
TH1F *hw1 = new TH1F("hw1","Njet30",Xbin,Xmin,Xmax);
W1.Draw("Njet30>>hw1");
double w1n = hw1->GetEntries();
hw1->Scale(1./w1n);
hw1->SetLineColor(3);
hw1->SetLineWidth(2);
TH1F *hw2 = new TH1F("hw2","Njet30",Xbin,Xmin,Xmax);
W2.Draw("Njet30>>hw2");
double w2n = hw2->GetEntries();
hw2->Scale(1./w2n);
hw2->SetLineColor(4);
hw2->SetLineWidth(2);
TH1F *hw3 = new TH1F("hw3","Njet30",Xbin,Xmin,Xmax);
W3.Draw("Njet30>>hw3");
double w3n = hw3->GetEntries();
hw3->Scale(1./w3n);
hw3->SetLineColor(5);
hw3->SetLineWidth(2);
TH1F *hw4 = new TH1F("hw4","Njet30",Xbin,Xmin,Xmax);
W4.Draw("Njet30>>hw4");
double w4n = hw4->GetEntries();
hw4->Scale(1./w4n);
hw4->SetLineColor(6);
hw4->SetLineWidth(2);

/////////////////////////////////////////////////////////////

TCanvas *c1 = new TCanvas("c1"," ",1400,1000);
c1->SetFillColor(10);
c1->SetBorderSize(2);

//c1->Divide(3,3);

//c1->cd(1);
c1->SetGridx();
c1->SetGridy();
//c1->Logy();
hw1->Draw();
ht->Draw("same");
hs->Draw("same");
hw2->Draw("same");
hw3->Draw("same");
hw4->Draw("same");
 TLegend *lg = new TLegend(.65,.65,0.99,1.);
 lg->SetFillColor(0);
 lg->AddEntry(ht,"TTbar","l");
 lg->AddEntry(hw1,"Wjets 100<HT<200","l");
 lg->AddEntry(hw2,"Wjets 200<HT<400","l");
 lg->AddEntry(hw3,"Wjets 400<HT<600","l");
 lg->AddEntry(hw4,"Wjets 600<HT<800","l");
 lg->AddEntry(hs,"Signal (350,270)","l");
 lg->Draw();
 c1->SaveAs("ShapeNjet30_350-270.pdf");

}
