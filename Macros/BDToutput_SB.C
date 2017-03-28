#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "string"
#include "vector"

// Compare normalized BDToutput shapes for S & B
void BDToutput_SB(){

  gStyle->SetOptStat(0);

  TChain SgnFile("bdttree");
  SgnFile.Add("SET3118_ttLO_Zinv/T2DegStop_300_270_bdt.root");

  TChain TopFile("bdttree");
  TopFile.Add("SET3118_ttLO_Zinv/TTJets_bdt.root");

  TChain W1File("bdttree");
  W1File.Add("SET3118_ttLO_Zinv/Wjets_100to200_bdt.root");
  TChain W2File("bdttree");
  W2File.Add("SET3118_ttLO_Zinv/Wjets_200to400_bdt.root");
  TChain W3File("bdttree");
  W3File.Add("SET3118_ttLO_Zinv/Wjets_400to600_bdt.root");
  TChain W4File("bdttree");
  W4File.Add("SET3118_ttLO_Zinv/Wjets_600to800_bdt.root");
  TChain W5File("bdttree");
  W5File.Add("SET3118_ttLO_Zinv/Wjets_800to1200_bdt.root");
  TChain W6File("bdttree");
  W6File.Add("SET3118_ttLO_Zinv/Wjets_1200to2500_bdt.root");
  TChain W7File("bdttree");
  W7File.Add("SET3118_ttLO_Zinv/Wjets_2500toInf_bdt.root");

  TChain Z1File("bdttree");
  Z1File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT100to200_bdt.root");
  TChain Z2File("bdttree");
  Z2File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT200to400_bdt.root");
  TChain Z3File("bdttree");
  Z3File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT400to600_bdt.root");
  TChain Z4File("bdttree");
  Z4File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT600to800_bdt.root");
  TChain Z5File("bdttree");
  Z5File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT800to1200_bdt.root");
  TChain Z6File("bdttree");
  Z6File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT1200to2500_bdt.root");
  TChain Z7File("bdttree");
  Z7File.Add("SET3118_ttLO_Zinv/ZJetsToNuNu_HT2500toInf_bdt.root");


  const double Lumi = 12900.;

  // Lepton selection
  TCut lept = "LepPt<30.";
  //  TCut lept = "LepPt<10000.";
  TCut singlep = lept;

  // Jets
  TCut ISRjet = "Jet1Pt > 110.";
  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";

  // MET
  //  TCut met    = "Met>280.";
  TCut met    = "Met>280.";

  TCut HT = "HT30 > 200.";

  TCut bjet = "JetHBpt > 30.";

  // Preselection
  //  TCut presel = singlep;
  //  TCut presel = singlep+ISRjet;
  //  TCut presel = singlep+ISRjet+dphij1j2;
  TCut presel = singlep+ISRjet+dphij1j2+met+HT;

  TCut Coup = presel;

  const Int_t Xbin=25;
  const Double_t Xmin=-0.50;
  const Double_t Xmax=0.50;

  TH1F *s1bdt = new TH1F("s1bdt","BDToutput",Xbin,Xmin,Xmax);
  s1bdt->Sumw2();
  SgnFile.Draw("BDT>>s1bdt","weight*splitFactor"*Coup);
  s1bdt->Scale(Lumi);
  s1bdt->SetLineColor(4);
  TH1F *tbdt = new TH1F("tbdt","BDToutput",Xbin,Xmin,Xmax);
  tbdt->Sumw2();
  TopFile.Draw("BDT>>tbdt","weight*splitFactor"*Coup);
  tbdt->Scale(Lumi);
  tbdt->SetLineColor(3);
  TH1F *w1bdt = new TH1F("w1bdt","BDToutput",Xbin,Xmin,Xmax);
  w1bdt->Sumw2();
  W1File.Draw("BDT>>w1bdt","weight*splitFactor"*Coup);
  w1bdt->Scale(Lumi);
  TH1F *w2bdt = new TH1F("w2bdt","BDToutput",Xbin,Xmin,Xmax);
  w2bdt->Sumw2();
  W2File.Draw("BDT>>w2bdt","weight*splitFactor"*Coup);
  w2bdt->Scale(Lumi);
  TH1F *w3bdt = new TH1F("w3bdt","BDToutput",Xbin,Xmin,Xmax);
  w3bdt->Sumw2();
  W3File.Draw("BDT>>w3bdt","weight*splitFactor"*Coup);
  w3bdt->Scale(Lumi);
  TH1F *w4bdt = new TH1F("w4bdt","BDToutput",Xbin,Xmin,Xmax);
  w4bdt->Sumw2();
  W4File.Draw("BDT>>w4bdt","weight*splitFactor"*Coup);
  w4bdt->Scale(Lumi);
  TH1F *w5bdt = new TH1F("w5bdt","BDToutput",Xbin,Xmin,Xmax);
  w5bdt->Sumw2();
  W5File.Draw("BDT>>w5bdt","weight*splitFactor"*Coup);
  w5bdt->Scale(Lumi);
  TH1F *w6bdt = new TH1F("w6bdt","BDToutput",Xbin,Xmin,Xmax);
  w6bdt->Sumw2();
  W6File.Draw("BDT>>w6bdt","weight*splitFactor"*Coup);
  w6bdt->Scale(Lumi);
  TH1F *w7bdt = new TH1F("w7bdt","BDToutput",Xbin,Xmin,Xmax);
  w7bdt->Sumw2();
  W7File.Draw("BDT>>w7bdt","weight*splitFactor"*Coup);
  w7bdt->Scale(Lumi);

  TH1F *z1bdt = new TH1F("z1bdt","BDToutput",Xbin,Xmin,Xmax);
  z1bdt->Sumw2();
  Z1File.Draw("BDT>>z1bdt","weight*splitFactor"*Coup);
  z1bdt->Scale(Lumi);
  TH1F *z2bdt = new TH1F("z2bdt","BDToutput",Xbin,Xmin,Xmax);
  z2bdt->Sumw2();
  Z2File.Draw("BDT>>z2bdt","weight*splitFactor"*Coup);
  z2bdt->Scale(Lumi);
  TH1F *z3bdt = new TH1F("z3bdt","BDToutput",Xbin,Xmin,Xmax);
  z3bdt->Sumw2();
  Z3File.Draw("BDT>>z3bdt","weight*splitFactor"*Coup);
  z3bdt->Scale(Lumi);
  TH1F *z4bdt = new TH1F("z4bdt","BDToutput",Xbin,Xmin,Xmax);
  z4bdt->Sumw2();
  Z4File.Draw("BDT>>z4bdt","weight*splitFactor"*Coup);
  z4bdt->Scale(Lumi);
  TH1F *z5bdt = new TH1F("z5bdt","BDToutput",Xbin,Xmin,Xmax);
  z5bdt->Sumw2();
  Z5File.Draw("BDT>>z5bdt","weight*splitFactor"*Coup);
  z5bdt->Scale(Lumi);
  TH1F *z6bdt = new TH1F("z6bdt","BDToutput",Xbin,Xmin,Xmax);
  z6bdt->Sumw2();
  Z6File.Draw("BDT>>z6bdt","weight*splitFactor"*Coup);
  z6bdt->Scale(Lumi);
  TH1F *z7bdt = new TH1F("z7bdt","BDToutput",Xbin,Xmin,Xmax);
  z7bdt->Sumw2();
  Z7File.Draw("BDT>>z7bdt","weight*splitFactor"*Coup);
  z7bdt->Scale(Lumi);

  TH1F *wjetsbdt = static_cast<TH1F*>(w1bdt->Clone("Background"));
  wjetsbdt->Add(w2bdt);
  wjetsbdt->Add(w3bdt);
  wjetsbdt->Add(w4bdt);
  wjetsbdt->Add(w5bdt);
  wjetsbdt->Add(w6bdt);
  wjetsbdt->Add(w7bdt);
  wjetsbdt->SetLineColor(2);

  TH1F *zinvbdt = static_cast<TH1F*>(z1bdt->Clone("Background"));
  zinvbdt->Add(z2bdt);
  zinvbdt->Add(z3bdt);
  zinvbdt->Add(z4bdt);
  zinvbdt->Add(z5bdt);
  zinvbdt->Add(z6bdt);
  zinvbdt->Add(z7bdt);
  zinvbdt->SetLineColor(1);

  tbdt->SetMinimum(1.e-01);
  tbdt->SetMaximum(2.e03);

  TCanvas *c1 = new TCanvas("c1"," ",1100,700);
  c1->SetBorderSize(2);
  c1->SetFillColor(10);
  c1->SetBorderSize(2);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogy();

  tbdt->Draw("l");
  wjetsbdt->Draw("lsame");
  zinvbdt->Draw("lsame");
  s1bdt->Draw("lsame");
  TLegend *lg = new TLegend(.65,.65,0.99,1.);
  lg->SetFillColor(0);
  lg->AddEntry(s1bdt,"Signal (300,270)","l");
  lg->AddEntry(tbdt,"tt","l");
  lg->AddEntry(wjetsbdt,"W+jets","l");
  lg->AddEntry(zinvbdt,"Z \\rightarrow \\nu\\nu","l");
  lg->Draw();

  std::string PlotTitle = "plots/SET3118_ttLO_Zinv/BDToutput_SB.png";
  c1->SaveAs(PlotTitle.c_str());

}
