// Macro calculating the yield of different channels
// The yield is provided with its:
// JES central-value, weighted by all efficiencies (Trigger efficiencies / MC scaling / ...)
// +1sigma(JES) value, weighted by all efficiencies
// -1sigma(JES) value, weighted by all efficiencies
//
// The user has to:
// * provide the luminosity Lumi
// * chose which MC contributions to use
// * modify the section "Define analysis" to modify analysis cuts

#include <iomanip>
#include <fstream>
#include "TMath.h"
#include <iostream>
//#include <cmath>
//#include <ctgmath>
#include "TStyle.h"
#include "TCut.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TEventList.h"
#include "SingleLept_Dev.h"

using namespace std;

void SingleLept_Dev (){
  
  gStyle->SetOptStat(000000);

  //////////////////////////////////////////////////////////////////////////
  // Give luminosity (positive value ==> weighting whereas <=0 value ==> unweighting
  // Luminosity [pb]
  const double Lumi = 12900.;

  //////////////////////////////////////////////////////////////////////////
  // Which MC contributions to use
  bool use_SIGNAL = true;
  bool use_DY1          = true;
  bool use_DY2          = true;
  bool use_DY3          = true;
  bool use_DY4          = true;
  bool use_DY5          = true;
  bool use_DY6          = true;
  bool use_DY7          = true;
  bool use_DY8          = true;
  bool use_W1Jets       = true;
  bool use_W2Jets       = true;
  bool use_W3Jets       = true;
  bool use_W4Jets       = true;
  bool use_W5Jets       = true;
  bool use_W6Jets       = true;
  bool use_W7Jets       = true;
  bool use_Zinv1        = true;
  bool use_Zinv2        = true;
  bool use_Zinv3        = true;
  bool use_Zinv4        = true;
  bool use_Zinv5        = true;
  bool use_Zinv6        = true;
  bool use_Zinv7        = true;
  bool use_WW           = true;
  bool use_WZ           = true;
  bool use_ZZ           = true;
  bool use_SingleToptW  = true;
  bool use_SingleToptC  = true;
  bool use_SingleTopsC  = false;
  bool use_SingleTopBtW = true;
  bool use_SingleTopBtC = true;
  bool use_SingleTopBsC = false;
  bool use_TT           = true;
  bool use_Data         = false;

  //////////////////////////////////////////////////////////////////////////
  // Define analysis

  // MC Triggers
  TCut trgM   = "Event>-999.";
  // Data Triggers
  TCut trgD   = "Event>-999.";

  // Lepton selection
  TCut lept = "LepPt<30.";
  //  TCut lept = "LepPt<10000.";
  //  TCut lept = "(5.<LepPt)&&(LepPt<12.)";
  //  TCut lept = "(12.<LepPt)&&(LepPt<20.)";
  //  TCut lept = "(20.<LepPt)&&(LepPt<30.)";
  TCut singlep = lept;
  TCut etal1 = "abs(LepEta[0])<1.5";
  TCut etal2 = "abs(LepEta[0])<2.4";
  TCut chgl = "LepChg[0]<0.";

  // Jets
  TCut ISRjet = "Jet1Pt>110.";
  TCut ISRjet1 = "Jet1Pt > 100.";
  TCut ISRjet2 = "Jet1Pt > 325.";
  //  TCut dphij1j2 = "DPhiJet1Jet2 < 2.5";
  //  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";
  //  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||((30.<Jet2Pt)&&(Jet2Pt<60.))";
  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";
  //  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)";
  TCut Nhard = "Njet30 < 3";

  // MET
  //  TCut met    = "Met>280.";
  TCut met    = "Met>280.";

  // HT
  TCut HT = "HT30 > 200.";

  // MT
  TCut mta = "mt<60";
  TCut mtb = "(60<mt)&&(mt<95)";
  TCut mtc = "95<mt";

  // btag
  TCut bjet = "JetHBpt > 30.";
  TCut Bloose = "NbLoose30>0";
  TCut NoBloose = "NbLoose30==0";
  TCut Bhard = "NbTight30>0";
  TCut NoBhard = "NbTight30==0";

  // Preselection
  //  TCut presel = trgM;
  //  TCut presel = singlep;
  //  TCut presel = singlep+ISRjet;
  //  TCut presel = singlep+ISRjet+dphij1j2;
  //  TCut presel = singlep+ISRjet+dphij1j2+met;
  TCut presel = singlep+ISRjet+dphij1j2+met+HT;
  //  TCut presel = singlep+ISRjet+dphij1j2+met+bjet;

  // BDT cut
  TCut BDT = "BDT>0.36";

  // Signal regions a la SUS-16-031
  TCut SR1a = met+HT+ISRjet1+Nhard+NoBloose+lept+etal1+chgl+NoBloose+HT+mta;
  TCut SR1b = met+HT+ISRjet1+Nhard+NoBloose+lept+etal1+chgl+NoBloose+HT+mtb;
  TCut SR1c = met+HT+ISRjet1+Nhard+NoBloose+lept+etal1+NoBloose+HT+mtc;
  TCut SR2  = met+ISRjet2+Nhard+Bloose+lept+etal2;

  // 3 parallel JES varied final analysis
  //  TCut Coup   = presel;
  //  TCut Coupjp = presel;
  //  TCut Coupjm = presel;
  //  TCut Coup   = SR2;
  //  TCut Coupjp = SR2;
  //  TCut Coupjm = SR2;
  TCut Coup   = presel+BDT;
  TCut Coupjp = presel+BDT;
  TCut Coupjm = presel+BDT;


  //////////////////////////////////////////////////////////////////////////
  // Contribution / Stat. error / JES plus - minus for different channels
  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(2);
  
  // Data
  TH1F *mdata = new TH1F("mdata","Data",200,-100,100);
  if(use_Data){
    TH1F *mdem = new TH1F("mdem","Data",200,-100,100);
    TChain l1("bdttree");
    //  pb-1
    l1.Add("");
    l1.Add("");
    l1.Draw("nVert>>mdem",trgD+Coup);
    mdata->Add(mdem);
    mdata->Sumw2();
    mdata->SetMarkerStyle(8);

  }

  // Two reference signal points
  double SIGNAL=0, eSIGNAL=0, SIGNALp=0, SIGNALm=0;
  if ( use_SIGNAL ) {
    GetYield(Lumi,"SIGNAL",trgM+Coup,trgM+Coupjp,trgM+Coupjm, SIGNAL, eSIGNAL, SIGNALp, SIGNALm);
  }   

  // 0- TTbar
  double TT=0, eTT=0, TTp=0, TTm=0;
  if ( use_TT ) {
    GetYield(Lumi,"TTbar",trgM+Coup,trgM+Coupjp,trgM+Coupjm, TT, eTT, TTp, TTm);
  }   

  // 2- W+jets
  double W1j=0, eW1j=0, W1jp=0, W1jm=0;
  if ( use_W1Jets ) {
    GetYield(Lumi,"W1Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W1j, eW1j, W1jp, W1jm);
  }
  double W2j=0, eW2j=0, W2jp=0, W2jm=0;
  if ( use_W2Jets ) {
    GetYield(Lumi,"W2Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W2j, eW2j, W2jp, W2jm);
  }
  double W3j=0, eW3j=0, W3jp=0, W3jm=0;
  if ( use_W3Jets ) {
    GetYield(Lumi,"W3Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W3j, eW3j, W3jp, W3jm);
  }
  double W4j=0, eW4j=0, W4jp=0, W4jm=0;
  if ( use_W4Jets ) {
    GetYield(Lumi,"W4Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W4j, eW4j, W4jp, W4jm);
  }
    double W5j=0, eW5j=0, W5jp=0, W5jm=0;
  if ( use_W5Jets ) {
    GetYield(Lumi,"W5Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W5j, eW5j, W5jp, W5jm);
  }
  double W6j=0, eW6j=0, W6jp=0, W6jm=0;
  if ( use_W6Jets ) {
    GetYield(Lumi,"W6Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W6j, eW6j, W6jp, W6jm);
  }
  double W7j=0, eW7j=0, W7jp=0, W7jm=0;
  if ( use_W7Jets ) {
    GetYield(Lumi,"W7Jets",trgM+Coup,trgM+Coupjp,trgM+Coupjm, W7j, eW7j, W7jp, W7jm);
  }
  
  // 4- WW
  double WW=0, eWW=0, WWp=0, WWm=0;
  if ( use_WW ) {
    GetYield(Lumi,"WW",trgM+Coup,trgM+Coupjp,trgM+Coupjm, WW, eWW, WWp, WWm);
  }
  
  // 5- WZ
  double WZ=0, eWZ=0, WZp=0, WZm=0;
  if ( use_WZ ) {
    GetYield(Lumi,"WZ",trgM+Coup,trgM+Coupjp,trgM+Coupjm, WZ, eWZ, WZp, WZm);
  }
  
  // 7- SingleToptW
  double sTw=0, esTw=0, sTwp=0, sTwm=0;
  if ( use_SingleToptW ) {
    GetYield(Lumi,"SingleToptW",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTw, esTw, sTwp, sTwm);
  }
  
  // 8- SingleToptChannel
  double sTtC=0, esTtC=0, sTtCp=0, sTtCm=0;
  if ( use_SingleToptC) {
    GetYield(Lumi,"SingleToptChannel",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTtC, esTtC, sTtCp, sTtCm);
  }
  
  // 9- SingleTopsChannel
  double sTsC=0, esTsC=0, sTsCp=0, sTsCm=0;
  if ( use_SingleTopsC) {
    GetYield(Lumi,"SingleTopsChannel",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTsC, esTsC, sTsCp, sTsCm);
  }
  
   // 7- SingleToptW
  double sTBw=0, esTBw=0, sTBwp=0, sTBwm=0;
  if ( use_SingleTopBtW ) {
    GetYield(Lumi,"SingleTopBtW",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTBw, esTBw, sTBwp, sTBwm);
  }
  
  // 8- SingleToptChannel
  double sTBtC=0, esTBtC=0, sTBtCp=0, sTBtCm=0;
  if ( use_SingleTopBtC) {
    GetYield(Lumi,"SingleTopBtChannel",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTBtC, esTBtC, sTBtCp, sTBtCm);
  }
  
  // 9- SingleTopsChannel
  double sTBsC=0, esTBsC=0, sTBsCp=0, sTBsCm=0;
  if ( use_SingleTopBsC) {
    GetYield(Lumi,"SingleTopBsChannel",trgM+Coup,trgM+Coupjp,trgM+Coupjm, sTBsC, esTBsC, sTBsCp, sTBsCm);
  }
  
  // 1- DY1
  double DY1=0, eDY1=0, DY1p=0, DY1m=0;
  if ( use_DY1 ) {
    GetYield(Lumi,"DY1",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY1, eDY1, DY1p, DY1m);
  }
  
  // 12- DY2
  double DY2=0, eDY2=0, DY2p=0, DY2m=0;
  if ( use_DY2 ) {
    GetYield(Lumi,"DY2",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY2, eDY2, DY2p, DY2m);
  }
  
  // 13- DY3
  double DY3=0, eDY3=0, DY3p=0, DY3m=0;
  if ( use_DY3 ) {
    GetYield(Lumi,"DY3",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY3, eDY3, DY3p, DY3m);
  }
  
  // 14- DY4
  double DY4=0, eDY4=0, DY4p=0, DY4m=0;
  if ( use_DY4 ) {
    GetYield(Lumi,"DY4",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY4, eDY4, DY4p, DY4m);
  }
  
  // 1- DY5
  double DY5=0, eDY5=0, DY5p=0, DY5m=0;
  if ( use_DY5 ) {
    GetYield(Lumi,"DY5",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY5, eDY5, DY5p, DY5m);
  }
  
  // 12- DY6
  double DY6=0, eDY6=0, DY6p=0, DY6m=0;
  if ( use_DY6 ) {
    GetYield(Lumi,"DY6",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY6, eDY6, DY6p, DY6m);
  }
  
  // 13- DY7
  double DY7=0, eDY7=0, DY7p=0, DY7m=0;
  if ( use_DY7 ) {
    GetYield(Lumi,"DY7",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY7, eDY7, DY7p, DY7m);
  }
  
  // 14- DY8
  double DY8=0, eDY8=0, DY8p=0, DY8m=0;
  if ( use_DY8 ) {
    GetYield(Lumi,"DY8",trgM+Coup,trgM+Coupjp,trgM+Coupjm, DY8, eDY8, DY8p, DY8m);
  }

  // 16- Zinv1
  double Zinv1=0, eZinv1=0, Zinv1p=0, Zinv1m=0;
  if ( use_Zinv1 ) {
    GetYield(Lumi,"Zinv1",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv1, eZinv1, Zinv1p, Zinv1m);
  }
  
  // 15- Zinv2
  double Zinv2=0, eZinv2=0, Zinv2p=0, Zinv2m=0;
  if ( use_Zinv2 ) {
    GetYield(Lumi,"Zinv2",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv2, eZinv2, Zinv2p, Zinv2m);
  }
  
  // 17- Zinv3
  double Zinv3=0, eZinv3=0, Zinv3p=0, Zinv3m=0;
  if ( use_Zinv3 ) {
    GetYield(Lumi,"Zinv3",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv3, eZinv3, Zinv3p, Zinv3m);
  }

  // 17- Zinv4
  double Zinv4=0, eZinv4=0, Zinv4p=0, Zinv4m=0;
  if ( use_Zinv4 ) {
    GetYield(Lumi,"Zinv4",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv4, eZinv4, Zinv4p, Zinv4m);
  }

  // 17- Zinv5
  double Zinv5=0, eZinv5=0, Zinv5p=0, Zinv5m=0;
  if ( use_Zinv5 ) {
    GetYield(Lumi,"Zinv5",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv5, eZinv5, Zinv5p, Zinv5m);
  }

  // 17- Zinv6
  double Zinv6=0, eZinv6=0, Zinv6p=0, Zinv6m=0;
  if ( use_Zinv6 ) {
    GetYield(Lumi,"Zinv6",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv6, eZinv6, Zinv6p, Zinv6m);
  }

  // 17- Zinv7
  double Zinv7=0, eZinv7=0, Zinv7p=0, Zinv7m=0;
  if ( use_Zinv7 ) {
    GetYield(Lumi,"Zinv7",trgM+Coup,trgM+Coupjp,trgM+Coupjm, Zinv7, eZinv7, Zinv7p, Zinv7m);
  }

  // 18- ZZ
  double ZZ=0, eZZ=0, ZZp=0, ZZm=0;
  if ( use_ZZ ) {
    GetYield(Lumi,"ZZ",trgM+Coup,trgM+Coupjp,trgM+Coupjm, ZZ, eZZ, ZZp, ZZm);
  }

  //19- Total background
  double bckg  = W1j+W2j+W3j+W4j+W5j+W6j+W7j+ WW+WZ+ZZ +sTw+sTtC+sTsC+sTBw+sTBtC+sTBsC +DY1+DY2+DY3+DY4+DY5+DY6+DY7+DY8 
    +Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7 +TT;
  double bckgp = W1jp+W2jp+W3jp+W4jp+W5jp+W6jp+W7jp+WWp+WZp+ZZp+sTwp+sTtCp+sTsCp+DY1p+DY2p+DY3p+DY4p+DY5p+DY6p+DY7p+DY8p
    +Zinv1p+Zinv2p+Zinv3p+Zinv4p+Zinv5p+Zinv6p+Zinv7p+TTp;
  double bckgm = W1jm+W2jm+W3jm+W4jm+W5jm+W6jm+W7jm+WWm+WZm+ZZm+sTwm+sTtCm+sTsCm+DY1m+DY2m+DY3m+DY4m+DY5m+DY6m+DY7m+DY8m
    +Zinv1p+Zinv2p+Zinv3p+Zinv4p+Zinv5p+Zinv6p+Zinv7m+TTm;
  double ebckg =
    sqrt(pow(eDY1,2)+pow(eW1j,2)+pow(eW2j,2)+pow(eW3j,2)+pow(eW4j,2)+pow(eW5j,2)+pow(eW6j,2)+pow(eW7j,2)+pow(eWW,2)+pow(eWZ,2)+
	 pow(esTw,2)+pow(esTtC,2)+pow(esTsC,2)+pow(esTBw,2)+pow(esTBtC,2)+pow(esTBsC,2)+pow(eDY2,2)+pow(eDY3,2)+pow(eDY4,2)+
	 pow(eZinv1,2)+pow(eZinv2,2)+pow(eZinv3,2)+pow(eZinv4,2)+pow(eZinv5,2)+pow(eZinv6,2)+pow(eZinv7,2)+pow(eTT,2));

  cout << " " <<  endl;
  cout << "SIGNAL     : " << SIGNAL << " +- " << eSIGNAL << " (stat) +" << SIGNALp-SIGNAL << " -" << SIGNAL-SIGNALm << " (jes)" << "\n";
  cout << "BACKGROUND : " << bckg << " +- " << ebckg << " (stat) +" << bckgp-bckg << " -" << bckg-bckgm << " (jes)"  << "\n";
  if (use_Data) {
    cout << "Data : " << mdata->GetEntries() << "\n";
  }
  cout << "---------------------------------------------------------" << "\n"; 

  if ( use_W1Jets ) {
    cout << "W1Jets  : " << W1j << " +- " << eW1j << " (stat) +" << W1jp-W1j << " -" << W1j-W1jm << " (jes)" << "\n";
  }
  if ( use_W2Jets ) {
    cout << "W2Jets  : " << W2j << " +- " << eW2j << " (stat) +" << W2jp-W2j << " -" << W2j-W2jm << " (jes)" << "\n";
  }
  if ( use_W3Jets ) {
    cout << "W3Jets  : " << W3j << " +- " << eW3j << " (stat) +" << W3jp-W3j << " -" << W3j-W3jm << " (jes)" << "\n";
  }
  if ( use_W4Jets ) {
    cout << "W4Jets  : " << W4j << " +- " << eW4j << " (stat) +" << W4jp-W4j << " -" << W4j-W4jm << " (jes)" << "\n";
  }
  if ( use_W5Jets ) {
    cout << "W5Jets  : " << W5j << " +- " << eW5j << " (stat) +" << W5jp-W5j << " -" << W5j-W5jm << " (jes)" << "\n";
  }
  if ( use_W6Jets ) {
    cout << "W6Jets  : " << W6j << " +- " << eW6j << " (stat) +" << W6jp-W6j << " -" << W6j-W6jm << " (jes)" << "\n";
  }
  if ( use_W7Jets ) {
    cout << "W7Jets  : " << W7j << " +- " << eW7j << " (stat) +" << W7jp-W7j << " -" << W7j-W7jm << " (jes)" << "\n";
  }
  if ( use_TT ) {
    cout << "TTbar   : " << TT << " +- " << eTT << " (stat) +" << TTp-TT << " -" << TT-TTm << " (jes)" << "\n";
  }
  if ( use_Zinv1 ) {
    cout << "Zinv1   : " << Zinv1 << " +- " << eZinv1 << " (stat) +" << Zinv1p-Zinv1 << " -" << Zinv1-Zinv1m << " (jes)" << "\n";
  }
  if ( use_Zinv2 ) {
    cout << "Zinv2   : " << Zinv2 << " +- " << eZinv2 << " (stat) +" << Zinv2p-Zinv2 << " -" << Zinv2-Zinv2m << " (jes)" << "\n";
  }
  if ( use_Zinv3 ) {
    cout << "Zinv3   : " << Zinv3 << " +- " << eZinv3 << " (stat) +" << Zinv3p-Zinv3 << " -" << Zinv3-Zinv3m << " (jes)" << "\n";
  }
  if ( use_Zinv4 ) {
    cout << "Zinv4   : " << Zinv4 << " +- " << eZinv4 << " (stat) +" << Zinv4p-Zinv4 << " -" << Zinv4-Zinv4m << " (jes)" << "\n";
  }
  if ( use_Zinv5 ) {
    cout << "Zinv5   : " << Zinv5 << " +- " << eZinv5 << " (stat) +" << Zinv5p-Zinv5 << " -" << Zinv5-Zinv5m << " (jes)" << "\n";
  }
  if ( use_Zinv6 ) {
    cout << "Zinv6   : " << Zinv6 << " +- " << eZinv6 << " (stat) +" << Zinv6p-Zinv6 << " -" << Zinv6-Zinv6m << " (jes)" << "\n";
  }
  if ( use_Zinv7 ) {
    cout << "Zinv7   : " << Zinv7 << " +- " << eZinv7 << " (stat) +" << Zinv7p-Zinv7 << " -" << Zinv7-Zinv7m << " (jes)" << "\n";
  }
  if ( use_WW || use_WZ || use_ZZ ) {
    cout << "VV      : " << WW+WZ+ZZ << " +- " << sqrt(pow(eWW,2)+pow(eWZ,2)+pow(eZZ,2)) << " (stat) +" << (WWp-WW)+(WZp-WZ)+(ZZp-ZZ) << " -" << (WW-WWm)+(WZ-WZm)+(ZZ-ZZm) << " (jes)" << "\n";
  }
  if ( use_SingleToptW || use_SingleToptC || use_SingleTopsC || use_SingleTopBtW || use_SingleTopBtC || use_SingleTopBsC ) {
    cout << "S Top   : " << sTw+sTtC+sTsC+sTBw+sTBtC+sTBsC << " +- " 
	 << sqrt( pow(esTw,2)+pow(esTtC,2)+pow(esTsC,2)+pow(esTBw,2)+pow(esTBtC,2)+pow(esTBsC,2) ) << " (stat) +" 
	 << (sTwp-sTw)+(sTtCp-sTtC)+(sTsCp-sTsC)+(sTBwp-sTBw)+(sTBtCp-sTBtC)+(sTBsCp-sTBsC) << " -" 
	 << (sTw-sTwm)+(sTtC-sTtCm)+(sTsC-sTsCm)+(sTBw-sTBwm)+(sTBtC-sTBtCm)+(sTBsC-sTBsCm) << " (jes)" << "\n";
  }
  if ( use_DY1 || use_DY2 || use_DY3 || use_DY4 || use_DY5 || use_DY6 || use_DY7 || use_DY8 ) {
    cout << "DY      : " << DY1+DY2+DY3+DY4+DY5+DY6+DY7+DY8 << " +- " 
	 << sqrt( pow(eDY1,2)+pow(eDY2,2)+pow(eDY3,2)+pow(eDY4,2)+pow(eDY5,2)+pow(eDY6,2)+pow(eDY7,2)+pow(eDY8,2) ) << " (stat) +" 
	 << (DY1p-DY1)+(DY2p-DY2)+(DY3p-DY3)+(DY4p-DY4) << " -" 
	 << (DY1-DY1m)+(DY2-DY2m)+(DY3-DY3m)+(DY4-DY4m) << " (jes)" << "\n";
  }

  cout << "---------------------------------------------------------" << "\n"; 
  
// cout << "\n"; 

// Choose which signal should be looked upon
  double Sgn = SIGNAL;
// FOM
  double f = 0.2;
  // http://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf
  double SysB = f * bckg;
  double FOMn1 = (Sgn+bckg)*(bckg+pow(SysB,2));
  double FOMd1 = pow(bckg,2) + (Sgn+bckg)*pow(SysB,2);
  double FOM1 = (Sgn+bckg)*log(FOMn1/FOMd1);
  double FOMn2 = pow(SysB,2)*Sgn;
  double FOMd2 = bckg*(bckg+pow(SysB,2));
  double FOM2 = pow(bckg/SysB,2)*log(1.+(FOMn2/FOMd2));
  double FOM = sqrt(2.*(FOM1 - FOM2));

  double sS2 = Sgn;
  double sB2 = pow(ebckg,2);
  double f2 = pow(f,2);
  double lnu = 1. + (f2*(Sgn+bckg));
  double lde = (1. + (bckg*f2))*(Sgn+bckg);
  double lognum = log(lnu/lde);
  double numbr = bckg*(1.+(bckg*f2))*(sS2+sB2)*lognum;
  double numb =  (2.*Sgn*sB2) + numbr;
  double nnum = (pow(Sgn,2)*sB2) + ( bckg*(1.+(bckg*f2))*lognum*numb );
  double lnd = (1.+(bckg*f2))*(bckg+Sgn);
  double ldd = bckg*(1. + (f2*(bckg+Sgn)));
  double nden1 = f2*(Sgn+bckg)*log(lnd/ldd);
  double nden2 = log( 1. + (Sgn*f2/(1.+(bckg*f2))) );
  double nden = nden1 - nden2;
  double num = f * sqrt(nnum/nden);
  double den = sqrt(2.)*bckg*(1.+(bckg*f2));
  double dFOM = num/den;

  cout << "FOM = " << FOM << "\n";
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n"; 

  // Create file with yields
  string fileName =  "CutDevTable";
  ofstream yieldFile;

  string FOMname = "FOM";
  string Sname   = "Signal";
  string Bname   = "Background";
  string W1name  = "Wjets [100,200]";
  string W2name  = "Wjets [200,400]";
  string W3name  = "Wjets [400,600]";
  string W4name  = "Wjets [600,800]";
  string W5name  = "Wjets [800,1200]";
  string W6name  = "Wjets [1200,2500]";
  string W7name  = "Wjets [2500,Inf]";
  string WJname  = "Wjets";
  string TTname  = "tt";
  string Z1name  = "Zinv [100,200]";
  string Z2name  = "Zinv [200,400]";
  string Z3name  = "Zinv [400,600]";
  string Z4name  = "Zinv [600,800]";
  string Z5name  = "Zinv [800,1200]";
  string Z6name  = "Zinv [1200,2500]";
  string Z7name  = "Zinv [2500,Inf]";
  string ZIname  = "Zinv";
  string WWname  = "WW";
  string WZname  = "WZ";
  string ZZname  = "ZZ";
  string VVname  = "VV";
  string STname  = "Single tops";
  string DYname  = "DY";

  yieldFile.open((fileName + ".tex").c_str());
  yieldFile << "\\documentclass{article}" << std::endl;
  yieldFile << "\\usepackage[utf8]{inputenc}" << std::endl;
  yieldFile << "\\usepackage{cancel}" << std::endl;
  yieldFile << "\\begin{document}" << std::endl;
  yieldFile << "\\begin{table}[!h]" << std::endl;
  yieldFile << "\\centering" << std::endl;
  string col = "ll";
  yieldFile << "\\begin{tabular}{" << col << "}" << std::endl;
  yieldFile << "\\hline" << std::endl;
  yieldFile << "\\hline" << std::endl;

  yieldFile << FOMname           << " & " << std::setprecision(3) << FOM  << "\\\\" << std::endl;
  yieldFile << Sname             << " & " << std::setprecision(3) << SIGNAL << " $\\pm$ " << eSIGNAL << "\\\\" << std::endl;
  yieldFile << Bname             << " & " << std::setprecision(3) << bckg << " $\\pm$ " << ebckg << "\\\\" << std::endl;
  yieldFile << "\\hline" << std::endl;
  yieldFile << WJname            << "&" << std::setprecision(3) << W1j+W2j+W3j+W4j+W5j+W6j+W7j << " $\\pm$ " << std::setprecision(2) <<
    sqrt(pow(eW1j,2)+pow(eW2j,2)+pow(eW3j,2)+pow(eW4j,2)+pow(eW5j,2)+pow(eW6j,2)+pow(eW7j,2)) << "\\\\" << std::endl;
  /*
  yieldFile << W1name            << "&" << std::setprecision(2) << W1j << "\\\\" << std::endl;
  yieldFile << W2name            << "&" << std::setprecision(2) << W2j << "\\\\" << std::endl;
  yieldFile << W3name            << "&" << std::setprecision(2) << W3j << "\\\\" << std::endl;
  yieldFile << W4name            << "&" << std::setprecision(2) << W4j << "\\\\" << std::endl;
  yieldFile << W5name            << "&" << std::setprecision(2) << W5j << "\\\\" << std::endl;
  yieldFile << W6name            << "&" << std::setprecision(2) << W6j << "\\\\" << std::endl;
  yieldFile << W7name            << "&" << std::setprecision(2) << W7j << "\\\\" << std::endl;
  */
  yieldFile << TTname            << "&" << std::setprecision(3) << TT  << "\\\\" << std::endl;
  yieldFile << ZIname            << "&" << std::setprecision(3) << Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7  << " $\\pm$ " << std::setprecision(2) <<
    sqrt(pow(eZinv1,2)+pow(eZinv2,2)+pow(eZinv3,2)+pow(eZinv4,2)+pow(eZinv5,2)+pow(eZinv6,2)+pow(eZinv7,2)) << "\\\\" << std::endl;
  /*
  yieldFile << Z1name            << "&" << std::setprecision(2) << Zinv1  << " $\\pm$ " << eZinv1 << "\\\\" << std::endl;
  yieldFile << Z2name            << "&" << std::setprecision(2) << Zinv2  << " $\\pm$ " << eZinv2 << "\\\\" << std::endl;
  yieldFile << Z3name            << "&" << std::setprecision(2) << Zinv3  << " $\\pm$ " << eZinv3 << "\\\\" << std::endl;
  yieldFile << Z4name            << "&" << std::setprecision(2) << Zinv4  << " $\\pm$ " << eZinv4 << "\\\\" << std::endl;
  yieldFile << Z5name            << "&" << std::setprecision(2) << Zinv5  << " $\\pm$ " << eZinv5 << "\\\\" << std::endl;
  yieldFile << Z6name            << "&" << std::setprecision(2) << Zinv6  << " $\\pm$ " << eZinv6 << "\\\\" << std::endl;
  yieldFile << Z7name            << "&" << std::setprecision(2) << Zinv7  << " $\\pm$ " << eZinv7 << "\\\\" << std::endl;
  */
  yieldFile << VVname            << "&" << std::setprecision(3) << WW+WZ+ZZ  << " $\\pm$ " << std::setprecision(2) << 
    sqrt(pow(eWW,2)+pow(eWZ,2)+pow(eZZ,2)) << "\\\\" << std::endl;
  yieldFile << STname            << "&" << std::setprecision(3) << sTw+sTtC+sTsC+sTBw+sTBtC+sTBsC << " $\\pm$ " << std::setprecision(2) <<
    sqrt( pow(esTw,2)+pow(esTtC,2)+pow(esTsC,2)+pow(esTBw,2)+pow(esTBtC,2)+pow(esTBsC,2) ) << "\\\\" << std::endl;
  yieldFile << DYname            << "&" << std::setprecision(3) << DY1+DY2+DY3+DY4+DY5+DY6+DY7 << " $\\pm$ " << std::setprecision(2) << 
    sqrt( pow(eDY1,2)+pow(eDY2,2)+pow(eDY3,2)+pow(eDY4,2)+pow(eDY5,2)+pow(eDY6,2)+pow(eDY7,2)+pow(eDY8,2) ) << "\\\\" << std::endl;

  yieldFile << "\\hline" << std::endl;
  yieldFile << "\\hline" << std::endl;
  yieldFile << "\\end{tabular}" << std::endl;
  //  yieldFile << "\\caption{FOM}" << std::endl;
  yieldFile << "\\end{table}" << std::endl;
  yieldFile << "\\end{document}" << std::endl;

}


void GetYield(
	      double lumi,     // Input: Luminosity
	      string pc,       // Input: Name of the channel
	      TCut Coupure,    // Input: Analysis cut
	      TCut Coupurejp,  // Input: +1sigma(JES) analysis cut
	      TCut Coupurejm,  // Input: -1sigma(JES) analysis cut
	      double& Yield,   // Output: Weighted yield of the channel
	      double& eYield,  // Output: Error on the wighted yield
	      double& YieldP,  // Output: Weighted yield of the channel +1sigma(JES)
	      double& YieldM  // Output: Weighted yield of the channel -1sigma(JES)
	      ){ 

  // Which file to read
  TChain ProcFile("bdttree");

  // Directory of roottuples
  string BaseDir = "SET3112_ttLO_Zinv/";

  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3112_DM10/T2DegStop_300_290_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3112_DM20/T2DegStop_300_280_bdt.root");  }
  if (pc=="SIGNAL"){    ProcFile.Add("SET3112_ttLO_Zinv/T2DegStop_300_270_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3113_DM40/T2DegStop_300_260_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3113_DM50/T2DegStop_300_250_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3113_DM60/T2DegStop_300_240_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3113_DM70/T2DegStop_300_230_bdt.root");  }
  //  if (pc=="SIGNAL"){    ProcFile.Add("SET3113_DM80/T2DegStop_350_270_bdt.root");  }
  if (pc=="TTbar"){    ProcFile.Add("SET3112_ttLO_Zinv/TTJets_LO_bdt.root");  }
  if (pc=="DY1"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M50_HT100to200_bdt.root");  }
  if (pc=="DY2"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M50_HT200to400_bdt.root");  }
  if (pc=="DY3"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M50_HT400to600_bdt.root");  }
  if (pc=="DY4"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M50_HT600toInf_bdt.root");  }
  if (pc=="DY5"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M5to50_HT100to200_bdt.root");  }
  if (pc=="DY6"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M5to50_HT200to400_bdt.root");  }
  if (pc=="DY7"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M5to50_HT400to600_bdt.root");  }
  if (pc=="DY8"){    ProcFile.Add("SET3112_ttLO_Zinv/DYJetsToLL_M5to50_HT600toInf_bdt.root");  }
  if (pc=="W1Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_100to200_bdt.root");  }
  if (pc=="W2Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_200to400_bdt.root");  }
  if (pc=="W3Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_400to600_bdt.root");  }
  if (pc=="W4Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_600to800_bdt.root");  }
  if (pc=="W5Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_800to1200_bdt.root");  }
  if (pc=="W6Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_1200to2500_bdt.root");  }
  if (pc=="W7Jets"){    ProcFile.Add("SET3112_ttLO_Zinv/Wjets_2500toInf_bdt.root");  }
  if (pc=="Zinv1"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT100to200_bdt.root");  }
  if (pc=="Zinv2"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT200to400_bdt.root");  }
  if (pc=="Zinv3"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT400to600_bdt.root");  }
  if (pc=="Zinv4"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT600to800_bdt.root");  }
  if (pc=="Zinv5"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT800to1200_bdt.root");  }
  if (pc=="Zinv6"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT1200to2500_bdt.root");  }
  if (pc=="Zinv7"){    ProcFile.Add("SET3112_ttLO_Zinv/ZJetsToNuNu_HT2500toInf_bdt.root");  }
  if (pc=="WW"){    ProcFile.Add("SET3112_ttLO_Zinv/WW_bdt.root");  }
  if (pc=="WZ"){    ProcFile.Add("SET3112_ttLO_Zinv/WZ_bdt.root");  }
  if (pc=="ZZ"){    ProcFile.Add("SET3112_ttLO_Zinv/ZZ_bdt.root");  }
  if (pc=="SingleToptW"){    ProcFile.Add("SET3112_ttLO_Zinv/T_tWch_bdt.root");  }
  if (pc=="SingleToptChannel"){    ProcFile.Add("SET3112_ttLO_Zinv/T_tch_bdt.root");  }
  if (pc=="SingleTopsChannel"){    ProcFile.Add("");  }
  if (pc=="SingleTopBtW"){    ProcFile.Add("SET3112_ttLO_Zinv/TBar_tWch_bdt.root");  }
  if (pc=="SingleTopBtChannel"){    ProcFile.Add("SET3112_ttLO_Zinv/TBar_tch_bdt.root");  }
  if (pc=="SingleTopBsChannel"){    ProcFile.Add("");  }

  /*
  // Get the numerator
  TH1D *nbcut = new TH1D("nbcut","nbcut",200,-100,100);
  TH1D *nbcutjp = new TH1D("nbcutjp","nbcutjp",200,-100,100);
  TH1D *nbcutjm = new TH1D("nbcutjm","nbcutjm",200,-100,100);
  ProcFile.Draw("Njet>>nbcut",Coupure);    // Pick-up a variable which is always filled
  ProcFile.Draw("Njet>>nbcutjp",Coupurejp);
  ProcFile.Draw("Njet>>nbcutjm",Coupurejm);
  double NumWeight = nbcut->GetEntries();
  double NumWeightP = nbcutjp->GetEntries();
  double NumWeightM = nbcutjm->GetEntries();

  // Get the weight
  float xsec, fileff;
  ProcFile.SetBranchAddress("XS",&xsec);
  ProcFile.GetEntry(0);
  ProcFile.SetBranchAddress("filterEfficiency",&fileff);
  ProcFile.GetEntry(0);
  xsec *= fileff;
  // Get the denominator
  TH1D *nbt = new TH1D("nbt","nbt",10,0.,10000000.);
  ProcFile.Draw("Nevt>>nbt");
  double *Den = ProcFile.GetV1();
  // Attention: Only valid for splitted MC samples
  NumWeight *= 2.;
  NumWeightP *= 2.;
  NumWeightM *= 2.;
  // Calculate yields (weighted/unweighted)
  Yield  = lumi*xsec*NumWeight/Den[0];
  eYield = lumi*xsec*uneff(NumWeight,Den[0]);
  YieldP = lumi*xsec*NumWeightP/Den[0];
  YieldM = lumi*xsec*NumWeightM/Den[0];
  */

  // Get the weighted numerator, where division by denominator has been done
  TH1D *nbcut = new TH1D("nbcut","nbcut",200,-100,100);
  TH1D *nbcutjp = new TH1D("nbcutjp","nbcutjp",200,-100,100);
  TH1D *nbcutjm = new TH1D("nbcutjm","nbcutjm",200,-100,100);
  ProcFile.Draw("Njet>>nbcut","weight*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutjp","weight*splitFactor"*Coupurejp);
  ProcFile.Draw("Njet>>nbcutjm","weight*splitFactor"*Coupurejm);

  double NumWeight = 0.;
  double NumWeightP = 0.;
  double NumWeightM = 0.;
  for(int i=0; i<=(nbcut->GetNbinsX()+1); i++){ NumWeight += nbcut->GetBinContent(i); }
  for(int i=0; i<=(nbcutjp->GetNbinsX()+1); i++){ NumWeightP += nbcutjp->GetBinContent(i); }
  for(int i=0; i<=(nbcutjm->GetNbinsX()+1); i++){ NumWeightM += nbcutjm->GetBinContent(i); }

  // Get the numerator
  TH1D *njcut = new TH1D("njcut","njcut",200,-100,100);
  ProcFile.Draw("Njet>>njcut",Coupure);
  double Num = njcut->GetEntries();
  // Get the denominator
  TH1D *nbt = new TH1D("nbt","nbt",10,0.,10000000.);
  ProcFile.Draw("Nevt>>nbt");
  double *Den = ProcFile.GetV1();

  // Calculate yields 
  Yield  = lumi*NumWeight;
  eYield = lumi*uneff(Num,Den[0]);
  YieldP = lumi*NumWeightP;
  YieldM = lumi*NumWeightM;

  // Delete histograms
  delete nbt;
  delete nbcut;
  delete nbcutjp;
  delete nbcutjm;
  delete njcut;

  NumWeight = 0;
  NumWeightP = 0;
  NumWeightM = 0;

}


double eff(double a, double b){ 
  double af = double(a);
  double bf = double(b);   
  double effi = af/bf;
  return effi;
}
double uneff(double a, double b){
  double af = double(a);
  double bf = double(b);
  double r = af/bf;
  double unc = sqrt(af + (r*r*bf) )/bf;
  return unc;
}

double Max(double a, double b){
  if (a>b) return a;
  else return b;
}

double AbsV(double a){
  return sqrt(pow(a,2));
}
