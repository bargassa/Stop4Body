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
#include "TStyle.h"
#include "TCut.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TEventList.h"
#include "SingleLept_CutFlow.h"

using namespace std;

void SingleLept_CutFlow (){
  
  gStyle->SetOptStat(000000);

  //////////////////////////////////////////////////////////////////////////
  // Give luminosity (positive value ==> weighting whereas <=0 value ==> unweighting
  // Luminosity
  const double Lumi = 35866.0;

  //////////////////////////////////////////////////////////////////////////
  // Define analysis

  // Data Triggers
  TCut trgD   = "Event>-999.";

  // Lepton selection
  TCut Promp = "isPrompt==1";
  TCut singlep = "(isTight==1)";
  //  TCut lept = "LepPt<30.";
  TCut lept = "LepPt<10000.";
  TCut etal1 = "abs(LepEta[0])<1.5";
  TCut etal2 = "abs(LepEta[0])<2.4";
  TCut chgl = "LepChg[0]<0.";

  // Jets
  TCut ISRjet = "Jet1Pt > 110.";
  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";

  // MET
  TCut met    = "Met>280.";

  // HT
  TCut HT = "HT > 200.";

  // MT
  TCut mta = "mt<60";
  TCut mtb = "(60<mt)&&(mt<95)";
  TCut mtc = "95<mt";

  // btag
  TCut bjet = "JetHBpt > 30.";
  TCut Bloose = "NbLoose>0";
  TCut NoBloose = "NbLoose==0";
  TCut Bhard = "NbTight>0";
  TCut NoBhard = "NbTight==0";

  // Preselection
  TCut presel = singlep+lept+ISRjet+dphij1j2+met+HT;

  // BDT cut
  TCut BDT = "BDT>0.44";


  //////////////////////////////////////////////////////////////////////////
  
  // Prepare the .tex settings
  ofstream     tablesFile    ("SingleptCutFlow.tex");
  tablesFile.setf(ios::fixed);
  tablesFile.precision(1);
  tablesFile << "\\documentclass[amsmath,amssymb]{revtex4}" << endl;
  tablesFile << "\\oddsidemargin -0.5in" << endl;
  tablesFile << "\\usepackage{longtable}" << endl;
  tablesFile << "\\usepackage{color}" << endl;
  tablesFile << "\\begin{document}" << endl;
  tablesFile << "\\clearpage" << endl << endl << endl;
  tablesFile << "\\begin{table}[!htbp]" << endl;
  tablesFile << "\\begin{center}" << endl;
  tablesFile << "\\begin{tabular}{|l|cccc|c|c|}" << endl;
  tablesFile << "\\hline" << endl;
  tablesFile << "\\hline" << endl;
  tablesFile << "Selection      &      &                  &                            &        & Total      &   Signal  \\\\" << endl;
  tablesFile << "level & W+jets & $t\\overline{t}$ & $Z \\rightarrow \\nu \\nu$ & Other         & Background &  \\\\" << endl;
  tablesFile << "\\hline" << endl;

  // Loop over cuts
  TString CutName[10];
  int nbcutmax  = 2;
  
  for(int ic=0; ic<nbcutmax; ++ic) {
    
    // Cut specifics
    TCut Coup,Coupjp,Coupjm,CoupB,CoupBjp,CoupBjm;
    if (ic==0){
      Coup   = presel;
      Coupjp = presel;
      Coupjm = presel;
      CoupB   = presel+Promp;
      CoupBjp = presel+Promp;
      CoupBjm = presel+Promp;
      CutName[ic] = "Preselection";
    }
    if (ic==1){
      Coup   = presel+BDT;
      Coupjp = presel+BDT;
      Coupjm = presel+BDT;
      CoupB   = presel+Promp+BDT;
      CoupBjp = presel+Promp+BDT;
      CoupBjm = presel+Promp+BDT;
      CutName[ic] = "$BDT \\geq 0.44$";
    }

    cout << " " <<  endl;
    cout << " " <<  endl;
    cout << "selection " << CutName[ic] << endl;
    //////////////////////////////////////////////////////////////////////////
    // Contribution / Stat. error / JES plus - minus for different channels
    cout.setf(ios::floatfield,ios::fixed);
    cout<<setprecision(2);

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
  bool use_DY9          = true;
  bool use_DY10         = true;
  bool use_DY11         = true;
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
  bool use_TTgj         = true;
  bool use_TTwl         = true;
  bool use_TTwqq        = true;
  bool use_TTzl         = true;
  bool use_TTzl1m10     = true;
  bool use_TTzqq        = true;
  bool use_Data         = false;


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


  // Signal point
  double SIGNAL=0, eSIGNAL=0, SIGNALp=0, SIGNALm=0;
  if ( use_SIGNAL ) {    GetYield(Lumi,"SIGNAL",Coup,Coupjp,Coupjm, SIGNAL, eSIGNAL, SIGNALp, SIGNALm);  }   

  //  TTbar
  double TT=0, eTT=0, TTp=0, TTm=0;
  if ( use_TT ) {    GetYield(Lumi,"TTbar",CoupB,CoupBjp,CoupBjm, TT, eTT, TTp, TTm);  }   
  // W+jets
  double W1j=0, eW1j=0, W1jp=0, W1jm=0;
  if ( use_W1Jets ) {    GetYield(Lumi,"W1Jets",CoupB,CoupBjp,CoupBjm, W1j, eW1j, W1jp, W1jm);  }
  double W2j=0, eW2j=0, W2jp=0, W2jm=0;
  if ( use_W2Jets ) {    GetYield(Lumi,"W2Jets",CoupB,CoupBjp,CoupBjm, W2j, eW2j, W2jp, W2jm);  }
  double W3j=0, eW3j=0, W3jp=0, W3jm=0;
  if ( use_W3Jets ) {    GetYield(Lumi,"W3Jets",CoupB,CoupBjp,CoupBjm, W3j, eW3j, W3jp, W3jm);  }
  double W4j=0, eW4j=0, W4jp=0, W4jm=0;
  if ( use_W4Jets ) {    GetYield(Lumi,"W4Jets",CoupB,CoupBjp,CoupBjm, W4j, eW4j, W4jp, W4jm);  }
  double W5j=0, eW5j=0, W5jp=0, W5jm=0;
  if ( use_W5Jets ) {    GetYield(Lumi,"W5Jets",CoupB,CoupBjp,CoupBjm, W5j, eW5j, W5jp, W5jm);  }
  double W6j=0, eW6j=0, W6jp=0, W6jm=0;
  if ( use_W6Jets ) {    GetYield(Lumi,"W6Jets",CoupB,CoupBjp,CoupBjm, W6j, eW6j, W6jp, W6jm);  }
  double W7j=0, eW7j=0, W7jp=0, W7jm=0;
  if ( use_W7Jets ) {    GetYield(Lumi,"W7Jets",CoupB,CoupBjp,CoupBjm, W7j, eW7j, W7jp, W7jm);  }
  // WW
  double WW=0, eWW=0, WWp=0, WWm=0;
  if ( use_WW ) {    GetYield(Lumi,"WW",CoupB,CoupBjp,CoupBjm, WW, eWW, WWp, WWm);  }
  // WZ
  double WZ=0, eWZ=0, WZp=0, WZm=0;
  if ( use_WZ ) {    GetYield(Lumi,"WZ",CoupB,CoupBjp,CoupBjm, WZ, eWZ, WZp, WZm);  }
  // SingleToptW
  double sTw=0, esTw=0, sTwp=0, sTwm=0;
  if ( use_SingleToptW ) {    GetYield(Lumi,"SingleToptW",CoupB,CoupBjp,CoupBjm, sTw, esTw, sTwp, sTwm);  }
  // SingleToptChannel
  double sTtC=0, esTtC=0, sTtCp=0, sTtCm=0;
  if ( use_SingleToptC) {    GetYield(Lumi,"SingleToptChannel",CoupB,CoupBjp,CoupBjm, sTtC, esTtC, sTtCp, sTtCm);  }
  // SingleTopsChannel
  double sTsC=0, esTsC=0, sTsCp=0, sTsCm=0;
  if ( use_SingleTopsC) {    GetYield(Lumi,"SingleTopsChannel",CoupB,CoupBjp,CoupBjm, sTsC, esTsC, sTsCp, sTsCm);  }
  // SingleToptW
  double sTBw=0, esTBw=0, sTBwp=0, sTBwm=0;
  if ( use_SingleTopBtW ) {    GetYield(Lumi,"SingleTopBtW",CoupB,CoupBjp,CoupBjm, sTBw, esTBw, sTBwp, sTBwm);  }
  // SingleToptChannel
  double sTBtC=0, esTBtC=0, sTBtCp=0, sTBtCm=0;
  if ( use_SingleTopBtC) {    GetYield(Lumi,"SingleTopBtChannel",CoupB,CoupBjp,CoupBjm, sTBtC, esTBtC, sTBtCp, sTBtCm);  }
  // SingleTopsChannel
  double sTBsC=0, esTBsC=0, sTBsCp=0, sTBsCm=0;
  if ( use_SingleTopBsC) {    GetYield(Lumi,"SingleTopBsChannel",CoupB,CoupBjp,CoupBjm, sTBsC, esTBsC, sTBsCp, sTBsCm);  }
  // DY1
  double DY1=0, eDY1=0, DY1p=0, DY1m=0;
  if ( use_DY1 ) {    GetYield(Lumi,"DY1",CoupB,CoupBjp,CoupBjm, DY1, eDY1, DY1p, DY1m);  }
  // DY2
  double DY2=0, eDY2=0, DY2p=0, DY2m=0;
  if ( use_DY2 ) {    GetYield(Lumi,"DY2",CoupB,CoupBjp,CoupBjm, DY2, eDY2, DY2p, DY2m);  }
  // DY3
  double DY3=0, eDY3=0, DY3p=0, DY3m=0;
  if ( use_DY3 ) {    GetYield(Lumi,"DY3",CoupB,CoupBjp,CoupBjm, DY3, eDY3, DY3p, DY3m);  }
  // DY4
  double DY4=0, eDY4=0, DY4p=0, DY4m=0;
  if ( use_DY4 ) {    GetYield(Lumi,"DY4",CoupB,CoupBjp,CoupBjm, DY4, eDY4, DY4p, DY4m);  }
  // DY5
  double DY5=0, eDY5=0, DY5p=0, DY5m=0;
  if ( use_DY5 ) {    GetYield(Lumi,"DY5",CoupB,CoupBjp,CoupBjm, DY5, eDY5, DY5p, DY5m);  }
  // DY6
  double DY6=0, eDY6=0, DY6p=0, DY6m=0;
  if ( use_DY6 ) {    GetYield(Lumi,"DY6",CoupB,CoupBjp,CoupBjm, DY6, eDY6, DY6p, DY6m);  }
  // DY7
  double DY7=0, eDY7=0, DY7p=0, DY7m=0;
  if ( use_DY7 ) {    GetYield(Lumi,"DY7",CoupB,CoupBjp,CoupBjm, DY7, eDY7, DY7p, DY7m);  }
  // DY8
  double DY8=0, eDY8=0, DY8p=0, DY8m=0;
  if ( use_DY8 ) {    GetYield(Lumi,"DY8",CoupB,CoupBjp,CoupBjm, DY8, eDY8, DY8p, DY8m);  }
  // DY9
  double DY9=0, eDY9=0, DY9p=0, DY9m=0;
  if ( use_DY9 ) {    GetYield(Lumi,"DY9",CoupB,CoupBjp,CoupBjm, DY9, eDY9, DY9p, DY9m);  }
  // DY10
  double DY10=0, eDY10=0, DY10p=0, DY10m=0;
  if ( use_DY10 ) {    GetYield(Lumi,"DY10",CoupB,CoupBjp,CoupBjm, DY10, eDY10, DY10p, DY10m);  }
  // DY11
  double DY11=0, eDY11=0, DY11p=0, DY11m=0;
  if ( use_DY11 ) {    GetYield(Lumi,"DY11",CoupB,CoupBjp,CoupBjm, DY11, eDY11, DY11p, DY11m);  }
  // Zinv1
  double Zinv1=0, eZinv1=0, Zinv1p=0, Zinv1m=0;
  if ( use_Zinv1 ) {    GetYield(Lumi,"Zinv1",CoupB,CoupBjp,CoupBjm, Zinv1, eZinv1, Zinv1p, Zinv1m);  }
  // Zinv2
  double Zinv2=0, eZinv2=0, Zinv2p=0, Zinv2m=0;
  if ( use_Zinv2 ) {    GetYield(Lumi,"Zinv2",CoupB,CoupBjp,CoupBjm, Zinv2, eZinv2, Zinv2p, Zinv2m);  }
  // Zinv3
  double Zinv3=0, eZinv3=0, Zinv3p=0, Zinv3m=0;
  if ( use_Zinv3 ) {    GetYield(Lumi,"Zinv3",CoupB,CoupBjp,CoupBjm, Zinv3, eZinv3, Zinv3p, Zinv3m);  }
  // Zinv4
  double Zinv4=0, eZinv4=0, Zinv4p=0, Zinv4m=0;
  if ( use_Zinv4 ) {    GetYield(Lumi,"Zinv4",CoupB,CoupBjp,CoupBjm, Zinv4, eZinv4, Zinv4p, Zinv4m);  }
  // Zinv5
  double Zinv5=0, eZinv5=0, Zinv5p=0, Zinv5m=0;
  if ( use_Zinv5 ) {    GetYield(Lumi,"Zinv5",CoupB,CoupBjp,CoupBjm, Zinv5, eZinv5, Zinv5p, Zinv5m);  }
  // Zinv6
  double Zinv6=0, eZinv6=0, Zinv6p=0, Zinv6m=0;
  if ( use_Zinv6 ) {    GetYield(Lumi,"Zinv6",CoupB,CoupBjp,CoupBjm, Zinv6, eZinv6, Zinv6p, Zinv6m);  }
  // Zinv7
  double Zinv7=0, eZinv7=0, Zinv7p=0, Zinv7m=0;
  if ( use_Zinv7 ) {    GetYield(Lumi,"Zinv7",CoupB,CoupBjp,CoupBjm, Zinv7, eZinv7, Zinv7p, Zinv7m);  }
  // ZZ
  double ZZ=0, eZZ=0, ZZp=0, ZZm=0;
  if ( use_ZZ ) {    GetYield(Lumi,"ZZ",CoupB,CoupBjp,CoupBjm, ZZ, eZZ, ZZp, ZZm);  }
  // TTX
  double TTgj=0, eTTgj=0, TTgjp=0, TTgjm=0;
  if ( use_TTgj ) { GetYield(Lumi,"TTgj",CoupB,CoupBjp,CoupBjm, TTgj, eTTgj, TTgjp, TTgjm); }   
  double TTwl=0, eTTwl=0, TTwlp=0, TTwlm=0;
  if ( use_TTwl ) { GetYield(Lumi,"TTwl",CoupB,CoupBjp,CoupBjm, TTwl, eTTwl, TTwlp, TTwlm); }   
  double TTwqq=0, eTTwqq=0, TTwqqp=0, TTwqqm=0;
  if ( use_TTwqq ) { GetYield(Lumi,"TTwqq",CoupB,CoupBjp,CoupBjm, TTwqq, eTTwqq, TTwqqp, TTwqqm); }   
  double TTzl=0, eTTzl=0, TTzlp=0, TTzlm=0;
  if ( use_TTzl ) { GetYield(Lumi,"TTzl",CoupB,CoupBjp,CoupBjm, TTzl, eTTzl, TTzlp, TTzlm); }   
  double TTzl1m10=0, eTTzl1m10=0, TTzl1m10p=0, TTzl1m10m=0;
  if ( use_TTzl1m10 ) { GetYield(Lumi,"TTzl1m10",CoupB,CoupBjp,CoupBjm, TTzl1m10, eTTzl1m10, TTzl1m10p, TTzl1m10m); }   
  double TTzqq=0, eTTzqq=0, TTzqqp=0, TTzqqm=0;
  if ( use_TTzqq ) { GetYield(Lumi,"TTzqq",CoupB,CoupBjp,CoupBjm, TTzqq, eTTzqq, TTzqqp, TTzqqm); }   

  //19- Total background
  double bckg  = W1j+W2j+W3j+W4j+W5j+W6j+W7j+ WW+WZ+ZZ +sTw+sTtC+sTsC+sTBw+sTBtC+sTBsC +DY1+DY2+DY3+DY4+DY5+DY6+DY7+DY8+DY9+DY10+DY11 
    +Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7 +TT +TTgj+TTwl+TTwqq+TTzl+TTzl1m10+TTzqq;
  double bckgp = W1jp+W2jp+W3jp+W4jp+W5jp+W6jp+W7jp+WWp+WZp+ZZp+sTwp+sTtCp+sTsCp+DY1p+DY2p+DY3p+DY4p+DY5p+DY6p+DY7p+DY8p+DY9p+DY10p+DY11p
    +Zinv1p+Zinv2p+Zinv3p+Zinv4p+Zinv5p+Zinv6p+Zinv7p+TTp+TTgjp+TTwlp+TTwqqp+TTzlp+TTzl1m10p+TTzqqp;
  double bckgm = W1jm+W2jm+W3jm+W4jm+W5jm+W6jm+W7jm+WWm+WZm+ZZm+sTwm+sTtCm+sTsCm+DY1m+DY2m+DY3m+DY4m+DY5m+DY6m+DY7m+DY8m+DY9m+DY10m+DY11m
    +Zinv1p+Zinv2p+Zinv3p+Zinv4p+Zinv5p+Zinv6p+Zinv7m+TTm+TTgjm+TTwlm+TTwqqm+TTzlm+TTzl1m10m+TTzqqm;
  double ebckg =
    sqrt(pow(eDY1,2)+pow(eW1j,2)+pow(eW2j,2)+pow(eW3j,2)+pow(eW4j,2)+pow(eW5j,2)+pow(eW6j,2)+pow(eW7j,2)+pow(eWW,2)+pow(eWZ,2)+
	 pow(esTw,2)+pow(esTtC,2)+pow(esTsC,2)+pow(esTBw,2)+pow(esTBtC,2)+pow(esTBsC,2)+pow(eDY2,2)+pow(eDY3,2)+pow(eDY4,2)+
	 pow(eZinv1,2)+pow(eZinv2,2)+pow(eZinv3,2)+pow(eZinv4,2)+pow(eZinv5,2)+pow(eZinv6,2)+pow(eZinv7,2)+pow(eTT,2)+
	 +pow(eTTgj,2)+pow(eTTwl,2)+pow(eTTwqq,2)+pow(eTTzl,2)+pow(eTTzl1m10,2)+pow(eTTzqq,2));



  tablesFile << CutName[ic] << " & "
	     << W1j+W2j+W3j+W4j+W5j+W6j+W7j          << " & "
	     << TT          << " & "
	     << Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7          << " & "
	     << bckg-TT-(W1j+W2j+W3j+W4j+W5j+W6j+W7j)-(Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7)         << " & "
	     << bckg << " $\\pm$ " << ebckg << " & "
    //             << mdata->GetEntries() << " & "
	     << SIGNAL << " $\\pm$ " << eSIGNAL << " \\\\" << endl;

  }

tablesFile << "\\hline " << endl;
tablesFile << "\\hline" << endl;
tablesFile << "\\end{tabular}" << endl;
tablesFile << "\\caption{xxx}" << endl;
//tablesFile << "\\label{tab:SingleLeptCutFlow}" << endl;
tablesFile << "\\end{center}" << endl;
tablesFile << "\\end{table}" << endl;
tablesFile << endl << endl << "\\end{document}" << endl;
tablesFile.close();


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

  TChain ProcFile("bdttree");
  // Directory of roottuples
  //  string BaseDir = "/lstore/cms/bargassa/Stop4body/SET2102_DM30/";
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt80/";
  string FFile = " ";

  // File to read
  if (pc=="SIGNAL"){     FFile = BaseDir + "T2DegStop_500_420_bdt.root";  }
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200_bdt.root";  }
  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400_bdt.root";  }
  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600_bdt.root";  }
  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800_bdt.root";  }
  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200_bdt.root";  }
  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500_bdt.root";  }
  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf_bdt.root";  }
  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200_bdt.root";  }
  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400_bdt.root";  }
  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600_bdt.root";  }
  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf_bdt.root";  }
  if (pc=="W1Jets"){     FFile = BaseDir + "Wjets_100to200_bdt.root";  }
  if (pc=="W2Jets"){     FFile = BaseDir + "Wjets_200to400_bdt.root";  }
  if (pc=="W3Jets"){     FFile = BaseDir + "Wjets_400to600_bdt.root";  }
  if (pc=="W4Jets"){     FFile = BaseDir + "Wjets_600to800_bdt.root";  }
  if (pc=="W5Jets"){     FFile = BaseDir + "Wjets_800to1200_bdt.root";  }
  if (pc=="W6Jets"){     FFile = BaseDir + "Wjets_1200to2500_bdt.root";  }
  if (pc=="W7Jets"){     FFile = BaseDir + "Wjets_2500toInf_bdt.root";  }
  if (pc=="Zinv1"){     FFile = BaseDir + "ZJetsToNuNu_HT100to200_bdt.root";  }
  if (pc=="Zinv2"){     FFile = BaseDir + "ZJetsToNuNu_HT200to400_bdt.root";  }
  if (pc=="Zinv3"){     FFile = BaseDir + "ZJetsToNuNu_HT400to600_bdt.root";  }
  if (pc=="Zinv4"){     FFile = BaseDir + "ZJetsToNuNu_HT600to800_bdt.root";  }
  if (pc=="Zinv5"){     FFile = BaseDir + "ZJetsToNuNu_HT800to1200_bdt.root";  }
  if (pc=="Zinv6"){     FFile = BaseDir + "ZJetsToNuNu_HT1200to2500_bdt.root";  }
  if (pc=="Zinv7"){     FFile = BaseDir + "ZJetsToNuNu_HT2500toInf_bdt.root";  }
  if (pc=="WW"){     FFile = BaseDir + "WW_bdt.root";  }
  if (pc=="WZ"){     FFile = BaseDir + "WZ_bdt.root";  }
  if (pc=="ZZ"){     FFile = BaseDir + "ZZ_bdt.root";  }
  if (pc=="SingleToptW"){     FFile = BaseDir + "T_tWch_ext_bdt.root";  }
  if (pc=="SingleToptChannel"){     FFile = BaseDir + "T_tch_powheg_bdt.root";  }
  if (pc=="SingleTopBtW"){     FFile = BaseDir + "TBar_tWch_ext_bdt.root";  }
  if (pc=="SingleTopBtChannel"){     FFile = BaseDir + "TBar_tch_powheg_bdt.root";  }
  ProcFile.Add((FFile).c_str());


  // Get the weighted numerator, where division by denominator has been done
  TH1D *nbcut = new TH1D("nbcut","nbcut",200,-100,100);
  nbcut->Sumw2();
  TH1D *nbcutjp = new TH1D("nbcutjp","nbcutjp",200,-100,100);
  TH1D *nbcutjm = new TH1D("nbcutjm","nbcutjm",200,-100,100);
  ProcFile.Draw("Njet>>nbcut","weight*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutjp","weight*splitFactor"*Coupurejp);
  ProcFile.Draw("Njet>>nbcutjm","weight*splitFactor"*Coupurejm);

  double NumWeight = 0.;
  double NumErr2 = 0.;
  double NumWeightP = 0.;
  double NumWeightM = 0.;
  for(int i=0; i<=(nbcut->GetNbinsX()+1); i++){ 
    NumWeight += nbcut->GetBinContent(i); 
    double NumErr = nbcut->GetBinError(i);
    NumErr2 += (NumErr*NumErr);
  }
  for(int i=0; i<=(nbcutjp->GetNbinsX()+1); i++){ NumWeightP += nbcutjp->GetBinContent(i); }
  for(int i=0; i<=(nbcutjm->GetNbinsX()+1); i++){ NumWeightM += nbcutjm->GetBinContent(i); }

  /*
  // Get the numerator
  TH1D *njcut = new TH1D("njcut","njcut",200,-100,100);
  ProcFile.Draw("Njet>>njcut",Coupure);
  double Num = njcut->GetEntries();
  // Get the denominator
  TH1D *nbt = new TH1D("nbt","nbt",10,0.,10000000.);
  ProcFile.Draw("Nevt>>nbt");
  double *Den = ProcFile.GetV1();
  */

  // Calculate yields 
  Yield  = lumi*NumWeight;
  //  eYield = lumi*uneff(Num,Den[0]);
  eYield = lumi*sqrt(NumErr2);
  YieldP = lumi*NumWeightP;
  YieldM = lumi*NumWeightM;

  // Delete histograms
  delete nbcut;
  delete nbcutjp;
  delete nbcutjm;
  //  delete njcut;
  //  delete nbt;

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
