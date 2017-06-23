// The user has to:
// * provide the luminosity Lumi
// * chose which MC contributions to use
// * modify the section "Define analysis" to modify analysis cuts

#include <iomanip>
#include <fstream>
#include "TMath.h"
#include <iostream>
#include <vector>
#include "string"
#include "TStyle.h"
#include "TCut.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TEventList.h"
#include "Tmc.h"

using namespace std;

void Tmc (){
  
  gStyle->SetOptStat(000000);

  //////////////////////////////////////////////////////////////////////////
  // Give luminosity (positive value ==> weighting whereas <=0 value ==> unweighting
  // Luminosity [pb]
  const double Lumi = 35866.0;

  //////////////////////////////////////////////////////////////////////////
  // Which MC contributions to use
  bool use_W1Jets = true;
  bool use_W2Jets = true;
  bool use_W3Jets = true;
  bool use_W4Jets = true;
  bool use_W5Jets = true;
  bool use_W6Jets = true;
  bool use_W7Jets = true;
  bool use_TT     = true;

  //////////////////////////////////////////////////////////////////////////
  cout.setf(ios::floatfield,ios::fixed);

  // Get Delta_M to cover
  string idm = " ";
  cout << "Which Delta_M to cover ?" << std::endl ;
  cin >> idm ;

  //////////////////////////////////////////////////////////////////////////
  // Define analysis

  TCut trgD   = "Event>-999.";

  TCut Promp = "isPrompt==1";
  //  TCut singlep = "(isTight==1)&&((nGoodEl + nGoodMu) <= 2)";
  TCut singlep = "(isTight==1)";
  TCut lept = "LepPt<10000.";
  TCut BDT = "BDT>-999.";
  if (idm=="10"){
    lept = "LepPt<30.";
    BDT = "BDT>0.31";
  }
  if (idm=="20"){
    lept = "LepPt<30.";
    BDT = "BDT>0.39";
  }
  if (idm=="30"){
    lept = "LepPt<30.";
    BDT = "BDT>0.47";
  }
  if (idm=="40"){
    lept = "LepPt<30.";
    BDT = "BDT>0.48";
  }
  if (idm=="50"){
    lept = "LepPt<30.";
    BDT = "BDT>0.45";
  }
  if (idm=="60"){
    lept = "LepPt<30.";
    BDT = "BDT>0.50";
  }
  if (idm=="70"){ BDT = "BDT>0.46"; }
  if (idm=="80"){ BDT = "BDT>0.44"; }

  TCut ISRjet = "Jet1Pt>110.";
  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";
  TCut met    = "Met>280.";
  TCut HT = "HT > 200.";
  TCut NoBloose = "NbLoose==0";
  TCut Bhard = "NbTight>0";
  TCut BD = "BDT < 0.2";

  TCut presel = singlep+lept+ISRjet+dphij1j2+met+HT;

  // For background: Take only the prompt component
  TCut CoupB   = Promp+presel+BDT;
  TCut CoupTT  = Promp+presel+Bhard+BD;
  TCut CoupWJ  = Promp+presel+NoBloose+BD;

  //////////////////////////////////////////////////////////////////////////

  // TTbar
  double TT=0, eTT=0, TTcr=0, eTTcr=0;
  double sISR1pTTj=0, sISR2pTTj=0, sISR3pTTj=0, sISR4pTTj=0, sISR5pTTj=0, sISR6pTTj=0, sISR7pTTj=0;
  double sISR1mTTj=0, sISR2mTTj=0, sISR3mTTj=0, sISR4mTTj=0, sISR5mTTj=0, sISR6mTTj=0, sISR7mTTj=0;
  if ( use_TT ) { 
    GetYield(Lumi,idm,"TTbar",CoupB, TT, eTT);
    GetYield(Lumi,idm,"TTbar",CoupTT, TTcr, eTTcr);
    //    GetSysISR(Lumi,idm,"TTbar",CoupB,
    //	      sISR1pTTj,sISR2pTTj,sISR3pTTj,sISR4pTTj,sISR5pTTj,sISR6pTTj,sISR7pTTj,
    //	      sISR1mTTj,sISR2mTTj,sISR3mTTj,sISR4mTTj,sISR5mTTj,sISR6mTTj,sISR7mTTj); 
  }
  // W+jets
  double W1j=0, eW1j=0, W1jcr=0, eW1jcr=0;
  double sISR1pW1j=0, sISR2pW1j=0, sISR3pW1j=0, sISR4pW1j=0, sISR5pW1j=0, sISR6pW1j=0, sISR7pW1j=0;
  double sISR1mW1j=0, sISR2mW1j=0, sISR3mW1j=0, sISR4mW1j=0, sISR5mW1j=0, sISR6mW1j=0, sISR7mW1j=0;
  if ( use_W1Jets ) { 
    GetYield(Lumi,idm,"W1Jets",CoupB, W1j, eW1j);
    GetYield(Lumi,idm,"W1Jets",CoupWJ, W1jcr, eW1jcr);
    //    GetSysISR(Lumi,idm,"W1j",CoupB,
    //	      sISR1pW1j,sISR2pW1j,sISR3pW1j,sISR4pW1j,sISR5pW1j,sISR6pW1j,sISR7pW1j,
    //	      sISR1mW1j,sISR2mW1j,sISR3mW1j,sISR4mW1j,sISR5mW1j,sISR6mW1j,sISR7mW1j); 
  }
  double W2j=0, eW2j=0, W2jcr=0, eW2jcr=0;
  double sISR1pW2j=0, sISR2pW2j=0, sISR3pW2j=0, sISR4pW2j=0, sISR5pW2j=0, sISR6pW2j=0, sISR7pW2j=0;
  double sISR1mW2j=0, sISR2mW2j=0, sISR3mW2j=0, sISR4mW2j=0, sISR5mW2j=0, sISR6mW2j=0, sISR7mW2j=0;
  if ( use_W2Jets ) { 
    GetYield(Lumi,idm,"W2Jets",CoupB, W2j, eW2j);
    GetYield(Lumi,idm,"W2Jets",CoupWJ, W2jcr, eW2jcr);
    //    GetSysISR(Lumi,idm,"W2j",CoupB,
    //	      sISR1pW2j,sISR2pW2j,sISR3pW2j,sISR4pW2j,sISR5pW2j,sISR6pW2j,sISR7pW2j,
    //	      sISR1mW2j,sISR2mW2j,sISR3mW2j,sISR4mW2j,sISR5mW2j,sISR6mW2j,sISR7mW2j); 
  }
  double W3j=0, eW3j=0, W3jcr=0, eW3jcr=0;
  double sISR1pW3j=0, sISR2pW3j=0, sISR3pW3j=0, sISR4pW3j=0, sISR5pW3j=0, sISR6pW3j=0, sISR7pW3j=0;
  double sISR1mW3j=0, sISR2mW3j=0, sISR3mW3j=0, sISR4mW3j=0, sISR5mW3j=0, sISR6mW3j=0, sISR7mW3j=0;
  if ( use_W3Jets ) { 
    GetYield(Lumi,idm,"W3Jets",CoupB, W3j, eW3j);
    GetYield(Lumi,idm,"W3Jets",CoupWJ, W3jcr, eW3jcr);
    //    GetSysISR(Lumi,idm,"W3j",CoupB,
    //	      sISR1pW3j,sISR2pW3j,sISR3pW3j,sISR4pW3j,sISR5pW3j,sISR6pW3j,sISR7pW3j,
    //	      sISR1mW3j,sISR2mW3j,sISR3mW3j,sISR4mW3j,sISR5mW3j,sISR6mW3j,sISR7mW3j); 
  }
  double W4j=0, eW4j=0, W4jcr=0, eW4jcr=0;
  double sISR1pW4j=0, sISR2pW4j=0, sISR3pW4j=0, sISR4pW4j=0, sISR5pW4j=0, sISR6pW4j=0, sISR7pW4j=0;
  double sISR1mW4j=0, sISR2mW4j=0, sISR3mW4j=0, sISR4mW4j=0, sISR5mW4j=0, sISR6mW4j=0, sISR7mW4j=0;
  if ( use_W4Jets ) { 
    GetYield(Lumi,idm,"W4Jets",CoupB, W4j, eW4j);
    GetYield(Lumi,idm,"W4Jets",CoupWJ, W4jcr, eW4jcr);
    //    GetSysISR(Lumi,idm,"W4j",CoupB,
    //	      sISR1pW4j,sISR2pW4j,sISR3pW4j,sISR4pW4j,sISR5pW4j,sISR6pW4j,sISR7pW4j,
    //	      sISR1mW4j,sISR2mW4j,sISR3mW4j,sISR4mW4j,sISR5mW4j,sISR6mW4j,sISR7mW4j); 
  }
  double W5j=0, eW5j=0, W5jcr=0, eW5jcr=0;
  double sISR1pW5j=0, sISR2pW5j=0, sISR3pW5j=0, sISR4pW5j=0, sISR5pW5j=0, sISR6pW5j=0, sISR7pW5j=0;
  double sISR1mW5j=0, sISR2mW5j=0, sISR3mW5j=0, sISR4mW5j=0, sISR5mW5j=0, sISR6mW5j=0, sISR7mW5j=0;
  if ( use_W5Jets ) { 
    GetYield(Lumi,idm,"W5Jets",CoupB, W5j, eW5j);
    GetYield(Lumi,idm,"W5Jets",CoupWJ, W5jcr, eW5jcr);
    //    GetSysISR(Lumi,idm,"W5j",CoupB,
    //	      sISR1pW5j,sISR2pW5j,sISR3pW5j,sISR4pW5j,sISR5pW5j,sISR6pW5j,sISR7pW5j,
    //	      sISR1mW5j,sISR2mW5j,sISR3mW5j,sISR4mW5j,sISR5mW5j,sISR6mW5j,sISR7mW5j); 
  }
  double W6j=0, eW6j=0, W6jcr=0, eW6jcr=0;
  double sISR1pW6j=0, sISR2pW6j=0, sISR3pW6j=0, sISR4pW6j=0, sISR5pW6j=0, sISR6pW6j=0, sISR7pW6j=0;
  double sISR1mW6j=0, sISR2mW6j=0, sISR3mW6j=0, sISR4mW6j=0, sISR5mW6j=0, sISR6mW6j=0, sISR7mW6j=0;
  if ( use_W6Jets ) { 
    GetYield(Lumi,idm,"W6Jets",CoupB, W6j, eW6j);
    GetYield(Lumi,idm,"W6Jets",CoupWJ, W6jcr, eW6jcr);
    //    GetSysISR(Lumi,idm,"W6j",CoupB,
    //	      sISR1pW6j,sISR2pW6j,sISR3pW6j,sISR4pW6j,sISR5pW6j,sISR6pW6j,sISR7pW6j,
    //	      sISR1mW6j,sISR2mW6j,sISR3mW6j,sISR4mW6j,sISR5mW6j,sISR6mW6j,sISR7mW6j); 
  }
  double W7j=0, eW7j=0, W7jcr=0, eW7jcr=0;
  double sISR1pW7j=0, sISR2pW7j=0, sISR3pW7j=0, sISR4pW7j=0, sISR5pW7j=0, sISR6pW7j=0, sISR7pW7j=0;
  double sISR1mW7j=0, sISR2mW7j=0, sISR3mW7j=0, sISR4mW7j=0, sISR5mW7j=0, sISR6mW7j=0, sISR7mW7j=0;
  if ( use_W7Jets ) { 
    GetYield(Lumi,idm,"W7Jets",CoupB, W7j, eW7j);
    GetYield(Lumi,idm,"W7Jets",CoupWJ, W7jcr, eW7jcr);
    //    GetSysISR(Lumi,idm,"W7j",CoupB,
    //	      sISR1pW7j,sISR2pW7j,sISR3pW7j,sISR4pW7j,sISR5pW7j,sISR6pW7j,sISR7pW7j,
    //	      sISR1mW7j,sISR2mW7j,sISR3mW7j,sISR4mW7j,sISR5mW7j,sISR6mW7j,sISR7mW7j); 
  }

  // Backgrounds
  double WjT = W1j+W2j+W3j+W4j+W5j+W6j+W7j;
  double WjTcr = W1jcr+W2jcr+W3jcr+W4jcr+W5jcr+W6jcr+W7jcr;

  cout << "T_MC(Wjets) = " << eff(WjT,WjTcr) << std::endl;
  cout << "T_MC(ttbar) = " << eff(TT,TTcr)   << std::endl;
  //  cout << "T_MC_ISR1(Wjets) = " << eff(,) << std::endl;


}


void GetYield(
	      double lumi,     // Input: Luminosity
	      string dm,
	      string pc,       // Input: Name of the channel
	      TCut Coupure,    // Input: Analysis cut
	      double& Yield,   // Output: Weighted yield of the channel
	      double& eYield   // Output: Error on the weighted yield
	      ){ 

  TChain ProcFile("bdttree");
  // Directory of roottuples
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt" + dm + "/";
  //  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_SysVar/";
  string FFile = " ";

  // File to read
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  if (pc=="W1Jets"){     FFile = BaseDir + "Wjets_100to200_bdt.root";  }
  if (pc=="W2Jets"){     FFile = BaseDir + "Wjets_200to400_bdt.root";  }
  if (pc=="W3Jets"){     FFile = BaseDir + "Wjets_400to600_bdt.root";  }
  if (pc=="W4Jets"){     FFile = BaseDir + "Wjets_600to800_bdt.root";  }
  if (pc=="W5Jets"){     FFile = BaseDir + "Wjets_800to1200_bdt.root";  }
  if (pc=="W6Jets"){     FFile = BaseDir + "Wjets_1200to2500_bdt.root";  }
  if (pc=="W7Jets"){     FFile = BaseDir + "Wjets_2500toInf_bdt.root";  }
  ProcFile.Add((FFile).c_str());
  
  // Get the weighted numerator, where division by denominator has been done
  TH1D *nbcut = new TH1D("nbcut","nbcut",200,-100,100);
  nbcut->Sumw2();
  ProcFile.Draw("Njet>>nbcut","weight*splitFactor"*Coupure);
  //
  double NumWeight = 0., NumErr2 = 0.;
  for(int i=0; i<=(nbcut->GetNbinsX()+1); i++){ 
    NumWeight += nbcut->GetBinContent(i); 
    double NumErr = nbcut->GetBinError(i);
    NumErr2 += (NumErr*NumErr);
  }
  // Calculate yields 
  Yield  = lumi*NumWeight;
  eYield = lumi*sqrt(NumErr2);

  // Delete histograms
  delete nbcut;

}


void GetSysISR(
		  double lumi,     // Input: Luminosity
		  string dm,
		  string pc,       // Input: Name of the channel
		  TCut Coupure,    // Input: Analysis cut
		  double& YISR1p,  // 
		  double& YISR2p,  // 
		  double& YISR3p,  // 
		  double& YISR4p,  // 
		  double& YISR5p,  // 
		  double& YISR6p,  // 
		  double& YISR7p,  // 
		  double& YISR1m,  // 
		  double& YISR2m,  // 
		  double& YISR3m,  // 
		  double& YISR4m,  // 
		  double& YISR5m,  // 
		  double& YISR6m,  // 
		  double& YISR7m   // 
		  ){ 

  TChain ProcFile("bdttree");
  // Directory of roottuples
  //  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt" + dm + "/";
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_SysVar/";
  string FFile = " ";

  // File to read
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  if (pc=="W1Jets"){     FFile = BaseDir + "Wjets_100to200_bdt.root";  }
  if (pc=="W2Jets"){     FFile = BaseDir + "Wjets_200to400_bdt.root";  }
  if (pc=="W3Jets"){     FFile = BaseDir + "Wjets_400to600_bdt.root";  }
  if (pc=="W4Jets"){     FFile = BaseDir + "Wjets_600to800_bdt.root";  }
  if (pc=="W5Jets"){     FFile = BaseDir + "Wjets_800to1200_bdt.root";  }
  if (pc=="W6Jets"){     FFile = BaseDir + "Wjets_1200to2500_bdt.root";  }
  if (pc=="W7Jets"){     FFile = BaseDir + "Wjets_2500toInf_bdt.root";  }
  ProcFile.Add((FFile).c_str());
  
  // sysISR
  TH1D *sISR1p = new TH1D("sISR1p","sISR1p",200,-100,100);
  TH1D *sISR2p = new TH1D("sISR2p","sISR1p",200,-100,100);
  TH1D *sISR3p = new TH1D("sISR3p","sISR1p",200,-100,100);
  TH1D *sISR4p = new TH1D("sISR4p","sISR1p",200,-100,100);
  TH1D *sISR5p = new TH1D("sISR5p","sISR1p",200,-100,100);
  TH1D *sISR6p = new TH1D("sISR6p","sISR1p",200,-100,100);
  TH1D *sISR7p = new TH1D("sISR7p","sISR1p",200,-100,100);
  TH1D *sISR1m = new TH1D("sISR1m","sISR1m",200,-100,100);
  TH1D *sISR2m = new TH1D("sISR2m","sISR1m",200,-100,100);
  TH1D *sISR3m = new TH1D("sISR3m","sISR1m",200,-100,100);
  TH1D *sISR4m = new TH1D("sISR4m","sISR1m",200,-100,100);
  TH1D *sISR5m = new TH1D("sISR5m","sISR1m",200,-100,100);
  TH1D *sISR6m = new TH1D("sISR6m","sISR1m",200,-100,100);
  TH1D *sISR7m = new TH1D("sISR7m","sISR1m",200,-100,100);
  if ((pc=="W1Jets")||(pc=="W2Jets")||(pc=="W3Jets")||(pc=="W4Jets")
      ||(pc=="W5Jets")||(pc=="W6Jets")||(pc=="W7Jets")){
    ProcFile.Draw("Njet>>sISR1p","weight_EWKISRweight_Bin1_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR2p","weight_EWKISRweight_Bin2_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR3p","weight_EWKISRweight_Bin3_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR4p","weight_EWKISRweight_Bin4_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR5p","weight_EWKISRweight_Bin5_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR6p","weight_EWKISRweight_Bin6_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR7p","weight_EWKISRweight_Bin7_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR1m","weight_EWKISRweight_Bin1_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR2m","weight_EWKISRweight_Bin2_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR3m","weight_EWKISRweight_Bin3_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR4m","weight_EWKISRweight_Bin4_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR5m","weight_EWKISRweight_Bin5_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR6m","weight_EWKISRweight_Bin6_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR7m","weight_EWKISRweight_Bin7_Down*splitFactor"*Coupure);
  }
  else if (pc=="TTbar"){
    ProcFile.Draw("Njet>>sISR1p","weight_ISRweight_Bin1_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR2p","weight_ISRweight_Bin2_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR3p","weight_ISRweight_Bin3_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR4p","weight_ISRweight_Bin4_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR5p","weight_ISRweight_Bin5_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR6p","weight_ISRweight_Bin6_Up*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR1m","weight_ISRweight_Bin1_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR2m","weight_ISRweight_Bin2_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR3m","weight_ISRweight_Bin3_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR4m","weight_ISRweight_Bin4_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR5m","weight_ISRweight_Bin5_Down*splitFactor"*Coupure);
    ProcFile.Draw("Njet>>sISR6m","weight_ISRweight_Bin6_Down*splitFactor"*Coupure);
  }
  double ISR1p=0,ISR2p=0,ISR3p=0,ISR4p=0,ISR5p=0,ISR6p=0,ISR7p=0;
  double ISR1m=0,ISR2m=0,ISR3m=0,ISR4m=0,ISR5m=0,ISR6m=0,ISR7m=0;
  for(int i=0; i<=(sISR1p->GetNbinsX()+1); i++){ 
    ISR1p += sISR1p->GetBinContent(i);
    ISR2p += sISR2p->GetBinContent(i);
    ISR3p += sISR3p->GetBinContent(i);
    ISR4p += sISR4p->GetBinContent(i);
    ISR5p += sISR5p->GetBinContent(i);
    ISR6p += sISR6p->GetBinContent(i);
    ISR7p += sISR7p->GetBinContent(i);
    ISR1m += sISR1m->GetBinContent(i);
    ISR2m += sISR2m->GetBinContent(i);
    ISR3m += sISR3m->GetBinContent(i);
    ISR4m += sISR4m->GetBinContent(i);
    ISR5m += sISR5m->GetBinContent(i);
    ISR6m += sISR6m->GetBinContent(i);
    ISR7m += sISR7m->GetBinContent(i);
  }
  YISR1p = lumi*ISR1p;
  YISR2p = lumi*ISR2p;
  YISR3p = lumi*ISR3p;
  YISR4p = lumi*ISR4p;
  YISR5p = lumi*ISR5p;
  YISR6p = lumi*ISR6p;
  YISR7p = lumi*ISR7p;
  YISR1m = lumi*ISR1m;
  YISR2m = lumi*ISR2m;
  YISR3m = lumi*ISR3m;
  YISR4m = lumi*ISR4m;
  YISR5m = lumi*ISR5m;
  YISR6m = lumi*ISR6m;
  YISR7m = lumi*ISR7m;

  // Delete histograms
  delete sISR1p;
  delete sISR2p;
  delete sISR3p;
  delete sISR4p;
  delete sISR5p;
  delete sISR6p;
  delete sISR7p;
  delete sISR1m;
  delete sISR2m;
  delete sISR3m;
  delete sISR4m;
  delete sISR5m;
  delete sISR6m;
  delete sISR7m;

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
