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
#include "DataCard_Factory.h"

using namespace std;

void DataCard_Factory (){
  
  gStyle->SetOptStat(000000);

  //////////////////////////////////////////////////////////////////////////
  // Give luminosity (positive value ==> weighting whereas <=0 value ==> unweighting
  // Luminosity [pb]
  const double Lumi = 35866.0;

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
  bool use_W1Jets       = false;
  bool use_W2Jets       = false;
  bool use_W3Jets       = false;
  bool use_W4Jets       = false;
  bool use_W5Jets       = false;
  bool use_W6Jets       = false;
  bool use_W7Jets       = false;
  bool use_Zinv1        = false;
  bool use_Zinv2        = false;
  bool use_Zinv3        = false;
  bool use_Zinv4        = false;
  bool use_Zinv5        = false;
  bool use_Zinv6        = false;
  bool use_Zinv7        = false;
  bool use_WW           = true;
  bool use_WZ           = true;
  bool use_ZZ           = true;
  bool use_SingleToptW  = false;
  bool use_SingleToptC  = false;
  bool use_SingleTopBtW = false;
  bool use_SingleTopBtC = false;
  bool use_TT           = false;
  bool use_TTgj         = false;
  bool use_TTwl         = false;
  bool use_TTwqq        = false;
  bool use_TTzl         = false;
  bool use_TTzl1m10     = false;
  bool use_TTzqq        = false;
  bool use_Data         = false;

  //////////////////////////////////////////////////////////////////////////
  cout.setf(ios::floatfield,ios::fixed);
  //  cout<<setprecision(3);

  // Get Delta_M to cover and corresponding signal points
  string idm = " ";
  cout << "Which Delta_M to cover ?" << std::endl ;
  cin >> idm ;
  string dm = "DM" + idm;
  ifstream dmFile;
  dmFile.open(dm + ".txt",ifstream::in);
  string dmline;
  std::vector<std::string> dmlines;
  while (dmFile.good()) {
    dmFile >> dmline;
    dmlines.push_back(dmline);
  }
  dmFile.close();

  // Read the DD backgrounds
  string dd = "DD" + idm;
  ifstream ddFile;
  ddFile.open(dd + ".txt",ifstream::in);
  string ddline;
  std::vector<std::string> ddlines;
  while (ddFile.good()) {
    ddFile >> ddline;
    ddlines.push_back(ddline);
  }
  ddFile.close();
  double DDwj = std::stod(ddlines[0]);
  double eDDwj = std::stod(ddlines[1]);
  double DDtt = std::stod(ddlines[2]);
  double eDDtt = std::stod(ddlines[3]);
  double DDfk = std::stod(ddlines[4]);
  double eDDfk = std::stod(ddlines[5]);

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

  // 3 JES-varied serie of cuts on jet-related variables
  // Central JES
  TCut ISRjet = "Jet1Pt>110.";
  TCut dphij1j2 = "(DPhiJet1Jet2 < 2.5)||(Jet2Pt<60.)";
  TCut met    = "Met>280.";
  TCut HT = "HT > 200.";
  // +1 JES
  TCut ISRjetjp = "Jet1Pt>110.";
  TCut metjp    = "Met>280.";
  TCut HTjp = "HT > 200.";
  // -1 JES
  TCut ISRjetjm = "Jet1Pt>110.";
  TCut metjm    = "Met>280.";
  TCut HTjm = "HT > 200.";

  // 3 JES-varied preselection
  TCut presel = singlep+lept+ISRjet+dphij1j2+met+HT;
  TCut preseljp = singlep+lept+ISRjetjp+dphij1j2+metjp+HTjp;
  TCut preseljm = singlep+lept+ISRjetjm+dphij1j2+metjm+HTjm;

  // 3 JES-varied final selection
  TCut Coup   = presel;
  TCut Coupjp = preseljp;
  TCut Coupjm = preseljm;
  // For background: Take only the prompt component
  TCut CoupB   = Promp+presel;
  TCut CoupBjp = Promp+preseljp;
  TCut CoupBjm = Promp+preseljm;

  //////////////////////////////////////////////////////////////////////////

  // TTbar
  double TT=0, eTT=0, TTp=0, TTm=0, TTpup=0, TTpum=0;
  if ( use_TT ) { GetYield(Lumi,idm,"TTbar",CoupB,CoupBjp,CoupBjm, TT, eTT, TTp, TTm, TTpup, TTpum); }   
  // W+jets
  double W1j=0, eW1j=0, W1jp=0, W1jm=0, W1jpup=0, W1jpum=0;
  if ( use_W1Jets ) { GetYield(Lumi,idm,"W1Jets",CoupB,CoupBjp,CoupBjm, W1j, eW1j, W1jp, W1jm, W1jpup, W1jpum); }
  double W2j=0, eW2j=0, W2jp=0, W2jm=0, W2jpup=0, W2jpum=0;
  if ( use_W2Jets ) { GetYield(Lumi,idm,"W2Jets",CoupB,CoupBjp,CoupBjm, W2j, eW2j, W2jp, W2jm, W2jpup, W2jpum); }
  double W3j=0, eW3j=0, W3jp=0, W3jm=0, W3jpup=0, W3jpum=0;
  if ( use_W3Jets ) { GetYield(Lumi,idm,"W3Jets",CoupB,CoupBjp,CoupBjm, W3j, eW3j, W3jp, W3jm, W3jpup, W3jpum); }
  double W4j=0, eW4j=0, W4jp=0, W4jm=0, W4jpup=0, W4jpum=0;
  if ( use_W4Jets ) { GetYield(Lumi,idm,"W4Jets",CoupB,CoupBjp,CoupBjm, W4j, eW4j, W4jp, W4jm, W4jpup, W4jpum); }
  double W5j=0, eW5j=0, W5jp=0, W5jm=0, W5jpup=0, W5jpum=0;
  if ( use_W5Jets ) { GetYield(Lumi,idm,"W5Jets",CoupB,CoupBjp,CoupBjm, W5j, eW5j, W5jp, W5jm, W5jpup, W5jpum); }
  double W6j=0, eW6j=0, W6jp=0, W6jm=0, W6jpup=0, W6jpum=0;
  if ( use_W6Jets ) { GetYield(Lumi,idm,"W6Jets",CoupB,CoupBjp,CoupBjm, W6j, eW6j, W6jp, W6jm, W6jpup, W6jpum); }
  double W7j=0, eW7j=0, W7jp=0, W7jm=0, W7jpup=0, W7jpum=0;
  if ( use_W7Jets ) { GetYield(Lumi,idm,"W7Jets",CoupB,CoupBjp,CoupBjm, W7j, eW7j, W7jp, W7jm, W7jpup, W7jpum); }
  // Zinv1
  double Zinv1=0, eZinv1=0, Zinv1p=0, Zinv1m=0, Zinv1pup=0, Zinv1pum=0;
  if ( use_Zinv1 ) { GetYield(Lumi,idm,"Zinv1",CoupB,CoupBjp,CoupBjm, Zinv1, eZinv1, Zinv1p, Zinv1m, Zinv1pup, Zinv1pum); }  
  // Zinv2
  double Zinv2=0, eZinv2=0, Zinv2p=0, Zinv2m=0, Zinv2pup=0, Zinv2pum=0;
  if ( use_Zinv2 ) { GetYield(Lumi,idm,"Zinv2",CoupB,CoupBjp,CoupBjm, Zinv2, eZinv2, Zinv2p, Zinv2m, Zinv2pup, Zinv2pum); }  
  // Zinv3
  double Zinv3=0, eZinv3=0, Zinv3p=0, Zinv3m=0, Zinv3pup=0, Zinv3pum=0;
  if ( use_Zinv3 ) { GetYield(Lumi,idm,"Zinv3",CoupB,CoupBjp,CoupBjm, Zinv3, eZinv3, Zinv3p, Zinv3m, Zinv3pup, Zinv3pum); }
  // Zinv4
  double Zinv4=0, eZinv4=0, Zinv4p=0, Zinv4m=0, Zinv4pup=0, Zinv4pum=0;
  if ( use_Zinv4 ) { GetYield(Lumi,idm,"Zinv4",CoupB,CoupBjp,CoupBjm, Zinv4, eZinv4, Zinv4p, Zinv4m, Zinv4pup, Zinv4pum); }
  // Zinv5
  double Zinv5=0, eZinv5=0, Zinv5p=0, Zinv5m=0, Zinv5pup=0, Zinv5pum=0;
  if ( use_Zinv5 ) { GetYield(Lumi,idm,"Zinv5",CoupB,CoupBjp,CoupBjm, Zinv5, eZinv5, Zinv5p, Zinv5m, Zinv5pup, Zinv5pum); }
  // Zinv6
  double Zinv6=0, eZinv6=0, Zinv6p=0, Zinv6m=0, Zinv6pup=0, Zinv6pum=0;
  if ( use_Zinv6 ) { GetYield(Lumi,idm,"Zinv6",CoupB,CoupBjp,CoupBjm, Zinv6, eZinv6, Zinv6p, Zinv6m, Zinv6pup, Zinv6pum); }
  // Zinv7
  double Zinv7=0, eZinv7=0, Zinv7p=0, Zinv7m=0, Zinv7pup=0, Zinv7pum=0;
  if ( use_Zinv7 ) { GetYield(Lumi,idm,"Zinv7",CoupB,CoupBjp,CoupBjm, Zinv7, eZinv7, Zinv7p, Zinv7m, Zinv7pup, Zinv7pum); }
  // WW
  double WW=0, eWW=0, WWp=0, WWm=0, WWpup=0, WWpum=0;
  if ( use_WW ) {     GetYield(Lumi,idm,"WW",CoupB,CoupBjp,CoupBjm, WW, eWW, WWp, WWm, WWpup, WWpum);   }
  // WZ
  double WZ=0, eWZ=0, WZp=0, WZm=0, WZpup=0, WZpum=0;
  if ( use_WZ ) {     GetYield(Lumi,idm,"WZ",CoupB,CoupBjp,CoupBjm, WZ, eWZ, WZp, WZm, WZpup, WZpum);   }
  // ZZ
  double ZZ=0, eZZ=0, ZZp=0, ZZm=0, ZZpup=0, ZZpum=0;
  if ( use_ZZ ) {     GetYield(Lumi,idm,"ZZ",CoupB,CoupBjp,CoupBjm, ZZ, eZZ, ZZp, ZZm, ZZpup, ZZpum);   }
  // SingleToptW
  double sTw=0, esTw=0, sTwp=0, sTwm=0, sTwpup=0, sTwpum=0;
  if ( use_SingleToptW ) { GetYield(Lumi,idm,"SingleToptW",CoupB,CoupBjp,CoupBjm, sTw, esTw, sTwp, sTwm, sTwpup, sTwpum); }
  // SingleToptChannel
  double sTtC=0, esTtC=0, sTtCp=0, sTtCm=0, sTtCpup=0, sTtCpum=0;
  if ( use_SingleToptC) { GetYield(Lumi,idm,"SingleToptChannel",CoupB,CoupBjp,CoupBjm, sTtC, esTtC, sTtCp, sTtCm, sTtCpup, sTtCpum); }
  // SingleToptW
  double sTBw=0, esTBw=0, sTBwp=0, sTBwm=0, sTBwpup=0, sTBwpum=0;
  if ( use_SingleTopBtW ) { GetYield(Lumi,idm,"SingleTopBtW",CoupB,CoupBjp,CoupBjm, sTBw, esTBw, sTBwp, sTBwm, sTBwpup, sTBwpum); }
  // SingleToptChannel
  double sTBtC=0, esTBtC=0, sTBtCp=0, sTBtCm=0, sTBtCpup=0, sTBtCpum=0;
  if ( use_SingleTopBtC) { GetYield(Lumi,idm,"SingleTopBtChannel",CoupB,CoupBjp,CoupBjm, sTBtC, esTBtC, sTBtCp, sTBtCm, sTBtCpup, sTBtCpum); }
  // DY1
  double DY1=0, eDY1=0, DY1p=0, DY1m=0, DY1pup=0, DY1pum=0;
  double sISR1pDY1=0, sISR2pDY1=0, sISR3pDY1=0, sISR4pDY1=0, sISR5pDY1=0, sISR6pDY1=0, sISR7pDY1=0;
  double sISR1mDY1=0, sISR2mDY1=0, sISR3mDY1=0, sISR4mDY1=0, sISR5mDY1=0, sISR6mDY1=0, sISR7mDY1=0;
  if ( use_DY1 ) { 
    GetYield(Lumi,idm,"DY1",CoupB,CoupBjp,CoupBjm, DY1, eDY1, DY1p, DY1m, DY1pup, DY1pum); 
    GetSysISR(Lumi,idm,"DY1",CoupB,
	      sISR1pDY1,sISR2pDY1,sISR3pDY1,sISR4pDY1,sISR5pDY1,sISR6pDY1,sISR7pDY1,
	      sISR1mDY1,sISR2mDY1,sISR3mDY1,sISR4mDY1,sISR5mDY1,sISR6mDY1,sISR7mDY1); 
  }
  // DY2
  double DY2=0, eDY2=0, DY2p=0, DY2m=0, DY2pup=0, DY2pum=0;
  double sISR1pDY2=0, sISR2pDY2=0, sISR3pDY2=0, sISR4pDY2=0, sISR5pDY2=0, sISR6pDY2=0, sISR7pDY2=0;
  double sISR1mDY2=0, sISR2mDY2=0, sISR3mDY2=0, sISR4mDY2=0, sISR5mDY2=0, sISR6mDY2=0, sISR7mDY2=0;
  if ( use_DY2 ) { 
    GetYield(Lumi,idm,"DY2",CoupB,CoupBjp,CoupBjm, DY2, eDY2, DY2p, DY2m, DY2pup, DY2pum); 
    GetSysISR(Lumi,idm,"DY2",CoupB,
	      sISR1pDY2,sISR2pDY2,sISR3pDY2,sISR4pDY2,sISR5pDY2,sISR6pDY2,sISR7pDY2,
	      sISR1mDY2,sISR2mDY2,sISR3mDY2,sISR4mDY2,sISR5mDY2,sISR6mDY2,sISR7mDY2); 
  }
  // DY3
  double DY3=0, eDY3=0, DY3p=0, DY3m=0, DY3pup=0, DY3pum=0;
  double sISR1pDY3=0, sISR2pDY3=0, sISR3pDY3=0, sISR4pDY3=0, sISR5pDY3=0, sISR6pDY3=0, sISR7pDY3=0;
  double sISR1mDY3=0, sISR2mDY3=0, sISR3mDY3=0, sISR4mDY3=0, sISR5mDY3=0, sISR6mDY3=0, sISR7mDY3=0;
  if ( use_DY3 ) { 
    GetYield(Lumi,idm,"DY3",CoupB,CoupBjp,CoupBjm, DY3, eDY3, DY3p, DY3m, DY3pup, DY3pum); 
    GetSysISR(Lumi,idm,"DY3",CoupB,
	      sISR1pDY3,sISR2pDY3,sISR3pDY3,sISR4pDY3,sISR5pDY3,sISR6pDY3,sISR7pDY3,
	      sISR1mDY3,sISR2mDY3,sISR3mDY3,sISR4mDY3,sISR5mDY3,sISR6mDY3,sISR7mDY3); 
  }
  // DY4
  double DY4=0, eDY4=0, DY4p=0, DY4m=0, DY4pup=0, DY4pum=0;
  double sISR1pDY4=0, sISR2pDY4=0, sISR3pDY4=0, sISR4pDY4=0, sISR5pDY4=0, sISR6pDY4=0, sISR7pDY4=0;
  double sISR1mDY4=0, sISR2mDY4=0, sISR3mDY4=0, sISR4mDY4=0, sISR5mDY4=0, sISR6mDY4=0, sISR7mDY4=0;
  if ( use_DY4 ) { 
    GetYield(Lumi,idm,"DY4",CoupB,CoupBjp,CoupBjm, DY4, eDY4, DY4p, DY4m, DY4pup, DY4pum); 
    GetSysISR(Lumi,idm,"DY4",CoupB,
	      sISR1pDY4,sISR2pDY4,sISR3pDY4,sISR4pDY4,sISR5pDY4,sISR6pDY4,sISR7pDY4,
	      sISR1mDY4,sISR2mDY4,sISR3mDY4,sISR4mDY4,sISR5mDY4,sISR6mDY4,sISR7mDY4); 
  }
  // DY5
  double DY5=0, eDY5=0, DY5p=0, DY5m=0, DY5pup=0, DY5pum=0;
  double sISR1pDY5=0, sISR2pDY5=0, sISR3pDY5=0, sISR4pDY5=0, sISR5pDY5=0, sISR6pDY5=0, sISR7pDY5=0;
  double sISR1mDY5=0, sISR2mDY5=0, sISR3mDY5=0, sISR4mDY5=0, sISR5mDY5=0, sISR6mDY5=0, sISR7mDY5=0;
  if ( use_DY5 ) { 
    GetYield(Lumi,idm,"DY5",CoupB,CoupBjp,CoupBjm, DY5, eDY5, DY5p, DY5m, DY5pup, DY5pum); 
    GetSysISR(Lumi,idm,"DY5",CoupB,
	      sISR1pDY5,sISR2pDY5,sISR3pDY5,sISR4pDY5,sISR5pDY5,sISR6pDY5,sISR7pDY5,
	      sISR1mDY5,sISR2mDY5,sISR3mDY5,sISR4mDY5,sISR5mDY5,sISR6mDY5,sISR7mDY5); 
  }
  // DY6
  double DY6=0, eDY6=0, DY6p=0, DY6m=0, DY6pup=0, DY6pum=0;
  double sISR1pDY6=0, sISR2pDY6=0, sISR3pDY6=0, sISR4pDY6=0, sISR5pDY6=0, sISR6pDY6=0, sISR7pDY6=0;
  double sISR1mDY6=0, sISR2mDY6=0, sISR3mDY6=0, sISR4mDY6=0, sISR5mDY6=0, sISR6mDY6=0, sISR7mDY6=0;
  if ( use_DY6 ) { 
    GetYield(Lumi,idm,"DY6",CoupB,CoupBjp,CoupBjm, DY6, eDY6, DY6p, DY6m, DY6pup, DY6pum); 
    GetSysISR(Lumi,idm,"DY6",CoupB,
	      sISR1pDY6,sISR2pDY6,sISR3pDY6,sISR4pDY6,sISR5pDY6,sISR6pDY6,sISR7pDY6,
	      sISR1mDY6,sISR2mDY6,sISR3mDY6,sISR4mDY6,sISR5mDY6,sISR6mDY6,sISR7mDY6); 
  }
  // DY7
  double DY7=0, eDY7=0, DY7p=0, DY7m=0, DY7pup=0, DY7pum=0;
  double sISR1pDY7=0, sISR2pDY7=0, sISR3pDY7=0, sISR4pDY7=0, sISR5pDY7=0, sISR6pDY7=0, sISR7pDY7=0;
  double sISR1mDY7=0, sISR2mDY7=0, sISR3mDY7=0, sISR4mDY7=0, sISR5mDY7=0, sISR6mDY7=0, sISR7mDY7=0;
  if ( use_DY7 ) { 
    GetYield(Lumi,idm,"DY7",CoupB,CoupBjp,CoupBjm, DY7, eDY7, DY7p, DY7m, DY7pup, DY7pum); 
    GetSysISR(Lumi,idm,"DY7",CoupB,
	      sISR1pDY7,sISR2pDY7,sISR3pDY7,sISR4pDY7,sISR5pDY7,sISR6pDY7,sISR7pDY7,
	      sISR1mDY7,sISR2mDY7,sISR3mDY7,sISR4mDY7,sISR5mDY7,sISR6mDY7,sISR7mDY7); 
  }
  // DY8
  double DY8=0, eDY8=0, DY8p=0, DY8m=0, DY8pup=0, DY8pum=0;
  double sISR1pDY8=0, sISR2pDY8=0, sISR3pDY8=0, sISR4pDY8=0, sISR5pDY8=0, sISR6pDY8=0, sISR7pDY8=0;
  double sISR1mDY8=0, sISR2mDY8=0, sISR3mDY8=0, sISR4mDY8=0, sISR5mDY8=0, sISR6mDY8=0, sISR7mDY8=0;
  if ( use_DY8 ) { 
    GetYield(Lumi,idm,"DY8",CoupB,CoupBjp,CoupBjm, DY8, eDY8, DY8p, DY8m, DY8pup, DY8pum); 
    GetSysISR(Lumi,idm,"DY8",CoupB,
	      sISR1pDY8,sISR2pDY8,sISR3pDY8,sISR4pDY8,sISR5pDY8,sISR6pDY8,sISR7pDY8,
	      sISR1mDY8,sISR2mDY8,sISR3mDY8,sISR4mDY8,sISR5mDY8,sISR6mDY8,sISR7mDY8); 
  }
  // DY9
  double DY9=0, eDY9=0, DY9p=0, DY9m=0, DY9pup=0, DY9pum=0;
  double sISR1pDY9=0, sISR2pDY9=0, sISR3pDY9=0, sISR4pDY9=0, sISR5pDY9=0, sISR6pDY9=0, sISR7pDY9=0;
  double sISR1mDY9=0, sISR2mDY9=0, sISR3mDY9=0, sISR4mDY9=0, sISR5mDY9=0, sISR6mDY9=0, sISR7mDY9=0;
  if ( use_DY9 ) { 
    GetYield(Lumi,idm,"DY9",CoupB,CoupBjp,CoupBjm, DY9, eDY9, DY9p, DY9m, DY9pup, DY9pum); 
    GetSysISR(Lumi,idm,"DY9",CoupB,
	      sISR1pDY9,sISR2pDY9,sISR3pDY9,sISR4pDY9,sISR5pDY9,sISR6pDY9,sISR7pDY9,
	      sISR1mDY9,sISR2mDY9,sISR3mDY9,sISR4mDY9,sISR5mDY9,sISR6mDY9,sISR7mDY9); 
  }
  // DY10
  double DY10=0, eDY10=0, DY10p=0, DY10m=0, DY10pup=0, DY10pum=0;
  double sISR1pDY10=0, sISR2pDY10=0, sISR3pDY10=0, sISR4pDY10=0, sISR5pDY10=0, sISR6pDY10=0, sISR7pDY10=0;
  double sISR1mDY10=0, sISR2mDY10=0, sISR3mDY10=0, sISR4mDY10=0, sISR5mDY10=0, sISR6mDY10=0, sISR7mDY10=0;
  if ( use_DY10 ) { 
    GetYield(Lumi,idm,"DY10",CoupB,CoupBjp,CoupBjm, DY10, eDY10, DY10p, DY10m, DY10pup, DY10pum); 
    GetSysISR(Lumi,idm,"DY10",CoupB,
	      sISR1pDY10,sISR2pDY10,sISR3pDY10,sISR4pDY10,sISR5pDY10,sISR6pDY10,sISR7pDY10,
	      sISR1mDY10,sISR2mDY10,sISR3mDY10,sISR4mDY10,sISR5mDY10,sISR6mDY10,sISR7mDY10); 
  }
  // DY11
  double DY11=0, eDY11=0, DY11p=0, DY11m=0, DY11pup=0, DY11pum=0;
  double sISR1pDY11=0, sISR2pDY11=0, sISR3pDY11=0, sISR4pDY11=0, sISR5pDY11=0, sISR6pDY11=0, sISR7pDY11=0;
  double sISR1mDY11=0, sISR2mDY11=0, sISR3mDY11=0, sISR4mDY11=0, sISR5mDY11=0, sISR6mDY11=0, sISR7mDY11=0;
  if ( use_DY11 ) { 
    GetYield(Lumi,idm,"DY11",CoupB,CoupBjp,CoupBjm, DY11, eDY11, DY11p, DY11m, DY11pup, DY11pum); 
    GetSysISR(Lumi,idm,"DY11",CoupB,
	      sISR1pDY11,sISR2pDY11,sISR3pDY11,sISR4pDY11,sISR5pDY11,sISR6pDY11,sISR7pDY11,
	      sISR1mDY11,sISR2mDY11,sISR3mDY11,sISR4mDY11,sISR5mDY11,sISR6mDY11,sISR7mDY11); 
  }
  // TTX
  double TTgj=0, eTTgj=0, TTgjp=0, TTgjm=0, TTgjpup=0, TTgjpum=0;
  if ( use_TTgj ) { GetYield(Lumi,idm,"TTgj",CoupB,CoupBjp,CoupBjm, TTgj, eTTgj, TTgjp, TTgjm, TTgjpup, TTgjpum); }   
  double TTwl=0, eTTwl=0, TTwlp=0, TTwlm=0, TTwlpup=0, TTwlpum=0;
  if ( use_TTwl ) { GetYield(Lumi,idm,"TTwl",CoupB,CoupBjp,CoupBjm, TTwl, eTTwl, TTwlp, TTwlm, TTwlpup, TTwlpum); }   
  double TTwqq=0, eTTwqq=0, TTwqqp=0, TTwqqm=0, TTwqqpup=0, TTwqqpum=0;
  if ( use_TTwqq ) { GetYield(Lumi,idm,"TTwqq",CoupB,CoupBjp,CoupBjm, TTwqq, eTTwqq, TTwqqp, TTwqqm, TTwqqpup, TTwqqpum); }   
  double TTzl=0, eTTzl=0, TTzlp=0, TTzlm=0, TTzlpup=0, TTzlpum=0;
  if ( use_TTzl ) { GetYield(Lumi,idm,"TTzl",CoupB,CoupBjp,CoupBjm, TTzl, eTTzl, TTzlp, TTzlm, TTzlpup, TTzlpum); }   
  double TTzl1m10=0, eTTzl1m10=0, TTzl1m10p=0, TTzl1m10m=0, TTzl1m10pup=0, TTzl1m10pum=0;
  if ( use_TTzl1m10 ) { GetYield(Lumi,idm,"TTzl1m10",CoupB,CoupBjp,CoupBjm, TTzl1m10, eTTzl1m10, TTzl1m10p, TTzl1m10m, TTzl1m10pup, TTzl1m10pum); }   
  double TTzqq=0, eTTzqq=0, TTzqqp=0, TTzqqm=0, TTzqqpup=0, TTzqqpum=0;
  if ( use_TTzqq ) { GetYield(Lumi,idm,"TTzqq",CoupB,CoupBjp,CoupBjm, TTzqq, eTTzqq, TTzqqp, TTzqqm, TTzqqpup, TTzqqpum); }   

  // Data
  TH1F *mdata = new TH1F("mdata","Data",200,-100,100);
  if(use_Data){
    TH1F *mdem = new TH1F("mdem","Data",200,-100,100);
    TChain l1("bdttree");
    l1.Draw("nVert>>mdem",trgD+Coup);
    mdata->Add(mdem);
    mdata->Sumw2();
    mdata->SetMarkerStyle(8);
  }
  // Backgrounds
  double WjT = W1j+W2j+W3j+W4j+W5j+W6j+W7j;
  double ZiT = Zinv1+Zinv2+Zinv3+Zinv4+Zinv5+Zinv6+Zinv7;
  double VVT = WW+WZ+ZZ;
  double VVTpup = WWpup+WZpup+ZZpup;
  double VVTpum = WWpum+WZpum+ZZpum;
  double VVTp = WWp+WZp+ZZp;
  double VVTm = WWm+WZm+ZZm;
  double STT = sTw+sTtC+sTBw+sTBtC;
  double STTpup = sTwpup+sTtCpup+sTBwpup+sTBtCpup;
  double STTpum = sTwpum+sTtCpum+sTBwpum+sTBtCpum;
  double STTp = sTwp+sTtCp+sTBwp+sTBtCp;
  double STTm = sTwm+sTtCm+sTBwm+sTBtCm;
  double DYT = DY1+DY2+DY3+DY4+DY5+DY6+DY7+DY8+DY9+DY10+DY11;
  double sISR1pDYT = sISR1pDY1+sISR1pDY2+sISR1pDY3+sISR1pDY4+sISR1pDY5+sISR1pDY6+sISR1pDY7+sISR1pDY8+sISR1pDY9+sISR1pDY10+sISR1pDY11;
  double sISR2pDYT = sISR1pDY2+sISR2pDY2+sISR2pDY3+sISR2pDY4+sISR2pDY5+sISR2pDY6+sISR2pDY7+sISR2pDY8+sISR2pDY9+sISR2pDY10+sISR2pDY11;
  double sISR3pDYT = sISR1pDY3+sISR3pDY2+sISR3pDY3+sISR3pDY4+sISR3pDY5+sISR3pDY6+sISR3pDY7+sISR3pDY8+sISR3pDY9+sISR3pDY10+sISR3pDY11;
  double sISR4pDYT = sISR1pDY4+sISR4pDY2+sISR4pDY3+sISR4pDY4+sISR4pDY5+sISR4pDY6+sISR4pDY7+sISR4pDY8+sISR4pDY9+sISR4pDY10+sISR4pDY11;
  double sISR5pDYT = sISR1pDY5+sISR5pDY2+sISR5pDY3+sISR5pDY4+sISR5pDY5+sISR5pDY6+sISR5pDY7+sISR5pDY8+sISR5pDY9+sISR5pDY10+sISR5pDY11;
  double sISR6pDYT = sISR1pDY6+sISR6pDY2+sISR6pDY3+sISR6pDY4+sISR6pDY5+sISR6pDY6+sISR6pDY7+sISR6pDY8+sISR6pDY9+sISR6pDY10+sISR6pDY11;
  double sISR7pDYT = sISR1pDY7+sISR7pDY2+sISR7pDY3+sISR7pDY4+sISR7pDY5+sISR7pDY6+sISR7pDY7+sISR7pDY8+sISR7pDY9+sISR7pDY10+sISR7pDY11;
  double sISR1mDYT = sISR1mDY1+sISR1mDY2+sISR1mDY3+sISR1mDY4+sISR1mDY5+sISR1mDY6+sISR1mDY7+sISR1mDY8+sISR1mDY9+sISR1mDY10+sISR1mDY11;
  double sISR2mDYT = sISR1mDY2+sISR2mDY2+sISR2mDY3+sISR2mDY4+sISR2mDY5+sISR2mDY6+sISR2mDY7+sISR2mDY8+sISR2mDY9+sISR2mDY10+sISR2mDY11;
  double sISR3mDYT = sISR1mDY3+sISR3mDY2+sISR3mDY3+sISR3mDY4+sISR3mDY5+sISR3mDY6+sISR3mDY7+sISR3mDY8+sISR3mDY9+sISR3mDY10+sISR3mDY11;
  double sISR4mDYT = sISR1mDY4+sISR4mDY2+sISR4mDY3+sISR4mDY4+sISR4mDY5+sISR4mDY6+sISR4mDY7+sISR4mDY8+sISR4mDY9+sISR4mDY10+sISR4mDY11;
  double sISR5mDYT = sISR1mDY5+sISR5mDY2+sISR5mDY3+sISR5mDY4+sISR5mDY5+sISR5mDY6+sISR5mDY7+sISR5mDY8+sISR5mDY9+sISR5mDY10+sISR5mDY11;
  double sISR6mDYT = sISR1mDY6+sISR6mDY2+sISR6mDY3+sISR6mDY4+sISR6mDY5+sISR6mDY6+sISR6mDY7+sISR6mDY8+sISR6mDY9+sISR6mDY10+sISR6mDY11;
  double sISR7mDYT = sISR1mDY7+sISR7mDY2+sISR7mDY3+sISR7mDY4+sISR7mDY5+sISR7mDY6+sISR7mDY7+sISR7mDY8+sISR7mDY9+sISR7mDY10+sISR7mDY11;
  double DYTpup = DY1pup+DY2pup+DY3pup+DY4pup+DY5pup+DY6pup+DY7pup+DY8pup+DY9pup+DY10pup+DY11pup;
  double DYTpum = DY1pum+DY2pum+DY3pum+DY4pum+DY5pum+DY6pum+DY7pum+DY8pum+DY9pum+DY10pum+DY11pum;
  double DYTp = DY1p+DY2p+DY3p+DY4p+DY5p+DY6p+DY7p+DY8p+DY9p+DY10p+DY11p;
  double DYTm = DY1m+DY2m+DY3m+DY4m+DY5m+DY6m+DY7m+DY8m+DY9m+DY10m+DY11m;
  double TTX = TTgj+TTwl+TTwqq+TTzl+TTzl1m10+TTzqq;
  double TTXpup = TTgjpup+TTwlpup+TTwqqpup+TTzlpup+TTzl1m10pup+TTzqqpup;
  double TTXpum = TTgjpum+TTwlpum+TTwqqpum+TTzlpum+TTzl1m10pum+TTzqqpum;
  double TTXp = TTgjp+TTwlp+TTwqqp+TTzlp+TTzl1m10p+TTzqqp;
  double TTXm = TTgjm+TTwlm+TTwqqm+TTzlm+TTzl1m10m+TTzqqm;
  double bckg  = DDwj + DDtt + DDfk + VVT + STT + DYT + TTX;

  double eBMC = sqrt(pow(eWW,2)+pow(eWZ,2)+pow(eZZ,2)+pow(esTw,2)+pow(esTtC,2)+pow(esTBw,2)+pow(esTBtC,2)+
		     pow(eDY1,2)+pow(eDY2,2)+pow(eDY3,2)+pow(eDY4,2)+pow(eDY5,2)+pow(eDY6,2)+pow(eDY7,2)+pow(eDY8,2)+pow(eDY9,2)+pow(eDY10,2)+pow(eDY11,2)+
		     pow(eTTgj,2)+pow(eTTwl,2)+pow(eTTwqq,2)+pow(eTTzl,2)+pow(eTTzl1m10,2)+pow(eTTzqq,2));
  cout << "Bckg(Other MC) = " << (VVT + STT + DYT + TTX) << " +- " << eBMC << std::endl;

  // Signal point to cover
  string Spt = " ";

  for (int i=0; i<dmlines.size()-1; i++){ 

    Spt = dmlines[i];
    double SIGNAL=0, eSIGNAL=0, SIGNALp=0, SIGNALm=0, SIGNALpup=0, SIGNALpum=0;
    double sISR1pS=0, sISR2pS=0, sISR3pS=0, sISR4pS=0, sISR5pS=0, sISR6pS=0, sISR7pS=0;
    double sISR1mS=0, sISR2mS=0, sISR3mS=0, sISR4mS=0, sISR5mS=0, sISR6mS=0, sISR7mS=0;
    double SIGNALidp=0, SIGNALidm=0, SIGNALisop=0, SIGNALisom=0;
    if ( use_SIGNAL ) { 
      GetYield(Lumi,idm,Spt,Coup,Coupjp,Coupjm, SIGNAL, eSIGNAL, SIGNALp, SIGNALm, SIGNALpup, SIGNALpum); 
      GetSysISR(Lumi,idm,Spt,Coup,
		sISR1pS,sISR2pS,sISR3pS,sISR4pS,sISR5pS,sISR6pS,sISR7pS,
		sISR1mS,sISR2mS,sISR3mS,sISR4mS,sISR5mS,sISR6mS,sISR7mS); 
      GetSysIDISO(Lumi,idm,Spt,Coup, SIGNALidp, SIGNALidm, SIGNALisop, SIGNALisom); 
    }
    cout << "Signal point: " << Spt << std::endl;

    // Create & write the DataCard file
    ofstream outFile;
    outFile.open(("DataCards/" + Spt + ".txt").c_str());
    outFile << "imax 1 number of channels" << std::endl;
    outFile << "jmax * number of backgrounds" << std::endl;
    outFile << "jmax * number of backgrounds" << std::endl;
    outFile << "----------------------------" << std::endl;
    outFile << "bin " << Spt << "_SR" << std::endl;
    outFile << "observation " << bckg << std::endl;
    outFile << "----------------------------" << std::endl;
    outFile << "bin " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << Spt << "_SR " << std::endl;
    outFile << "process Sgn Wj tt Fake VV ST DY TTX" << std::endl;
    outFile << "process 0 1 2 3 4 5 6 7" << std::endl;
    outFile << "rate       " << SIGNAL << " " << DDwj << " " << DDtt << " " << DDfk << " " << VVT << " " << STT << " " << DYT << " " << TTX << std::endl;
    outFile << "SSta   lnN " << 1.+(eSIGNAL/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "Wj     lnN - " << 1.+(eDDwj/DDwj) << " - - - - - -" << std::endl;
    outFile << "tt     lnN - - " << 1.+(eDDtt/DDtt) << " - - - - -" << std::endl;
    outFile << "Fake   lnN - - - " << 1.+(eDDfk/DDfk) << " - - - -" << std::endl;
    outFile << "VV     lnN - - - - 1.5 - - -" << std::endl;
    outFile << "ST     lnN - - - - - 1.5 - -" << std::endl;
    outFile << "DY     lnN - - - - - - 1.5 -" << std::endl;
    outFile << "TTX    lnN - - - - - - - 1.5" << std::endl;
    outFile << "Lum    lnN 1.025 - - - 1.025 1.025 1.025 1.025" << std::endl;
    outFile << "Trg    lnN 1.010 - - - 1.010 1.010 1.010 1.010" << std::endl;
    outFile << "PU     lnN " << (SIGNALpup/SIGNAL) << "/" << (SIGNALpum/SIGNAL) << " - - - " << (VVTpup/VVT) << "/" << (VVTpum/VVT) 
	    << " " << (STTpup/STT) << "/" << (STTpum/STT) << " " << (DYTpup/DYT) << "/" << (DYTpum/DYT) << " " << (TTXpup/TTX) 
	    << "/" << (TTXpum/TTX) << std::endl;
    //    outFile << "JESR   lnN 1.04 - - - 1.04 1.04 1.04 1.04" << std::endl;
    outFile << "JESR   lnN " << (SIGNALp/SIGNAL) << "/" << (SIGNALp/SIGNAL) << " - - - " << 
      (VVTp/VVT) << "/" << (VVTm/VVT) << " " << (STTp/STT) << "/" << (STTm/STT) << " " << (DYTp/DYT) << "/" << (DYTm/DYT) << 
      (TTXp/TTX) << "/" << (TTXm/TTX) << std::endl;
    //    outFile << "ID     lnN " << (SIGNALidp/SIGNAL)  << "/" << (SIGNALidm/SIGNAL) <<  " - - - 1.03 1.03 1.03 1.03" << std::endl;
    //    outFile << "ISO    lnN " << (SIGNALisop/SIGNAL) << "/" << (SIGNALisom/SIGNAL) << " - - - 1.02 1.02 1.02 1.02" << std::endl;
    outFile << "ID     lnN 1.03 - - - 1.03 1.03 1.03 1.03" << std::endl;
    outFile << "ISO    lnN 1.02 - - - 1.02 1.02 1.02 1.02" << std::endl;
    outFile << "SISR1  lnN " << (sISR1pS/SIGNAL) << "/" << (sISR1mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "SISR2  lnN " << (sISR2pS/SIGNAL) << "/" << (sISR2mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "SISR3  lnN " << (sISR3pS/SIGNAL) << "/" << (sISR3mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "SISR4  lnN " << (sISR4pS/SIGNAL) << "/" << (sISR4mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "SISR5  lnN " << (sISR5pS/SIGNAL) << "/" << (sISR5mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "SISR6  lnN " << (sISR6pS/SIGNAL) << "/" << (sISR6mS/SIGNAL) << " - - - - - - -" << std::endl;
    outFile << "DYISR1 lnN " << "- - - - - - " << (sISR1pDYT/DYT) << "/" << (sISR1mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR2 lnN " << "- - - - - - " << (sISR2pDYT/DYT) << "/" << (sISR2mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR3 lnN " << "- - - - - - " << (sISR3pDYT/DYT) << "/" << (sISR3mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR4 lnN " << "- - - - - - " << (sISR4pDYT/DYT) << "/" << (sISR4mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR5 lnN " << "- - - - - - " << (sISR5pDYT/DYT) << "/" << (sISR5mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR6 lnN " << "- - - - - - " << (sISR6pDYT/DYT) << "/" << (sISR6mDYT/DYT) << " -" << std::endl;
    outFile << "DYISR7 lnN " << "- - - - - - " << (sISR7pDYT/DYT) << "/" << (sISR7mDYT/DYT) << " -" << std::endl;
    outFile << "" << std::endl;

  }

  delete mdata;

}


void GetYield(
	      double lumi,     // Input: Luminosity
	      string dm,
	      string pc,       // Input: Name of the channel
	      TCut Coupure,    // Input: Analysis cut
	      TCut Coupurejp,  // Input: +1sigma(JES) analysis cut
	      TCut Coupurejm,  // Input: -1sigma(JES) analysis cut
	      double& Yield,   // Output: Weighted yield of the channel
	      double& eYield,  // Output: Error on the weighted yield
	      double& YieldP,  // Output: Weighted yield of the channel +1sigma(JES)
	      double& YieldM,  // Output: Weighted yield of the channel -1sigma(JES)
	      double& YPUP,    // Output: Weighted yield of the channel +1sigma(PU)
	      double& YPUM     // Output: Weighted yield of the channel -1sigma(PU)
	      ){ 

  TChain ProcFile("bdttree");
  // Directory of roottuples
  //  string BaseDir = "/lstore/cms/bargassa/Stop4body/SET2102_" + dm + "/";
  //  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt" + dm + "/";
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_SysVar/";
  string FFile = " ";

  // File to read
  //  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + "_bdt.root";  }
  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + ".root";  }
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  //  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200_bdt.root";  }
  //  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400_bdt.root";  }
  //  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600_bdt.root";  }
  //  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800_bdt.root";  }
  //  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200_bdt.root";  }
  //  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500_bdt.root";  }
  //  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf_bdt.root";  }
  //  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200_bdt.root";  }
  //  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400_bdt.root";  }
  //  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600_bdt.root";  }
  //  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf_bdt.root";  }
  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200.root";  }
  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400.root";  }
  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600.root";  }
  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800.root";  }
  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200.root";  }
  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500.root";  }
  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf.root";  }
  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200.root";  }
  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400.root";  }
  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600.root";  }
  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf.root";  }
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
  //  if (pc=="WW"){     FFile = BaseDir + "WW_bdt.root";  }
  //  if (pc=="WZ"){     FFile = BaseDir + "WZ_bdt.root";  }
  //  if (pc=="ZZ"){     FFile = BaseDir + "ZZ_bdt.root";  }
  if (pc=="WW"){     FFile = BaseDir + "WW.root";  }
  if (pc=="WZ"){     FFile = BaseDir + "WZ.root";  }
  if (pc=="ZZ"){     FFile = BaseDir + "ZZ.root";  }
  if (pc=="SingleToptW"){     FFile = BaseDir + "T_tWch_ext_bdt.root";  }
  if (pc=="SingleToptChannel"){     FFile = BaseDir + "T_tch_powheg_bdt.root";  }
  if (pc=="SingleTopBtW"){     FFile = BaseDir + "TBar_tWch_ext_bdt.root";  }
  if (pc=="SingleTopBtChannel"){     FFile = BaseDir + "TBar_tch_powheg_bdt.root";  }
  if (pc=="TTgj"){     FFile = BaseDir + "TTGJets_bdt.root";  }
  if (pc=="TTwl"){     FFile = BaseDir + "TTWToLNu_bdt.root";  }
  if (pc=="TTwqq"){     FFile = BaseDir + "TTWToQQ_bdt.root";  }
  if (pc=="TTzl"){     FFile = BaseDir + "TTZToLLNuNu_bdt.root";  }
  if (pc=="TTzl1m10"){     FFile = BaseDir + "TTZToLLNuNu_m1to10_bdt.root";  }
  if (pc=="TTzqq"){     FFile = BaseDir + "TTZToQQ_bdt.root";  }
  ProcFile.Add((FFile).c_str());
  
  // Get the weighted numerator, where division by denominator has been done
  TH1D *nbcut = new TH1D("nbcut","nbcut",200,-100,100);
  nbcut->Sumw2();
  TH1D *nbcutjp = new TH1D("nbcutjp","nbcutjp",200,-100,100);
  TH1D *nbcutjm = new TH1D("nbcutjm","nbcutjm",200,-100,100);
  TH1D *nbcutpup = new TH1D("nbcutpup","nbcutpup",200,-100,100);
  TH1D *nbcutpum = new TH1D("nbcutpum","nbcutpum",200,-100,100);
  ProcFile.Draw("Njet>>nbcut","weight*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutjp","weight*splitFactor"*Coupurejp);
  ProcFile.Draw("Njet>>nbcutjm","weight*splitFactor"*Coupurejm);
  ProcFile.Draw("Njet>>nbcutpup","weight_PU_Up*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutpum","weight_PU_Down*splitFactor"*Coupure);
  //
  double NumWeight = 0., NumErr2 = 0.;
  double NumWeightP = 0., NumWeightM = 0.;
  double NumWeightPUP = 0., NumWeightPUM = 0.;
  for(int i=0; i<=(nbcut->GetNbinsX()+1); i++){ 
    NumWeight += nbcut->GetBinContent(i); 
    double NumErr = nbcut->GetBinError(i);
    NumErr2 += (NumErr*NumErr);
    NumWeightPUP += nbcutpup->GetBinError(i);
    NumWeightPUM += nbcutpum->GetBinError(i);
  }
  for(int i=0; i<=(nbcutjp->GetNbinsX()+1); i++){ NumWeightP += nbcutjp->GetBinContent(i); }
  for(int i=0; i<=(nbcutjm->GetNbinsX()+1); i++){ NumWeightM += nbcutjm->GetBinContent(i); }
  // Calculate yields 
  Yield  = lumi*NumWeight;
  eYield = lumi*sqrt(NumErr2);
  YieldP = lumi*NumWeightP;
  YieldM = lumi*NumWeightM;
  YPUP   = lumi*NumWeightPUP;
  YPUM   = lumi*NumWeightPUM;

  // Delete histograms
  delete nbcut;
  delete nbcutjp;
  delete nbcutjm;
  delete nbcutpup;
  delete nbcutpum;

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
  //  string BaseDir = "/lstore/cms/bargassa/Stop4body/SET2102_" + dm + "/";
  //  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt" + dm + "/";
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_SysVar/";
  string FFile = " ";

  // File to read
  // File to read
  //  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + "_bdt.root";  }
  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + ".root";  }
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  //  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200_bdt.root";  }
  //  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400_bdt.root";  }
  //  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600_bdt.root";  }
  //  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800_bdt.root";  }
  //  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200_bdt.root";  }
  //  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500_bdt.root";  }
  //  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf_bdt.root";  }
  //  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200_bdt.root";  }
  //  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400_bdt.root";  }
  //  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600_bdt.root";  }
  //  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf_bdt.root";  }
  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200.root";  }
  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400.root";  }
  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600.root";  }
  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800.root";  }
  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200.root";  }
  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500.root";  }
  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf.root";  }
  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200.root";  }
  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400.root";  }
  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600.root";  }
  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf.root";  }
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
  if ((pc=="DY1")||(pc=="DY2")||(pc=="DY3")||(pc=="DY4")||(pc=="DY5")||
      (pc=="DY6")||(pc=="DY7")||(pc=="DY8")||(pc=="DY9")||(pc=="DY10")||(pc=="DY11")||
      (pc=="W1Jets")||(pc=="W2Jets")||(pc=="W3Jets")||(pc=="W4Jets")
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
  else if ((pc==pc)||(pc=="TTbar")){
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


void GetSysIDISO(
		 double lumi,     // Input: Luminosity
		 string dm,
		 string pc,       // Input: Name of the channel
		 TCut Coupure,    // Input: Analysis cut
		 double& YIDP,    // Output: Weighted yield of the channel +1sigma(ID)
		 double& YIDM,    // Output: Weighted yield of the channel -1sigma(ID)
		 double& YISOP,   // Output: Weighted yield of the channel +1sigma(ISO)
		 double& YISOM    // Output: Weighted yield of the channel -1sigma(ISO)
		 ){ 

  TChain ProcFile("bdttree");
  // Directory of roottuples
  //  string BaseDir = "/lstore/cms/bargassa/Stop4body/SET2102_" + dm + "/";
  //  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_test_bdt" + dm + "/";
  string BaseDir = "/lstore/cms/cbeiraod/Stop4Body/nTuples_v2017-06-05_SysVar/";
  string FFile = " ";

  // File to read
  //  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + "_bdt.root";  }
  if (pc==pc){     FFile = BaseDir + "T2DegStop_" + pc + ".root";  }
  if (pc=="TTbar"){     FFile = BaseDir + "TT_pow_bdt.root";  }
  //  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200_bdt.root";  }
  //  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400_bdt.root";  }
  //  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600_bdt.root";  }
  //  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800_bdt.root";  }
  //  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200_bdt.root";  }
  //  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500_bdt.root";  }
  //  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf_bdt.root";  }
  //  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200_bdt.root";  }
  //  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400_bdt.root";  }
  //  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600_bdt.root";  }
  //  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf_bdt.root";  }
  if (pc=="DY1"){    FFile = BaseDir + "DYJetsToLL_M50_HT100to200.root";  }
  if (pc=="DY2"){    FFile = BaseDir + "DYJetsToLL_M50_HT200to400.root";  }
  if (pc=="DY3"){     FFile = BaseDir + "DYJetsToLL_M50_HT400to600.root";  }
  if (pc=="DY4"){    FFile = BaseDir + "DYJetsToLL_M50_HT600to800.root";  }
  if (pc=="DY5"){    FFile = BaseDir + "DYJetsToLL_M50_HT800to1200.root";  }
  if (pc=="DY6"){     FFile = BaseDir + "DYJetsToLL_M50_HT1200to2500.root";  }
  if (pc=="DY7"){     FFile = BaseDir + "DYJetsToLL_M50_HT2500toInf.root";  }
  if (pc=="DY8"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT100to200.root";  }
  if (pc=="DY9"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT200to400.root";  }
  if (pc=="DY10"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT400to600.root";  }
  if (pc=="DY11"){     FFile = BaseDir + "DYJetsToLL_M5to50_HT600toInf.root";  }
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
  //  if (pc=="WW"){     FFile = BaseDir + "WW_bdt.root";  }
  //  if (pc=="WZ"){     FFile = BaseDir + "WZ_bdt.root";  }
  //  if (pc=="ZZ"){     FFile = BaseDir + "ZZ_bdt.root";  }
  if (pc=="WW"){     FFile = BaseDir + "WW.root";  }
  if (pc=="WZ"){     FFile = BaseDir + "WZ.root";  }
  if (pc=="ZZ"){     FFile = BaseDir + "ZZ.root";  }
  if (pc=="SingleToptW"){     FFile = BaseDir + "T_tWch_ext_bdt.root";  }
  if (pc=="SingleToptChannel"){     FFile = BaseDir + "T_tch_powheg_bdt.root";  }
  if (pc=="SingleTopBtW"){     FFile = BaseDir + "TBar_tWch_ext_bdt.root";  }
  if (pc=="SingleTopBtChannel"){     FFile = BaseDir + "TBar_tch_powheg_bdt.root";  }
  if (pc=="TTgj"){     FFile = BaseDir + "TTGJets_bdt.root";  }
  if (pc=="TTwl"){     FFile = BaseDir + "TTWToLNu_bdt.root";  }
  if (pc=="TTwqq"){     FFile = BaseDir + "TTWToQQ_bdt.root";  }
  if (pc=="TTzl"){     FFile = BaseDir + "TTZToLLNuNu_bdt.root";  }
  if (pc=="TTzl1m10"){     FFile = BaseDir + "TTZToLLNuNu_m1to10_bdt.root";  }
  if (pc=="TTzqq"){     FFile = BaseDir + "TTZToQQ_bdt.root";  }
  ProcFile.Add((FFile).c_str());
  
  // Get the weighted numerator, where division by denominator has been done
  TH1D *nbcutidp = new TH1D("nbcutidp","nbcutidp",200,-100,100);
  TH1D *nbcutidm = new TH1D("nbcutidm","nbcutidm",200,-100,100);
  TH1D *nbcutisop = new TH1D("nbcutisop","nbcutisop",200,-100,100);
  TH1D *nbcutisom = new TH1D("nbcutisom","nbcutisom",200,-100,100);
  ProcFile.Draw("Njet>>nbcutidp","weight_LeptonIDSF_AltCorr_Up*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutidm","weight_LeptonIDSF_AltCorr_Down*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutisop","weight_LeptonISOSF_AltCorr_Up*splitFactor"*Coupure);
  ProcFile.Draw("Njet>>nbcutisom","weight_LeptonISOSF_AltCorr_Down*splitFactor"*Coupure);
  //
  double NumWeightIDP = 0., NumWeightIDM = 0., NumWeightISOP = 0., NumWeightISOM = 0.;
  for(int i=0; i<=(nbcutidp->GetNbinsX()+1); i++){ 
    NumWeightIDP += nbcutidp->GetBinContent(i);
    NumWeightIDM += nbcutidm->GetBinContent(i);
    NumWeightISOP += nbcutisop->GetBinContent(i);
    NumWeightISOM += nbcutisom->GetBinContent(i);
  }
  // Calculate yields 
  YIDP   = lumi*NumWeightIDP;
  YIDM   = lumi*NumWeightIDM;
  YISOP  = lumi*NumWeightISOP;
  YISOM  = lumi*NumWeightISOM;

  // Delete histograms
  delete nbcutidp;
  delete nbcutidm;
  delete nbcutisop;
  delete nbcutisom;

}


double Max(double a, double b){
  if (a>b) return a;
  else return b;
}

double AbsV(double a){
  return sqrt(pow(a,2));
}
