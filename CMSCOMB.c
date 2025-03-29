//CMS COMBINATIONS: By Thabo Pilusa



#include "RooPlot.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h" 
#include "RooSimultaneous.h"	
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace RooStats;
using namespace RooFit;

void CMSCOMB()
 {
  

   
    // S - 2017 VBF
    // E - 2018 VBF
   
   
// read from the 2 fileS and create a ROOT tree
   TTree treeS("treeS","treeS");
   int nevtS= treeS.ReadFile("2017CMSVBF.txt","x:y");
   std::cout << "Read " << nevtS << " from the file " << std::endl;


   TTree treeE("treeE","treeE");
   int nevtE = treeE.ReadFile("2018CMSVBF.txt","x:y");
   std::cout << "Read " << nevtE << " from the file " << std::endl;

// Making the RooFit Models for both years
   RooWorkspace w("w");
  
   
// For 2017
  w.factory("x[65,120]");//The invariant mass range
  w.factory("y[8.03903,140]");//The invariant mass range
  w.factory("nbackgroundS[0,-100,100]");
   w.var("nbackgroundS")->setVal(nevtS);
   w.var("nbackgroundS")->setMin(0.1*nevtS);
   w.var("nbackgroundS")->setMax(100*nevtS);
  
//Exponential model



  /* w.factory("S[-0.0015, -100, 100]");
   w.factory("Exponential::bmodelS(x, S)");
   RooAbsPdf * bmodelS = w.pdf("bmodelS");*/
   
  w.factory("S1[-100, 100]");
  w.factory("S2[-100, 100]");
  w.factory("S3[-100, 100]");
   
  w.factory("expr::s('S1*x + S2*x*x + S3*x*x*x',x,S1, S2,S3)");
  w.factory("Exponential::bmodelS(s,1)");
  RooAbsPdf * bmodelS = w.pdf("bmodelS"); 
   
   
  
//Define the signal model
 w.factory("nsignalS[0, -100 ,100]");
//😄️Now define the mass, width(resolution/specs) and use them to make a the smodel pdf

 
 w.factory("massS[95, 95, 95]");
 w.factory("widthS[1.33, 1.33, 1.33]");
 w.factory("Gaussian::smodelS(x, massS,widthS)");
 RooAbsPdf * smodelS = w.pdf("smodelS");

// the sum of the background and the signal model

 w.factory("SUM::modelS(nbackgroundS*bmodelS, nsignalS * smodelS)");
 RooAbsPdf * modelS = w.pdf("modelS");
 
 
//............DO FOR 2018 2018..........................................................................................................................

  *w.factory("nbackgroundE[0,-100,100]");
  w.var("nbackgroundE")->setVal(nevtE);
  w.var("nbackgroundE")->setMin(0.1*nevtE);
  w.var("nbackgroundE")->setMax(100*nevtE);
   
   
  w.factory("E1[-100, 100]");
  w.factory("E2[-100, 100]");
  w.factory("E3[-100, 100]");
   
  w.factory("expr::e('E1*x + E2*x*x + E3*x*x*x',x,E1, E2,E3)");
  w.factory("Exponential::bmodelE(e,1)");
  RooAbsPdf * bmodelE = w.pdf("bmodelE"); 
 
 
   /*w.factory("E[0, -100, 1000]");
   w.factory("Exponential::bmodelE(x, E)");
   RooAbsPdf * bmodelE = w.pdf("bmodelE");*/
      
 
 
 w.factory("nsignalE[0, -100, 100]");
 
 w.factory("massE[95, 95, 95]");
 w.factory("widthE[1.21, 1.21, 1.21]");
 w.factory("Gaussian::smodelE(x, massE, widthE)");
 RooAbsPdf * smodelE = w.pdf("smodelE");
 
 w.factory("SUM::modelE(nbackgroundE * bmodelE, nsignalE * smodelE)");
 RooAbsPdf * modelE = w.pdf("modelE");
 
//Now creating RooDataSet for both cases..............................................................................................................................
 RooDataSet dataS("dataS","dataS",RooArgSet(*w.var("x"),*w.var("y")),WeightVar(*w.var("y")),Import(treeS));
 RooDataSet dataE("dataE", "dataE", RooArgSet(*w.var("x"), *w.var("y")),WeightVar(*w.var("y")),Import(treeE));


// Join both samples.......................................................................................................................................................
//Create categories to differentiate between both samples

//Now we make the mixed, combined dataset
RooCategory sample("sample", "Index Category");
sample.defineType("VBFS");
sample.defineType("VBFE");

RooDataSet mixedset("mixedset", "mixedset", RooArgSet(*w.var("x"), *w.var("y")), RooFit::Index(sample), RooFit::Import("VBFS", dataS), RooFit::Import("VBFE", dataE), RooFit::WeightVar(*w.var("y")));

 //RooDataSet mixedset("mixedset","mixedset", RooArgSet(*w.var("x"),*w.var("y")),WeightVar(*w.var("y")),Import("VBFS", dataS),Import("VBFE", dataE));

//Construct a Simultaneous pdf using category sample as index
 RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);


//Associate model with physics state and model_ctl with the control state😁️
 simPdf.addPdf(*modelS, "VBFS");
 simPdf.addPdf(*modelE, "VBFE");




//Do a sumultenous fit of model to data and model_ctl to data_ctl
  RooFitResult * r = simPdf.fitTo(mixedset, SumW2Error(true),Minimizer("Minuit2"),Save(true), Offset(true));
  r->Print();
  
//Plot model the data

  RooPlot * frame1 = w.var("x")->frame(Bins(55),Title("VBF 2017 data"));


  dataS.plotOn(frame1,DataError(RooAbsData::Poisson),MarkerColor(kBlack),LineColor(kBlack),XErrorSize(0),Name("Data17"));  
   
  RooPlot * frame2 = w.var("x")->frame(Bins(55),Title("VBF 2018 data"));


  dataE.plotOn(frame2,DataError(RooAbsData::Poisson),MarkerColor(kBlack),LineColor(kBlack),XErrorSize(0),Name("Data18"));  
  
  
   //.........................................................
   
   
   
   modelS->plotOn(frame1,Range(65, 120),Name("SignalbackgroundS"));
   modelS->plotOn(frame1, Components("bmodelS"),LineStyle(kDashed),Name("backgroundS"));
  RooCurve *bkg = (RooCurve*)frame1->getObject(2);

   modelS->plotOn(frame1, Range(65, 120),Name("BSMsignal"),Components("smodelS"),LineColor(kGreen));
   RooCurve *sig= (RooCurve*)frame1->getObject(3);
   
   
   
   modelE->plotOn(frame2,Range(65, 120),Name("SignalbackgroundE"));
   modelE->plotOn(frame2, Components("bmodelE"),LineStyle(kDashed),Name("backgroundE"));
  RooCurve *bkg1 = (RooCurve*)frame2->getObject(2);

   modelE->plotOn(frame2, Range(65, 120),Name("BSMsignal"),Components("smodelE"),LineColor(kGreen));
   RooCurve *sig1 = (RooCurve*)frame2->getObject(3);
   
   
   
   
   

/*
//   plot all data tagged as physics  sample     
mixedset.plotOn(frame1,  DataError(RooAbsData::Poisson),XErrorSize(0),Name("sample==sample::VBFS"));       
         
//plot  "physics" slice of simultaneous pdf
//fill row_hist_project the sample index category with data using ProjWData


//simPdf.plotOn(frame1, Slice(sample, "VBFS"),ProjWData(sample,mixedset)) ; 
dataS.plotOn(frame1, Slice(sample, "VBFS"),ProjWData(sample,mixedset), Range(65,120),Name("SignalBackground"));
dataS.plotOn(frame1, Components("bmodelS"),LineStyle(kDashed),Name("backgroundS"));
RooCurve *bkg = (RooCurve*)frame1->getObject(2);
dataS.plotOn(frame1, Slice(sample, "VBFS"),Range(65, 120), Name("BSM signalS"),Components("smodelS"),LineColor(kGreen));
RooCurve *sig = (RooCurve*)frame1->getObject(3);

// Now do the same for the control sample slice😪️
RooPlot * frame2 = w.var("x")->frame(Bins(55),Title("VBF 2018 "));

dataE.plotOn(frame2,  DataError(RooAbsData::Poisson),XErrorSize(0),Name("sample==sample::VBFE"));


dataE.plotOn(frame2, Slice(sample, "VBFE"),ProjWData(sample,mixedset),Range(65,120));
dataE.plotOn(frame2, Components("bmodelE"),LineStyle(kDashed),Name("backgroundE"));
 RooCurve *bkg1 = (RooCurve*)frame2->getObject(2);
dataE.plotOn(frame2, Slice(sample, "VBFE"),Range(65, 120), Name("BSM signalE"),Components("smodelE"),LineColor(kGreen));
 RooCurve *sig1 = (RooCurve*)frame2->getObject(3);
*/



   
  TCanvas* c = new TCanvas("Canvas", "Canvas",511,81,851,559);
  c->Divide(2);
  c->cd(1) ; 
  gPad->SetLeftMargin(0.25); 
  frame1->GetYaxis()->SetTitleOffset(1.4) ; 
  frame1 -> GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  frame1->Draw();
  c->cd(2) ; 
  gPad->SetLeftMargin(0.25); 
  frame2->GetYaxis()-> SetTitleOffset(1.4); 
  frame2 -> GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  frame2->Draw();

 
 
 
TLegend *leg1 = new TLegend(0.2,0.1,0.45,0.4);
leg1->SetFillColor(kWhite);
leg1->SetLineColor(kWhite);
leg1->AddEntry("Data18","Data", "EP");
leg1->AddEntry("SignalbackgroundE","Signal + background","LP");
leg1->AddEntry("backgroundE","Background only", "LP");
leg1->AddEntry("BSMsignal","BSM Signal", "LP");
leg1->Draw();


 
 
 RooRealVar* par0_fitresultS = (RooRealVar*) r->floatParsFinal().find("nbackgroundS");  
 double fitted_nbkgS = par0_fitresultS->getVal();
 
 RooRealVar* par1_fitresultS = (RooRealVar*) r->floatParsFinal().find("nsignalS");
 double fitted_nsignalS = par1_fitresultS->getVal();
 
  RooRealVar* par0_fitresultE = (RooRealVar*) r->floatParsFinal().find("nbackgroundE");  
 double fitted_nbkgE = par0_fitresultE->getVal();
 
 RooRealVar* par1_fitresultE = (RooRealVar*) r->floatParsFinal().find("nsignalE");
 double fitted_nsignalE = par1_fitresultE->getVal();
 
 
  RooRealVar x = *w.var("x");
  
  
  x.setRange("signal",94.25,97.75);
  RooAbsReal *integral_sS = smodelS->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)), Range("signal"));
  RooAbsReal *integral_bS = bmodelS->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)), Range("signal")); 
  
  RooAbsReal *integral_sE = smodelE->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)), Range("signal"));
  RooAbsReal *integral_bE = bmodelE->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)), Range("signal")); 
   
   
   double sS = integral_sS->getVal() * fitted_nsignalS;
   double bS    = integral_bS->getVal() * fitted_nbkgS;
   double sE = integral_sE->getVal() * fitted_nsignalE;
   double bE    = integral_bE->getVal() * fitted_nbkgE;
   
	   cout<< " N signal VBF 2017 = "<< sS <<endl;
	   cout<<" N background VBF 2017 = "  << bS << endl;
	   
	   cout<< " N signal VBF 2018 = "<<sE << endl;
	   cout<<" N background VBF 2018= "  << bE << endl; 
	   
   double s_S = 41.5 * ( integral_sS->getVal() + fitted_nsignalS/41.5);
   double s_E = 54.4 * ( integral_sE->getVal() + fitted_nsignalE/54.4); 
   
   double sigS = sqrt(2*((bS+ s_S)*log(1 + s_S/bS) - s_S));
   double sigE = sqrt(2*((bE+ s_E)*log(1 + s_E/bE) - s_E));
   
     cout << "Significance VBF 2017 = " << sqrt(2*((bS+sS)*log(1 + sS/bS) -sS)) << endl;
     cout << "Significance VBF 2018 = " << sqrt(2*((bE+sE)*log(1 + sE/bE) -sE)) << endl;
  
    
  double Sig = sqrt( sigS*sigS + sigE*sigE);
  
  double pvalue = 2*( 1- ROOT::Math::normal_cdf(Sig));   
    cout<<"significance = "<< Sig << endl;
    cout<<" p_value = " << 2*(1 - ROOT::Math::normal_cdf(Sig))<< endl;
    double mass = 95;
    TString  label = "mass";
    ofstream Files;
    Files.open("results/"+label+".txt");
    Files << mass << "   " << pvalue << endl;
    Files.close();
    
    
//AND ALL GOOD THINGS COME TO AN END
// By : Thabo Pilusa 
 //...............
 
 
}
