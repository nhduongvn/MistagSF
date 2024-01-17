#include <iostream>
#include <TROOT.h>
#include <TArrow.h>
#include <TLatex.h>
#include "TH1F.h"
#include "TMath.h"
#include "TStyle.h"
#include "Bases/CMS_lumi.C"
#include "Bases/LowStatAndRebin.C"


void AddHist(TH1F* hTot, TH1F* h, TH1F* h1, float c1=1, float c2=1) { //do bin-by-bin adding
  for (int i = 0; i <= hTot->GetNbinsX() + 1; i++) {
    hTot->SetBinContent(i, h->GetBinContent(i)*c1 + h1->GetBinContent(i)*c2) ;
    float err = TMath::Sqrt(TMath::Power(h->GetBinError(i)*c1 ,2) + TMath::Power(h1->GetBinError(i)*c2 ,2)) ;
    hTot->SetBinError(i, err) ;
  }
}

int plot(TString tagger="DeepFlavour", TString cut="40", int MinEta = 1, int MaxEta = 25, int myNbins = 50, TString period = "2022",Double_t minAlpha = 0.1, Double_t maxAlpha = 0.3, bool ImprovePlot = false)
{
  vector<TString> PthatRanges;
  if(period=="2022EE") PthatRanges.push_back("/QCD_PT-15to30_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-30to50_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-50to80_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-80to120_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-120to170_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-170to300_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-300to470_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/");
  PthatRanges.push_back("/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/");

  TString EtaRange = "";
  TString HEtaRange = "";
  char buffer[50];

  if( MinEta > MaxEta )
  {
    std::cout << "MinEta should be smaller than MaxEta " << endl;
    return 0;
  }
  TH1F *htemp = new TH1F("htemp","htemp",25,0,2.5); // Temp hist used to reproduce used binning...
  if( MinEta < 1 || MinEta > htemp->GetNbinsX() )
  {
    std::cout << "MinEta out of range. Should be between 1 and " << htemp->GetNbinsX() << " But found to be " << MinEta << endl;
    return 0;
  }
  if( MaxEta < 1 || MaxEta > htemp->GetNbinsX() )
  {
    std::cout << "MaxEta out of range. Should be between 1 and " << htemp->GetNbinsX() << " But found to be " << MaxEta << endl;
    return 0;
  }

  float EtaMin = htemp->GetBinLowEdge(MinEta) ;
  float EtaMax = htemp->GetBinLowEdge(MaxEta) + htemp->GetBinWidth(MaxEta) ;

  cout << "\n Eta range is: " << MinEta << "  " << MaxEta << endl ;
  cout << "\n Eta range is 1: " << EtaMin << "  " << EtaMax << endl ;

  HEtaRange = "#eta #in [";
  sprintf(buffer,"%3.1f",EtaMin);
  EtaRange += buffer;
  HEtaRange += buffer;
  EtaRange += "_";
  HEtaRange += ", ";
  sprintf(buffer,"%3.1f",EtaMax);
  EtaRange += buffer;
  HEtaRange += buffer;
  HEtaRange += "]";

  std::cout << "Will Use EtaRange = " << EtaRange << endl;

  Double_t Cuts[5]= {-1000., -1000., -1000., -1000., -1000.};
  //Double_t TaggerS[3] = {"L", "M", "T" };
  const char* TaggerS[5] = {"L", "M", "T","XT","XXT" };
  TLatex *tl[5];
  TArrow *ta[5];
 
  TString RunFileName = "";
  int nbin = myNbins;
  double binmin = -25., binmax = 25.;
  TString tagFile = "";
  if(tagger == "DeepFlavour")
  {
    tagFile = tagger + "L";
    binmin = -1.0;
    binmax =  1.0;
    if (period == "2022") {
      Cuts[0] = 0.0583;
      Cuts[1] = 0.3086;
      Cuts[2] = 0.7183;
      Cuts[3] = 0.8111;
      Cuts[4] = 0.9512;
    }
    //FIXME
    //if (period == "2022EE") {
      
    //}

  }
  else if(tagger == "PNet")
  {
    tagFile = tagger + "L";
    binmin = -1.0;
    binmax =  1.0;
    if (period == "2022") {
      Cuts[0] =0.047; 
      Cuts[1] =0.245;
      Cuts[2] =0.6734;
      Cuts[3] =0.7862;
      Cuts[4] =0.961;
    }
  }
  else if(tagger == "ParT")
  {
    tagFile = tagger + "L";
    binmin = -1.0;
    binmax =  1.0;
    if (period == "2022") {
      Cuts[0] =0.0849; 
      Cuts[1] =0.4319;
      Cuts[2] =0.8482;
      Cuts[3] =0.9151;
      Cuts[4] =0.9874;
    }
  }
  else 
  {
    cout << "\n !!!Warning: No tagger name found" ;
  }
 
  if( myNbins != 100 && myNbins != 50 && myNbins != 20 && myNbins !=10 )
  {
    cout << "Bad Nbins, choose a pair integer such that 100 divided by it is integer" << endl;
    return 0;
  }
 // *****************************************************************************
 
  Int_t stati=0;
  Bool_t  fit=0;
  Bool_t logy=1;
 
  Double_t HisBinWidth = (binmax-binmin)/ ( (Double_t) nbin );
  sprintf(buffer,"%4.2f",HisBinWidth);
  TString tHBW = "";
  tHBW += buffer;
 
  float sfK0 = 1.4;
  float sfLA = 1.5;
 
  char* hname  = "hData_PTagger"; 
  char* hlight = "hLightFlav_PTagger"; char* hcflav = "hCFlav_PTagger"; char* hnoflav = "hNoFlav_PTagger"; 
 
 // *****************************************************************************
 
  //TCanvas *c1 = new TCanvas("c1", "plots",200,10,1000,810);
  TCanvas *c1 = new TCanvas("c1", "plots",200,10,720,820);
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
 
 // *****************************************************************************
 
 // TPaveLabel *p01 = new TPaveLabel(0.2,0.93,0.8,0.97,
 //                   "HLT_11B_Jet300: Data 4_1_4  versus  MC QCD 3_11_3","br");
 // p01->SetBorderSize(0);
 // p01->SetFillStyle(0);
 // p01->SetTextAlign(13);
 // p01->SetTextFont(42);
 // p01->SetTextSize(0.7);
 // p01->Draw();
 
 TPad *pad1 = new TPad("pad1","This is pad1",0.00,0.36,1.00,1.00,21);
 TPad *rap1 = new TPad("rap1","This is rap1",0.00,0.00,1.00,0.34,21);
 
 pad1->SetFillColor(0.05);
 pad1->SetBorderMode(0);
 pad1->SetFrameFillColor(10);
 pad1->Draw();
 pad1->SetLogy(1);
 pad1->SetTopMargin(0.09);
 pad1->SetBottomMargin(0.08);
 pad1->SetRightMargin(0.07);
 pad1->SetLeftMargin(0.17);

 rap1->SetFillColor(0);
 rap1->SetBorderMode(0);
 rap1->SetFrameFillColor(10);
 rap1->Draw();
 rap1->SetLogy(0);
 rap1->SetTopMargin(0.04);
 rap1->SetBottomMargin(0.3);
 rap1->SetRightMargin(0.065);
 rap1->SetLeftMargin(0.17);

  gStyle->SetOptDate(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleW(0.4);
  gStyle->SetTitleH(0.09);
// gStyle->SetTitleX(0); // Set the position of the title box
// gStyle->SetTitleY(0.985); // Set the position of the title box
// gStyle->SetTitleStyle(Style_t style = 1001);
// gStyle->SetTitleBorderSize(2);
  gStyle->SetOptStat(stati);
// gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

  if (fit) {
    gStyle->SetOptFit(111);
    gStyle->SetStatW(0.5);
    gStyle->SetStatH(0.2);
  } else {
    gStyle->SetOptFit(0);
    gStyle->SetStatW(0.4);
    gStyle->SetStatH(0.3);
  }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pad1->cd();

  std::cout << "\n I am here" << std::endl;
  
  TString Lumi = "";
  TFile *f1 = NULL;
  TFile *f2 = NULL;
  TString SavDir = "";

  TString DataFileRoot = "";
  TString MCFileRoot = "";

  if (period == "2022") {
    //Lumi = "36.5 fb^{-1} (13 TeV, 25ns)";
    Lumi = "?? fb^{-1} (2022, 13.6 TeV)";
    DataFileRoot = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTV_subjet_2022/Output_"+period+"_noPSw_fixPuppi_4/Data/JetHT_JetMET/JetTree_data_FatJets_Subjets_sampleList_all_NW_" ;
    MCFileRoot = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTV_subjet_2022/Output_"+period+"_noPSw_fixPuppi_4/MC/QCD/";
    //SavDir = "plot_" ;
    //SavDir += period;
    SavDir = "Plots_"+period+"_All_noPSw_fixPuppi_4/" ;
  }
  else if (period == "2022EE") {
    //Lumi = "36.5 fb^{-1} (13 TeV, 25ns)";
    Lumi = "?? fb^{-1} (2022EE, 13.6 TeV)";
    DataFileRoot = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTV_subjet_2022/Output_"+period+"_noPSw_fixPuppi_4/Data/JetHT_JetMET/JetTree_data_FatJets_Subjets_sampleList_all_NW_" ;
    MCFileRoot = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTV_subjet_2022/Output_"+period+"_noPSw_fixPuppi_4/MC/QCD/";
    //SavDir = "plot_" ;
    //SavDir += period;
    SavDir = "Plots_"+period+"_All_noPSw_fixPuppi_4/" ;
  }

  else {
    std::cout << "\n Unknow run period" << std::endl;
    return 0;
  }

  gSystem->mkdir(SavDir,kTRUE);



 TString htitle = tagger+" Discriminator"; 

 TH1F* h0= new TH1F("h0","", nbin, binmin, binmax);
 TH1F* h1= new TH1F("h1","", nbin, binmin, binmax);
 TH1F* h2= new TH1F("h2","", nbin, binmin, binmax);
 TH1F* hl= new TH1F("hl","", nbin, binmin, binmax);
 TH1F* hc= new TH1F("hc","", nbin, binmin, binmax);
 TH1F* hno= new TH1F("hno","", nbin, binmin, binmax);
 TH1F* hv0= new TH1F("hv0","", nbin, binmin, binmax);
 h1->Sumw2(); 
 h2->Sumw2(); 

 std::cout << "\n I am here" << std::endl;
 //f1 = new TFile(DataFileRoot + cut+"_"+tagger+"_WeightStudy.root");
 f1 = new TFile(DataFileRoot + cut+"_"+tagger+".root");
 TH1F* h10 = FindRebin(f1, hname, EtaRange+"_Data", MinEta,MaxEta); //hname = hData_Tagger
 std::cout << "\n I am here 1" << std::endl;


 h10->Rebin(100/myNbins);
 //h1->Add(h10,h1,1,0);
 AddHist(h1,h10,h1,1,0) ;

// PVH ci dessous
 //TFile *f2 = NULL;
 TH1F* h20 = NULL;
 TH1F* hl0 = NULL;
 TH1F* hno0 = NULL;
 TH1F* hc0 = NULL;

 for( unsigned long iPthat = 0; iPthat < PthatRanges.size(); ++iPthat )
 {
   //f2 = new TFile(MCFileRoot+PthatRanges.at(iPthat)+"JetTree_mc_FatJets_Subjets_Data2015D_16Dec2015_v1_NW_"+ cut+"_"+tagger+".root");
   //TString pthatrange = PthatRanges.at(iPthat) ;
   //pthatrange.ReplaceAll("/", "") ;
   //pthatrange.ReplaceAll("-", "_") ;
   //Get file which store MC plot
   //f2 = new TFile(MCFileRoot+pthatrange+"_TuneCP5_13TeV_pythia8/JetTree_mc_FatJets_Subjets_sampleList_all_NW_"+ cut+"_"+tagger+".root");
   f2 = new TFile(MCFileRoot+PthatRanges.at(iPthat)+"/JetTree_mc_FatJets_Subjets_sampleList_all_NW_"+ cut+"_"+tagger+".root");
   if( f2 == NULL )
   {
     std::cout << "Error opening " << MCFileRoot << PthatRanges.at(iPthat) << "JetTree_mc_FatJets_Subjets_sampleList_all_NW_" <<  cut << "_" << tagger << ".root" << endl;
     getchar();
     return 0;
   }// else std::cout << DirOfMCFile<<PthatRanges.at(iPthat)<<RootOfMCFile << "0_" << CutName << ".root opened " << endl;
   f2->cd();


   if( iPthat == 0 )
   {
     h20 = FindRebin(f2,hname,EtaRange+"_DataTagger",MinEta,MaxEta);
     hl0 = FindRebin(f2,hlight, EtaRange+"_LightTagger",MinEta,MaxEta);
     hno0 = FindRebin(f2,hnoflav,EtaRange+"_noTagger",MinEta,MaxEta);
     hc0 = FindRebin(f2,hcflav,EtaRange+"_cTagger",MinEta,MaxEta);
   }
   else
   {

     h20->Add(FindRebin(f2,hname,EtaRange+"_DataTagger"+PthatRanges.at(iPthat),MinEta,MaxEta));
     hl0->Add(FindRebin(f2,hlight, EtaRange+"_LightTagger"+PthatRanges.at(iPthat),MinEta,MaxEta));
     hno0->Add(FindRebin(f2,hnoflav, EtaRange+"_noTagger"+PthatRanges.at(iPthat),MinEta,MaxEta));
     hc0->Add(FindRebin(f2,hcflav, EtaRange+"_cTagger"+PthatRanges.at(iPthat),MinEta,MaxEta));
   }
   h20->SetDirectory(0);
   hl0->SetDirectory(0);
   hno0->SetDirectory(0);
   hc0->SetDirectory(0);

   f2->Close();
   delete f2;
 }
// PVH au dessus

 std::cout << "\n I am here 3" << std::endl;
       //g1->cd();
 //TH1F* h20 = FindRebin(g1, hname, EtaRange+"_Data", MinEta,MaxEta);
 h20->Rebin(100/myNbins);
 //h2->Add(h2,h20,1,1.); 
 AddHist(h2,h2,h20,1,1) ;
 float norm = h1->Integral(0,nbin+1)/h2->Integral(0,nbin+1);
 h2->Add(h2,h2,norm,0);

 //TH1F* hl0 = FindRebin(g1, hname, EtaRange+"_Data", MinEta,MaxEta);
 hl0->Rebin(100/myNbins);

 //h0->Add(hl0,h0,norm,0);
 AddHist(h0,hl0,h0,norm,0) ;
 
 //hl->Add(hl0,hl,norm,norm);
 AddHist(hl,hl0,hl,norm, norm) ;
 
 //TH1F* hno0 = FindRebin(f1, hname, EtaRange+"_Data", MinEta,MaxEta);
 hno0->Rebin(100/myNbins);

 //hno->Add(hno0,hl,norm,1);
 AddHist(hno,hno0,hl,norm,1) ;
 //TH1F* hc0 = FindRebin(g1, hname, EtaRange+"_Data", MinEta,MaxEta);
 hc0->Rebin(100/myNbins);
 //hc->Add(hc0,hno,norm,1);
 AddHist(hc,hc0,hno,norm,1) ;

 TH1F* H1= new TH1F("H1","", nbin, binmin, binmax);
 TH1F* H2= new TH1F("H2","", nbin, binmin, binmax);
 TH1F* H3= new TH1F("H3","", nbin, binmin, binmax);
       for (int i=1; i<nbin/2+1; i++) {
         H1->SetBinContent(i,h2->GetBinContent(i));  
         H2->SetBinContent(i,hc->GetBinContent(i));  
         H3->SetBinContent(i,hno->GetBinContent(i));
       }  

       //Double_t RangeUser_min = 0.5*10**( (Double_t) ( (int)(log(h1->GetMinimum())/log(10.)) ) );
       //Double_t RangeUser_max = 2.*10.**( (Double_t) ( (int)(log(h1->GetMaximum())/log(10.)) ) );

       Double_t RangeUser_min = 0.1;
       if( h1->GetMinimum() > 0 )
       {
         RangeUser_min = 10*( (Double_t) ( (int)(log(h1->GetMinimum())/log(10.)) ) );
         RangeUser_min = ((int) (h1->GetMinimum()/RangeUser_min) ) * RangeUser_min/2.;
       }
       //if( RangeUser_min <= 0 ) RangeUser_min = 1;

       RangeUser_min = 10 ;

       Double_t RangeUser_max = 10.*( (Double_t) ( (int)(log(h1->GetMaximum())/log(10.)) ) );
       RangeUser_max = ((int) (h1->GetMaximum()/RangeUser_max) ) * 20. * RangeUser_max;
       RangeUser_max = TMath::Power(10, int(log(h1->GetMaximum())/log(10))+4) ;

       cout << h1->GetMinimum() << "\t" << RangeUser_min << endl;
       cout << h1->GetMaximum() << "\t" << RangeUser_max << endl;
       //getchar();
       h1->GetYaxis()->SetRangeUser(RangeUser_min,RangeUser_max);
       h1->Draw("E"); 
       Double_t yMin[5] = { 0, 0, 0, 0, 0 };
       Double_t yMax[5] = { 0, 0, 0, 0, 0 };
       h1->SetMarkerStyle(20);h1->SetMarkerSize(0.7); 
       h1->SetMarkerColor(kBlack);
       h1->SetLineColor(kBlack);
       h1->SetLineStyle(1);
       h1->SetLineWidth(1);
//        hv->Draw("Hsame"); 
//        hv->SetFillColor(kOrange);
       h2->SetMarkerStyle(0); h2->SetMarkerSize(0.);
       h2->Draw("Hsame"); 
       h2->SetFillColor(2);
       H1->SetMarkerStyle(0); H1->SetMarkerSize(0.);
       H1->Draw("Hsame"); 
       H1->SetFillColor(2);
       hc->SetMarkerStyle(0); hc->SetMarkerSize(0.);
       hc->Draw("Hsame"); 
       hc->SetFillColor(3);
       H2->SetMarkerStyle(0); H2->SetMarkerSize(0.);
       H2->Draw("Hsame"); 
       H2->SetFillColor(3);
       hno->SetMarkerStyle(0); hno->SetMarkerSize(0.);
       hno->Draw("Hsame"); 
       hno->SetFillColor(4);
       H3->SetMarkerStyle(0); H3->SetMarkerSize(0.);
       H3->Draw("Hsame"); 
       H3->SetFillColor(kAzure-4);
//        hl->Draw("Hsame"); 
//        hl->SetFillColor(kBlue-7);
//        hv0->Draw("Hsame"); 
//        hv0->SetFillColor(kBlue-7);
//        h0->Draw("Hsame"); 
//        h0->SetFillColor(kAzure-4);
       h1->SetTickLength(0.03, "YZ");
       h1->SetTickLength(-0.03,"X");
       h1->SetLabelOffset(0.023,"X");
       h1->SetLabelOffset(0.015,"Y");
       h1->SetLabelSize(0.07, "XYZ");
       h1->SetLabelFont(42, "XYZ"); h1->SetTitleFont(42, "XYZ");
       h1->SetTitleSize(0.08, "XYZ"); h1->SetTitleOffset(0.9,"Y");
       h1->GetXaxis()->SetTitle(htitle); h1->GetYaxis()->SetTitle("Jets / "+tHBW);
       h1->SetNdivisions(509,"XYZ");
       h1->Draw("Esame"); 

for( int i = 0; i < 5; ++i )
{
  if( Cuts[i] > -5 )
  {
    yMin[i] = exp( log(h1->GetMinimum())+ ( log(h1->GetMaximum())-log(h1->GetMinimum()) )* minAlpha );
    //yMax[i] = 0.1*h1->GetBinContent( (Int_t) ( (Cuts[i]-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(3))+1 );
    //yMax[i] = yMin[i]*(h1->GetMaximum()/(100*h1->GetMinimum()));
    yMax[i] = exp( log(h1->GetMinimum())+ ( log(h1->GetMaximum())-log(h1->GetMinimum()) )* maxAlpha );
    cout << i << "\t" << yMin[i]  << "\t" << yMax[i] << endl;
  }
  //ta[i] = new TArrow(0.,1.,0.,2.,0.03,"<");
  ta[i] = new TArrow(Cuts[i],yMin[i],Cuts[i],yMax[i],0.025,"<");
  tl[i] = new TLatex(Cuts[i],yMax[i]*1.5,TaggerS[i]);
  tl[i]->SetTextAlign(21);
}
       for( int i = 0; i < 5; ++i )
       {
         ta[i]->Draw();
         tl[i]->Draw();
         
         //TaggerS


       }
//        h1->SetMinimum(1); 
//        h1->SetMaximum(3000000);  

  //TLegend* leg = new TLegend(0.20,0.40,0.45,0.95);
  TLegend* leg = NULL;
    if( tagger == "DeepFlavour" ) leg = new TLegend(0.484637,0.631282,0.695531,0.880553);
    else {
      //leg = new TLegend(0.63,0.63,0.88,0.92);
      leg = new TLegend(0.681745,0.640201,0.893376,0.890534);
    }
    leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetTextAlign(32);
    if( tagger == "DeepFlavour" ) leg->SetTextAlign(12);
/*
    TString headerText = "#splitline{P_{T} > ";
    if( cut == "400" ) headerText += "450 GeV";
    headerText += "}{" + HEtaRange + "}";
*/
    
    TString cut1 = "" ;
    if (cut == "0") cut1 = "0 GeV" ; 
    if (cut == "40") cut1 = "50 GeV" ; 
    if (cut == "60") cut1 = "70 GeV" ; 
    if (cut == "80") cut1 = "100 GeV" ; 
    if (cut == "140") cut1 = "160 GeV" ; 
    if (cut == "200") cut1 = "220 GeV" ; 
    if (cut == "260") cut1 = "310 GeV" ; 
    if (cut == "320") cut1 = "360 GeV" ; 
    if (cut == "400") cut1 = "450 GeV" ; 
    if (cut == "450") cut1 = "500 GeV" ; 
    if (cut == "500") cut1 = "550 GeV" ; 
    TString headerText = "P_{T} > " + cut1 ;
    headerText += ", " + HEtaRange + "";

    //if( cut == "40" ) leg->SetHeader("P_T > 50 GeV\n" + HEtaRange);
    //if( cut == "320" ) leg->SetHeader("P_T > 360 GeV\n" + HEtaRange);
    leg->SetHeader(headerText);
    leg->AddEntry(h1," Data","PL");
    leg->AddEntry(h2," b","F");
    leg->AddEntry(hc," c","F");
    leg->AddEntry(hno,"udsg","F");
    leg->Draw();
    
  //include the official CMS label
  CMS_lumi((TPad*)c1->cd(1),Lumi,10,true);

  rap1->cd();
  rap1->SetGridy(1);
  TH1F* g0= new TH1F("g0","",nbin,binmin,binmax);
   g0->Divide(h1,h2,1,1);
   g0->Draw("E"); 
   g0->SetLineColor(1);
   g0->SetMarkerStyle(20);
   g0->SetMarkerColor(1);
   g0->SetMarkerSize(0.5);
   g0->SetTickLength(0.03, "YZ");
   g0->SetTickLength(0.03*3,"X");
   g0->SetLabelOffset(0.02,"X");
   g0->SetLabelOffset(0.007,"Y");
   g0->SetLabelSize(0.04*3, "XYZ");
   g0->SetLabelFont(42, "XYZ"); g0->SetTitleFont(42, "XYZ");
   g0->SetTitleSize(0.05*3, "XYZ"); 
   g0->SetTitleOffset(0.45,"Y");
   g0->SetTitleOffset(0.9,"X");
   g0->GetXaxis()->SetTitle(htitle); g0->GetYaxis()->SetTitle("Data/MC");
   g0->SetNdivisions(509,"X"); g0->SetNdivisions(504,"Y");
   g0->SetMinimum(0); g0->SetMaximum(2);  
   TF1 *Fun = new TF1("Fun","1.4",binmin, binmax);
    Fun->SetLineWidth(0.5);
    Fun->SetLineStyle(2);
    Fun->Draw("same");
    Fun = new TF1("Fun","1.2",binmin, binmax);
    Fun->SetLineWidth(0.5);
    Fun->SetLineStyle(2);
    Fun->Draw("same");
    Fun = new TF1("Fun","1",binmin, binmax);
    Fun->SetLineWidth(0.5);
    Fun->SetLineStyle(2);
    Fun->Draw("same");
    Fun = new TF1("Fun","0.8",binmin, binmax);
    Fun->SetLineWidth(0.5);
    Fun->SetLineStyle(2);
    Fun->Draw("same");
    Fun = new TF1("Fun","0.6",binmin, binmax);
    Fun->SetLineWidth(0.5);
    Fun->SetLineStyle(2);
    Fun->Draw("same");

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //return 1 ;
  c1->Update();
  if (ImprovePlot) c1->WaitPrimitive();
  std::cout << "\n Save tagger plot now" << std::endl;
  c1->SaveAs(SavDir + "/tagger_"+tagger+"_"+EtaRange.ReplaceAll(".","p")+"_PFCut_"+cut+".png");
  c1->SaveAs(SavDir + "/tagger_"+tagger+"_"+EtaRange.ReplaceAll(".","p")+"_PFCut_"+cut+".gif");
  c1->SaveAs(SavDir + "/tagger_"+tagger+"_"+EtaRange.ReplaceAll(".","p")+"_PFCut_"+cut+".C");
  c1->SaveAs(SavDir + "/tagger_"+tagger+"_"+EtaRange.ReplaceAll(".","p")+"_PFCut_"+cut+".pdf");
  c1->SaveAs(SavDir + "/tagger_"+tagger+"_"+EtaRange.ReplaceAll(".","p")+"_PFCut_"+cut+".root");

  return 1 ;
}
//plot(TString tagger="JPL", TString cut="300", int myNbins = 50)

