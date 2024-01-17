// Author: Pierre Van Hove 09/12/2015

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// LowStatAndRebin                                                      //
//                                                                      //
// Methods to get hist from file, rebin them and remove low stat bins	//
//                                                                      //
//////////////////////////////////////////////////////////////////////////
/*
#include <fstream>
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <TROOT.h>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TSystem.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
*/

#include "LowStatAndRebin.h"

//Double_t *xbins = NULL;
//Double_t nbins = 0;

void setBinVector(int size, Double_t *vec)
{
  nbins = size;
  xbins = new Double_t[size];
  for(int i = 0; i < size; ++i )
  {
    xbins[i] = vec[i];
  }
}

void setLowStatBinToZero(TH1F *h, std::vector<Int_t> &LowStatBins)
{
  if( LowStatBins.size() > 0 )
  {
    for( int iLS = 0; iLS < (int) LowStatBins.size(); ++iLS)
    {
      h->SetBinContent(LowStatBins.at(iLS),0);
      h->SetBinError(LowStatBins.at(iLS),0);
    }
  } else {
  }
  return;
}

TH1F *FindRebin(TFile *f, TString hNAme, TString NameSuffix, unsigned int MinEta, unsigned int MaxEta)
{
  f->cd();
 
  TH2F *temp = NULL;
  TH1F* ret = NULL;

  if( hNAme != "" )
  {
    //cout << "\n====================================================================" ;
    //cout << "\n Name is : " << f->GetName() << "  " << hNAme ;
    //cout << "\n====================================================================" ;
    temp = (TH2F*) gROOT->FindObject(hNAme);
    //TCanvas* cTest = new TCanvas("cTest") ;
    //temp->Draw() ;
    ret = (TH1F*) ( temp->ProjectionX(hNAme+NameSuffix, MinEta,MaxEta) );
    if( nbins > 0 ) ret = (TH1F*) ret->Rebin(nbins-1,ret->GetName(),xbins);
  }

  return ret;
}

TH1F *FindRebinSetLowStatToZero(TFile *f, TString hNAme, TString NameSuffix, unsigned int MinEta, unsigned int MaxEta, std::vector<Int_t> &LowStatBins )
{
  f->cd();
  TH1F *ret = (TH1F*) FindRebin(f,hNAme, NameSuffix, MinEta, MaxEta );
  setLowStatBinToZero(ret, LowStatBins);

  return ret;
}

void GetLowStatBins(TFile *f, TString hName, unsigned int MinEta, unsigned int MaxEta, Double_t xSigma, std::vector<Int_t> &LowStatBins)
{
  f->cd();
  TH1F *h = FindRebin(f,hName,"_temp",MinEta,MaxEta);

  LowStatBins.resize(0);
  Double_t BinError, BinContent;
  Int_t Nbins = h->GetNbinsX();
  for( int ibin = 0; ibin <= Nbins+1; ++ibin )
  {
    BinContent = h->GetBinContent(ibin);
    BinError = h->GetBinError(ibin);
    if( BinContent <= xSigma * BinError )
    {
      LowStatBins.push_back(ibin);
    }
  }
  return;
}
