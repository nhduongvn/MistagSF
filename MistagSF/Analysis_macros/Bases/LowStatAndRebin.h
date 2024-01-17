#ifndef LowStatAndRebin_H
#define LowStatAndRebin_H

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

//Double_t *xbins = NULL;
//int xbinsSize = 0;

Double_t *xbins = NULL;
Double_t nbins = 0;

void setBinVector(int size, Double_t *vec);
void setLowStatBinToZero(TH1F *h, std::vector<Int_t> &LowStatBins);
TH1F *FindRebin(TFile *f, TString hNAme, TString NameSuffix, unsigned int MinEta, unsigned int MaxEta);
TH1F *FindRebinSetLowStatToZero(TFile *f, TString hNAme, TString NameSuffix, unsigned int MinEta, unsigned int MaxEta, std::vector<Int_t> &LowStatBins );
void GetLowStatBins(TFile *f, TString hName, unsigned int MinEta, unsigned int MaxEta, Double_t xSigma, std::vector<Int_t> &LowStatBins);

#endif
