// Author: Pierre Van Hove 04/03/2014

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TH1Comp                                                              //
//                                                                      //
// A class for easy TH1F compareason.					//
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TH1Comp.h"

#include <iostream>
#include "TROOT.h"
#include "TObject.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

ClassImp(TH1Comp)

//-- Default Constructor -------------------------------------------------
TH1Comp::TH1Comp()
{
  _LegendTitle = "";
  _TH1Comp_debug = kFALSE ;
  _NewIsSubsetOfRef = kFALSE;
  _RefHis = NULL;
  _NewHis = NULL;
  _RebinnedRefHis = NULL;
  _RebinnedNewHis = NULL;
  _RatHis = NULL;
  _DifHis = NULL;
  _RelDifHis = NULL;
  _NRebinnedRefHis = NULL;
  _NRebinnedNewHis = NULL;
  _NRatHis = NULL;
  _RefCol = 1 ;
  _NewCol = 2 ;
  _opt = "same";
  _RefTitle = "Ref";
  _NewTitle = "New";
  _RatTitle = _NewTitle + " / " + _RefTitle ;
}
//-- Another constructor -------------------------------------------------
TH1Comp::TH1Comp(TH1F* hRef, TH1F* hNew, Color_t cRef, Color_t cNew, Bool_t IsSubset, TString LegendTitle, TString RefTitle, TString NewTitle)
{
  _LegendTitle = LegendTitle;
  _RefTitle = RefTitle;
  _NewTitle = NewTitle;
  _RatTitle = _NewTitle + " / " + _RefTitle ;
  if( hRef == NULL || hNew == NULL )
  {
    std::cout << "At least one of the histogram pointer is not existing " << std::endl;
    return;
  }
  _TH1Comp_debug = kFALSE ;
  _NewIsSubsetOfRef = IsSubset;
  _RebinnedRefHis = NULL;
  _RebinnedNewHis = NULL;
  _NRebinnedRefHis = NULL;
  _NRebinnedNewHis = NULL;
  _RefHis = (TH1F*) hRef->Clone();
  _NewHis = (TH1F*) hNew->Clone();
  _RefCol = cRef ;
  _NewCol = cNew ;
  _opt = "same";

  this->Update();
}
//-- Another constructor from files with names... ------------------------
//
TH1Comp::TH1Comp(TString fName1, TString fName2, TString hName, Color_t cRef, Color_t cNew, Bool_t IsSubset, TString LegendTitle, TString RefTitle, TString NewTitle)
{
  _LegendTitle = LegendTitle;
  _RefTitle = RefTitle;
  _NewTitle = NewTitle;
  _RatTitle = _NewTitle + " / " + _RefTitle ;
  TFile *f1 = new TFile(fName1);
  if( f1 == NULL )
  {
    std::cout << "Failed to open file " << fName1 << "." << std::endl;
    return;
  }
  TString chn = "Ref_"+hName;
  std::cout << "chn = " << chn << std::endl;
  TString h300Name = "";
  if( hName == "cMCRlightWeight_" || hName == "MCRlightWeight_" )
  {
    h300Name = hName + "Jet300";
  }
  else
  {
    h300Name = hName;
  }
  TH1F *hRef = (TH1F*) gROOT->FindObject(h300Name)->Clone(chn);
  std::cout << "Suceded to get chn = " << chn << std::endl;
  if( hRef == NULL )
  {
    std::cout << "Failed to get histo " << h300Name << " From file " << fName1 << "." << std::endl;
    return;
  }
  hRef->SetDirectory(0);
  f1->Close();

  TFile *f2 = new TFile(fName2);
  if( f2 == NULL )
  {
    std::cout << "Failed to open file " << fName2 << "." << std::endl;
    return;
  }
  chn = "New_"+hName;
  if( hName == "cMCRlightWeight_" || hName == "MCRlightWeight_" )
  {
    h300Name = hName + "300";
  }
  else
  {
    h300Name = hName;
  }
  TH1F *hNew = (TH1F*) gROOT->FindObject(h300Name)->Clone(chn);
  if( hNew == NULL )
  {
    std::cout << "Failed to get histo " << h300Name << " From file " << fName2 << "." << std::endl;
    return;
  }
  hNew->SetDirectory(0);
  f2->Close();

  _TH1Comp_debug = kFALSE ;
  _NewIsSubsetOfRef = IsSubset;
  _RebinnedRefHis = NULL;
  _RebinnedNewHis = NULL;
  _NRebinnedRefHis = NULL;
  _NRebinnedNewHis = NULL;
  _RefHis = (TH1F*) hRef->Clone();
  _NewHis = (TH1F*) hNew->Clone();
  _RefCol = cRef ;
  _NewCol = cNew ;
  _opt = "same";
  this->Update();

}
//
//-- Copy constructor ----------------------------------------------------
TH1Comp::TH1Comp(TH1Comp const& CompToCopy):TObject(CompToCopy)
{
  _RefTitle = CompToCopy._RefTitle;
  _NewTitle = CompToCopy._NewTitle;
  _RatTitle = CompToCopy._RatTitle;

  _TH1Comp_debug = CompToCopy._TH1Comp_debug ;
  _RefCol = CompToCopy._RefCol ;
  _NewCol = CompToCopy._NewCol ;
  _RebinnedRefHis = CompToCopy._RebinnedRefHis;
  _RebinnedNewHis = CompToCopy._RebinnedNewHis;
  _NRebinnedRefHis = CompToCopy._NRebinnedRefHis;
  _NRebinnedNewHis = CompToCopy._NRebinnedNewHis;
  _opt = CompToCopy._opt ;

  if( CompToCopy._RefHis == NULL )
   _RefHis = NULL;
  else
   _RefHis = (TH1F*) CompToCopy._RefHis->Clone();
  if( CompToCopy._NewHis == NULL )
   _NewHis = NULL;
  else
   _NewHis = (TH1F*) CompToCopy._NewHis->Clone();

  if( CompToCopy._RatHis == NULL )
   _RatHis = NULL;
  else
   _RatHis = (TH1F*) CompToCopy._RatHis->Clone();
  if( CompToCopy._NRatHis == NULL )
   _NRatHis = NULL;
  else
   _NRatHis = (TH1F*) CompToCopy._NRatHis->Clone();

  if( CompToCopy._DifHis == NULL )
   _DifHis = NULL;
  else
   _DifHis = (TH1F*) CompToCopy._DifHis->Clone();

  if( CompToCopy._RelDifHis == NULL )
   _RelDifHis = NULL;
  else
   _RelDifHis = (TH1F*) CompToCopy._RelDifHis->Clone();
}
//-- Destructor ----------------------------------------------------------
TH1Comp::~TH1Comp()
{
  if( _RefHis != NULL ) _RefHis->Delete();
  if( _NewHis != NULL ) _NewHis->Delete();
  if( _RebinnedRefHis != NULL ) _RebinnedRefHis->Delete();
  if( _RebinnedNewHis != NULL ) _RebinnedNewHis->Delete();
  if( _RatHis != NULL ) _RatHis->Delete();
  if( _DifHis != NULL ) _DifHis->Delete();
  if( _RelDifHis != NULL ) _RelDifHis->Delete();
  if( _NRebinnedRefHis != NULL ) _NRebinnedRefHis->Delete();
  if( _NRebinnedNewHis != NULL ) _NRebinnedNewHis->Delete();
  if( _NRatHis != NULL ) _NRatHis->Delete();
}
//------------------------------------------------------------------------
void TH1Comp::SetRef(TH1F* hRef)
{
  if(_TH1Comp_debug)
  this->Info("TH1Comp::SetRef(TH1F* hRef)",
             "processing");
  _RefHis = (TH1F*) hRef->Clone();
}
//------------------------------------------------------------------------
void TH1Comp::SetNew(TH1F* hNew)
{
  if(_TH1Comp_debug)
  this->Info("TH1Comp::SetNew(TH1F* hNew)",
             "processing");
  _NewHis = (TH1F*) hNew->Clone();
}
//------------------------------------------------------------------------
void TH1Comp::SetIsSubset(Bool_t IsSub)
{
  _NewIsSubsetOfRef = IsSub;
  this->Update();
}
//------------------------------------------------------------------------
void TH1Comp::SetHisStyle(TH1F* h)
{
  h->GetXaxis()->SetTitleFont(_RefHis->GetXaxis()->GetTitleFont());
  h->GetXaxis()->SetTitleSize(_RefHis->GetXaxis()->GetTitleSize());
  h->GetXaxis()->SetLabelFont(_RefHis->GetXaxis()->GetLabelFont());
  h->GetXaxis()->SetLabelSize(_RefHis->GetXaxis()->GetLabelSize());

  h->GetYaxis()->SetTitleFont(_RefHis->GetYaxis()->GetTitleFont());
  h->GetYaxis()->SetTitleSize(_RefHis->GetYaxis()->GetTitleSize());
  h->GetYaxis()->SetTitleOffset(_RefHis->GetYaxis()->GetTitleOffset());
  h->GetYaxis()->SetLabelFont(_RefHis->GetYaxis()->GetLabelFont());
  h->GetYaxis()->SetLabelSize(_RefHis->GetYaxis()->GetLabelSize());
/*
  h->GetXaxis()->SetTitleFont(62);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelSize(0.04);

  h->GetYaxis()->SetTitleFont(62);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetLabelFont(62);
  h->GetYaxis()->SetLabelSize(0.04);
*/
}
//------------------------------------------------------------------------
void TH1Comp::Update()
{
  Double_t *dRefBins = NULL;
  Double_t *dNewBins = NULL;
  Double_t *dOutBins = NULL;


  int nRefBins = 0;
  int nNewBins = 0;
  int nOutBins = 0;
  Bool_t bRef = kFALSE;
  Bool_t bNew = kFALSE;

  if(this->_TH1Comp_debug)
  this->Info("TH1Comp::Update()",
             "processing");
  if( this->_RefHis == NULL || this->_NewHis == NULL )
  {
    std::cout << "At least one of the histogram pointer is not existing " << std::endl;
    return;
  }
  this->_RefHis->SetLineColor(this->_RefCol);
  this->_RefHis->SetMarkerColor(this->_RefCol);
  this->_RefHis->SetMarkerStyle(21);
  this->_NewHis->SetLineColor(this->_NewCol);
  this->_NewHis->SetMarkerColor(this->_NewCol);
  this->_NewHis->SetMarkerStyle(20);

  this->_RefHis->Sumw2();
  this->_NewHis->Sumw2();
  // Also compute the ratio, diff and relative ratio histograms
  // Here a checks need to be done for compatibility
  nRefBins = this->_RefHis->GetNbinsX()+1;
  dRefBins = new Double_t [nRefBins];
  if( this->_RefHis->GetXaxis()->GetXbins()->GetSize() == 0 ) // Histos has constant binning
  {
    bRef = kTRUE;
    for(int i = 0; i < nRefBins; ++i)
    {
      dRefBins[i] = this->_RefHis->GetXaxis()->GetXmin() + i * this->_RefHis->GetXaxis()->GetBinWidth(i);
    }
  } else {						// Histo has variable binning
    const Double_t *arRef = this->_RefHis->GetXaxis()->GetXbins()->GetArray();
    for(int i = 0; i < nRefBins; ++i)
    {
      dRefBins[i] = arRef[i];
    }
  }
  nNewBins = this->_NewHis->GetNbinsX()+1;
  dNewBins = new Double_t [nNewBins];
  if( this->_NewHis->GetXaxis()->GetXbins()->GetSize() == 0 ) // Histos has constant binning
  {
    bNew = kTRUE;
    for(int i = 0; i < nNewBins; ++i)
    {
      dNewBins[i] = this->_NewHis->GetXaxis()->GetXmin() + i * this->_NewHis->GetXaxis()->GetBinWidth(i);
    }
  } else {						// Histo has variable binning
    const Double_t *arNew = this->_NewHis->GetXaxis()->GetXbins()->GetArray();
    for(int i = 0; i < nNewBins; ++i)
    {
      dNewBins[i] = arNew[i];
    }
  }

  UInt_t nMaxBins = nNewBins;
  if( nRefBins > nNewBins ) nMaxBins = nRefBins;
  dOutBins = new Double_t [nMaxBins];

  nOutBins = 0;
  int jmin = 0;
  Double_t precRef = (this->_RefHis->GetXaxis()->GetXmax()-this->_RefHis->GetXaxis()->GetXmin())/this->_RefHis->GetNbinsX();
  Double_t precNew = (this->_NewHis->GetXaxis()->GetXmax()-this->_NewHis->GetXaxis()->GetXmin())/this->_NewHis->GetNbinsX();
  Double_t BestPrec = min(precRef, precNew)/1000.;
  for(int i = 0; i < nRefBins; ++i)
  {
    for(int j = jmin; j < nNewBins; ++j)
    {
      if( (dRefBins[i]-dNewBins[j]) < BestPrec && (dRefBins[i]-dNewBins[j]) > -BestPrec )
      {
        dOutBins[nOutBins] = dRefBins[i];
        ++nOutBins;
        jmin = j;
        break;
      }
    }
  }
  if( nOutBins < 2 )
  {
    cout << "The two histograms are not compatible at all" << endl;
    cout << "less than two bin edges correspons" << endl;
    return;
  } else {
    cout << "Found " << nOutBins << " compatible bin edges so create hist with this bin number -1" << endl;
/*
    this->_RebinnedRefHis = new TH1F("_RebinnedRefHis", "_RebinnedRefHis", nOutBins-1, dOutBins);
    this->_RebinnedNewHis = new TH1F("_RebinnedNewHis", "_RebinnedNewHis", nOutBins-1, dOutBins);
    this->_NRebinnedRefHis = new TH1F("_NRebinnedRefHis", "_NRebinnedRefHis", nOutBins-1, dOutBins);
    this->_NRebinnedNewHis = new TH1F("_NRebinnedNewHis", "_NRebinnedNewHis", nOutBins-1, dOutBins);
*/

    this->_RebinnedRefHis = new TH1F("_RebinnedRefHis", this->_RefHis->GetTitle(), nOutBins-1, dOutBins);
    this->_RebinnedNewHis = new TH1F("_RebinnedNewHis", this->_RefHis->GetTitle(), nOutBins-1, dOutBins);
    this->_NRebinnedRefHis = new TH1F("_NRebinnedRefHis", this->_RefHis->GetTitle(), nOutBins-1, dOutBins);
    this->_NRebinnedNewHis = new TH1F("_NRebinnedNewHis", this->_RefHis->GetTitle(), nOutBins-1, dOutBins);

    this->_RebinnedRefHis->GetXaxis()->SetTitle(this->_RefHis->GetXaxis()->GetTitle());
    this->_RebinnedRefHis->SetLineColor(this->_RefCol);
    this->_RebinnedRefHis->SetMarkerColor(this->_RefCol);
    this->_RebinnedNewHis->GetXaxis()->SetTitle(this->_NewHis->GetXaxis()->GetTitle());
    this->_RebinnedNewHis->SetLineColor(this->_NewCol);
    this->_RebinnedNewHis->SetMarkerColor(this->_NewCol);

    this->_NRebinnedRefHis->GetXaxis()->SetTitle(this->_RefHis->GetXaxis()->GetTitle());
    this->_NRebinnedRefHis->SetLineColor(this->_RefCol);
    this->_NRebinnedRefHis->SetMarkerColor(this->_RefCol);
    this->_NRebinnedNewHis->GetXaxis()->SetTitle(this->_NewHis->GetXaxis()->GetTitle());
    this->_NRebinnedNewHis->SetLineColor(this->_NewCol);
    this->_NRebinnedNewHis->SetMarkerColor(this->_NewCol);

    Double_t iCont = 0;
    for( int i = 0; i < nOutBins+1; ++i)
    {
      Double_t cont = 0;
      Double_t err2 = 0;
      for( int j = 0; j < nRefBins+1; ++j)
      {
        if( (i == 0 && this->_RefHis->GetBinLowEdge(j+1) < dOutBins[0])
         || (i == nOutBins && this->_RefHis->GetBinLowEdge(j) >= dOutBins[nOutBins-1])
         || ( i > 0 && i < nOutBins && this->_RefHis->GetBinLowEdge(j) >= dOutBins[i-1]
                    && this->_RefHis->GetBinLowEdge(j) < dOutBins[i] ) )
        {
          cont += this->_RefHis->GetBinContent(j);
          err2 += (this->_RefHis->GetBinError(j)*this->_RefHis->GetBinError(j));
        }
      }
      iCont += cont;
      this->_RebinnedRefHis->SetBinContent(i, cont);
      this->_RebinnedRefHis->SetBinError(i, sqrt(err2));
      this->_NRebinnedRefHis->SetBinContent(i, cont);
      this->_NRebinnedRefHis->SetBinError(i, sqrt(err2));
      cont = err2 = 0;
    }
    _NRebinnedRefHis->Scale(1/iCont);
    
    iCont = 0;
    for( int i = 0; i < nOutBins+1; ++i)
    {
      Double_t cont = 0;
      Double_t err2 = 0;
      for( int j = 0; j < nNewBins+1; ++j)
      {
        if( (i == 0 && this->_NewHis->GetBinLowEdge(j+1) < dOutBins[0])
         || (i == nOutBins && this->_NewHis->GetBinLowEdge(j) >= dOutBins[nOutBins-1])
         || ( i > 0 && i < nOutBins && this->_NewHis->GetBinLowEdge(j) >= dOutBins[i-1]
                    && this->_NewHis->GetBinLowEdge(j) < dOutBins[i] ) )
        {
          cont += this->_NewHis->GetBinContent(j);
          err2 += (this->_NewHis->GetBinError(j)*this->_NewHis->GetBinError(j));
        }
      }
      iCont += cont;
      this->_RebinnedNewHis->SetBinContent(i, cont);
      this->_RebinnedNewHis->SetBinError(i, sqrt(err2));
      this->_NRebinnedNewHis->SetBinContent(i, cont);
      this->_NRebinnedNewHis->SetBinError(i, sqrt(err2));
      cont = err2 = 0;
    }
    _NRebinnedNewHis->Scale(1/iCont);

    this->_RatHis = (TH1F*) this->_RebinnedNewHis->Clone("RatHis");
    this->_NRatHis = (TH1F*) this->_NRebinnedNewHis->Clone("NRatHis");
    this->_DifHis = (TH1F*) this->_RebinnedNewHis->Clone("DifHis");
    this->_RelDifHis = (TH1F*) this->_RebinnedNewHis->Clone("RelDifHis");
  }
  this->_RatHis->Sumw2();
  this->_NRatHis->Sumw2();
  this->_DifHis->Sumw2();
  this->_RelDifHis->Sumw2();
  this->_RebinnedNewHis->Sumw2();
  this->_RebinnedRefHis->Sumw2();

  if( this->_NewIsSubsetOfRef )this->_RatHis->Divide(this->_RebinnedNewHis, this->_RebinnedRefHis, 1, 1, "B");
  else this->_RatHis->Divide(this->_RebinnedNewHis, this->_RebinnedRefHis, 1, 1);
  if( this->_NewIsSubsetOfRef )this->_NRatHis->Divide(this->_NRebinnedNewHis, this->_NRebinnedRefHis, 1, 1, "B");
  else this->_NRatHis->Divide(this->_NRebinnedNewHis, this->_NRebinnedRefHis, 1, 1);
  this->_DifHis->Add(this->_RebinnedNewHis, this->_RebinnedRefHis, 1, -1);
  this->_RelDifHis->Divide(this->_DifHis, this->_RebinnedRefHis);

  this->_RefHis->SetStats(0);
  this->_NewHis->SetStats(0);

  this->_RatHis->SetStats(0);
  this->_NRatHis->SetStats(0);
  this->_DifHis->SetStats(0);
  this->_RelDifHis->SetStats(0);
  this->_RebinnedNewHis->SetStats(0);
  this->_RebinnedRefHis->SetStats(0);

/*
  this->_RefHis->GetXaxis()->SetTitleFont(62);
  this->_RefHis->GetXaxis()->SetTitleSize(0.04);
  this->_RefHis->GetXaxis()->SetLabelFont(62);
  this->_RefHis->GetXaxis()->SetLabelSize(0.04);

  this->_RefHis->GetYaxis()->SetTitleFont(62);
  this->_RefHis->GetYaxis()->SetTitleSize(0.04);
  this->_RefHis->GetYaxis()->SetTitleOffset(1.2);
  this->_RefHis->GetYaxis()->SetLabelFont(62);
  this->_RefHis->GetYaxis()->SetLabelSize(0.04);
*/


  SetHisStyle(_RefHis);
  SetHisStyle(_NewHis);
  SetHisStyle(_RebinnedRefHis);
  SetHisStyle(_RebinnedNewHis);
  SetHisStyle(_RatHis);
  SetHisStyle(_DifHis);
  SetHisStyle(_RelDifHis);
  SetHisStyle(_NRatHis);		// Ratio between normalized histos
  SetHisStyle(_NRebinnedRefHis);	// Normalised his
  SetHisStyle(_NRebinnedNewHis);	// Normalised his
}
//------------------------------------------------------------------------
void TH1Comp::DrawBoth(TString option, Double_t xmin, Double_t xmax, TCanvas *tc , Color_t RatCol)
{
  //Color_t RatCol = 51;
  
  TCanvas *cp = tc;
  TLegend *tl1 = NULL;
  TLegend *tl2 = NULL;
  if(_TH1Comp_debug)
  this->Info("TH1Comp::DrawBoth()",
             "processing");

  if( option != "" ) _opt = option;

  if( cp == NULL ) // Needs to define a Canvas !
  {
    cp = (TCanvas*)gROOT->GetSelectedPad();
    if( cp == NULL )
    {
      //TCanvas *tc = new TCanvas("cp","cp",600,750);
      //tc->cd(1);
      //cp = (TPad*)gROOT->GetSelectedPad();
      cout << "cp is NULL" << endl;
      getchar();
      cp = new TCanvas("cp","cp",600,750);
    }
  }
  if( _RefHis == NULL || _NewHis == NULL )
  {
    std::cout << "At least one of the histogram pointer is not existing " << std::endl;
    return;
  }
  this->Update();

  cout << "_RefHis Title is " << _RefHis->GetTitle() << endl;
  cout << "_NewHis Title is " << _NewHis->GetTitle() << endl;
  if ( _opt.Contains("same")  && ! (_opt.Contains("rat") ))
  {
    _RefHis->Draw("same");
    _NewHis->Draw("same");
    tl1 = new TLegend(0.56,0.18,0.86,0.39);
    tl1->SetBorderSize(0);
    tl1->SetHeader(_LegendTitle);
    tl1->SetFillColor(0);
    tl1->AddEntry(_RefHis,_RefTitle,"PL");
    tl1->AddEntry(_NewHis,_NewTitle,"PL");
    tl1->Draw();
  } else if( _opt.Contains("both") ) {
    _RefHis->Draw();
    _NewHis->Draw("same");
  } else if( _opt.Contains("rat") ) {
    if( !(_opt.Contains("same")) )cp->Divide(1,2);
    cp->cd(1);
    tl1 = new TLegend(0.56,0.18,0.86,0.39);
    tl1->SetBorderSize(0);
    tl1->SetHeader(_LegendTitle);
    tl1->SetFillColor(0);
    tl1->AddEntry(_RefHis,_RefTitle,"PL");
    tl1->AddEntry(_NewHis,_NewTitle,"PL");
    //cp->cd(1)->SetLogy(1);
  _RefHis->SetMinimum(0);
  _RefHis->GetYaxis()->SetRangeUser(xmin,xmax);
  _RefHis->SetTitle("");
    if( (_opt.Contains("same")) ) _RefHis->Draw("same");
    else _RefHis->Draw();
    _NewHis->Draw("same");
  cout << "1_RefHis Title is " << _RefHis->GetTitle() << endl;
  cout << "1_NewHis Title is " << _NewHis->GetTitle() << endl;
    tl1->Draw();
    cp->cd(2);
    cp->cd(2)->SetGridy(1);
    //tl2 = new TLegend(0.25,0.32,0.68,0.48);
    //tl2->SetBorderSize(0);
    //tl2->SetFillColor(0);
    //tl2->SetFillColor(0);
    //tl2->AddEntry(_RatHis,_RatTitle,"PL");
  this->_RatHis->SetMarkerStyle(20);
    _RatHis->SetMarkerColor(RatCol);
    _RatHis->SetLineColor(RatCol);
    _RatHis->GetYaxis()->SetTitle(_RatTitle);
    //_RatHis->SetMinimum(0.0);
    //_RatHis->SetMaximum(2.0);
    _RatHis->GetYaxis()->SetRangeUser(xmin,xmax);
    _RatHis->SetStats(1);
    _RatHis->SetTitle("");
    //_RatHis->Fit("pol0","0");
    _RatHis->Fit("pol0","");
    if( (_opt.Contains("same")) ) _RatHis->Draw("same");
    else _RatHis->Draw();
    //TF1 *tf = (TF1*) _RatHis->GetFunction("pol0");
    //tf->SetLineWidth(2)
    //_RatHis->Draw("same");
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TPaveStats *st = (TPaveStats*) _RatHis->FindObject("stats");
    st->SetX1NDC(0.736043);
    st->SetY1NDC(0.741893);
    st->SetX2NDC(0.967687);
    st->SetY2NDC(0.882553);
    st->Draw();
    
    //tl2->Draw();
  } else if( _opt.Contains("nRat") ) {
    gStyle->SetOptStat(0);
    cp->Divide(1,2);
    cp->cd(1);
    //cp->cd(1)->SetLogy(1);
    //tl1 = new TLegend(0.56,0.18,0.86,0.39);
    //tl1->SetBorderSize(0);
    //tl1->SetHeader(_LegendTitle);
    //tl1->SetFillColor(0);
    //tl1->AddEntry(_RefHis,_RefTitle,"PL");
    //tl1->AddEntry(_NewHis,_NewTitle,"PL");

    tl1 = new TLegend(0.56,0.18,0.86,0.39);
    //tl1 = new TLegend(0.66,0.38,0.96,0.78);
    tl1->SetBorderSize(0);
    tl1->SetHeader(_LegendTitle);
    tl1->SetFillColor(0);
    tl1->AddEntry(_RefHis,_RefTitle,"PL");
    tl1->AddEntry(_NewHis,_NewTitle,"PL");
    _NRebinnedRefHis->GetXaxis()->SetRangeUser(xmin,xmax);
    _NRebinnedRefHis->Draw();
    _NRebinnedNewHis->Draw("same");
  cout << "1_RefHis Title is " << _NRebinnedRefHis->GetTitle() << endl;
  cout << "1_NewHis Title is " << _NRebinnedNewHis->GetTitle() << endl;
    tl1->Draw();
    cp->cd(2);
    cp->cd(2)->SetGridy(1);
    //tl2 = new TLegend(0.66,0.38,0.96,0.78);
    //tl2->SetFillColor(0);
    //tl2->AddEntry(_RatHis,_RatTitle,"PL");
    _NRatHis->SetMinimum(2.0);
    _NRatHis->SetMaximum(2.0);
    _NRatHis->SetStats(1);
    _NRatHis->SetTitle(_RatHis->GetTitle());
    _NRatHis->GetYaxis()->SetTitle(_RatTitle);
    _NRatHis->SetMarkerColor(RatCol);
    _NRatHis->GetXaxis()->SetRangeUser(xmin,xmax);
    _NRatHis->SetLineColor(RatCol);
    _NRatHis->Draw();
    _NRatHis->Fit("pol0");
/*
fX1NDC                        0.736043            X1 point in NDC coordinates
fY1NDC                        0.741893            Y1 point in NDC coordinates
fX2NDC                        0.967687            X2 point in NDC coordinates
fY2NDC                        0.882553   
*/
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TPaveStats *st = (TPaveStats*) _NRatHis->FindObject("stats");
    st->SetX1NDC(0.736043);
    st->SetY1NDC(0.741893);
    st->SetX2NDC(0.967687);
    st->SetY2NDC(0.882553);
    st->Draw();
    //tl2->Draw();
  } else if( _opt.Contains("dif") ) {
    cp->Divide(1,2);
    cp->cd(1);
    _RefHis->Draw();
    _NewHis->Draw("same");
    cp->cd(2);
    _DifHis->Draw();
  } else if( _opt.Contains("rel") ) {
    cp->Divide(1,2);
    cp->cd(1);
    _RefHis->Draw();
    _NewHis->Draw("same");
    cp->cd(2);
    _RelDifHis->Draw();
  }
}
//------------------------------------------------------------------------
void TH1Comp::DrawRat(Option_t *option)
{
 if( _RatHis == NULL )
 {
   cout << "Cannot Draw because _RatHis == NULL " << endl;
   return;
 }
 _RatHis->Draw(option);
}
//------------------------------------------------------------------------
void TH1Comp::DrawDif(Option_t *option)
{
 if( _DifHis == NULL )
 {
   cout << "Cannot Draw because _DifHis == NULL " << endl;
   return;
 }
 _DifHis->Draw(option);
}
//------------------------------------------------------------------------
void TH1Comp::DrawRelDif(Option_t *option)
{
 if( _RelDifHis == NULL )
 {
   cout << "Cannot Draw because _RelDifHis == NULL " << endl;
   return;
 }
 _RelDifHis->Draw(option);
}
//------------------------------------------------------------------------
