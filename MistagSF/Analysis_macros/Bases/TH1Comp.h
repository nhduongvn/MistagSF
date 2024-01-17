#ifndef TH1Comp_H
#define TH1Comp_H

/** @class TH1Comp
* @ A class to compare 1D histograms
*  
* @author Pierre Van Hove
* $Header$
*/

#include "TH1F.h"
#include "TCanvas.h"

class TH1Comp:public TObject{
  private :
  TH1F *_RefHis;
  TH1F *_NewHis;
  TH1F *_RebinnedRefHis;
  TH1F *_RebinnedNewHis;
  TH1F *_RatHis;
  TH1F *_DifHis;
  TH1F *_RelDifHis;
  TH1F *_NRatHis;		// Ratio between normalized histos
  TH1F *_NRebinnedRefHis;	// Normalised his
  TH1F *_NRebinnedNewHis;	// Normalised his

  TString _LegendTitle;
  TString _opt;
  Color_t _RefCol;
  Color_t _NewCol;
  Color_t _RatCol;
  Color_t _DifCol;
  Color_t _RelDifCol;
  Bool_t _TH1Comp_debug; 
  Bool_t _NewIsSubsetOfRef;
  TString _RatTitle;
  TString _RefTitle;
  TString _NewTitle;

  protected:

  public :
  // Accessible methods //
  TH1Comp();
  TH1Comp(TH1F* hRef, TH1F* hNew, Color_t cRef = 1, Color_t cNew = 2, Bool_t IsSubset = kFALSE, TString LegendTitle = "", TString RefTitle = "Ref", TString NewTitle="New");
  TH1Comp(TString fName1, TString fName2, TString hName, Color_t cRef = 1, Color_t cNew = 2, Bool_t IsSubset = kFALSE, TString LegendTitle = "", TString RefTitle = "Ref", TString NewTitle="New");
  TH1Comp(TH1Comp const& CompToCopy);
  virtual ~TH1Comp();

  TH1F* GetRefHist() {return _RefHis;}
  TH1F* GetNewHist() {return _NewHis;}
  TH1F* GetRebinnedRefHist() {return _RebinnedRefHis;}
  TH1F* GetRebinnedNewHist() {return _RebinnedNewHis;}
  TH1F* GetRatHist() {return _RatHis;}
  TH1F* GetDifHist() {return _DifHis;}
  TH1F* GetRelDifHist() {return _RelDifHis;}
  TH1F* GetNRebinnedRefHist() {return _NRebinnedRefHis;}
  TH1F* GetNRebinnedNewHist() {return _NRebinnedNewHis;}
  TH1F* GetNRatHist() {return _NRatHis;}
  Color_t GetRefCol() { return _RefCol; }
  Color_t GetNewCol() { return _NewCol; }

  void SetDebug(Bool_t adebug) {_TH1Comp_debug = adebug;}
  void SetRef(TH1F* hRef);
  void SetNew(TH1F* hNew);
  void SetHisStyle(TH1F* h);
  void SetRefCol(Color_t cRef = 1 ) { _RefCol = cRef; }
  void SetNewCol(Color_t cNew = 2 ) { _NewCol = cNew; }
  void SetIsSubset(Bool_t IsSub = kFALSE);
  void Update();

  void DrawBoth(TString option = "", Double_t xmin = 10000, Double_t xmax = -10000, TCanvas *tc = NULL, Color_t RatCol = 51);
  void DrawRef(Option_t *option = ""){ _RefHis->Draw(option); }
  void DrawNew(Option_t *option = ""){ _NewHis->Draw(option); }
  void DrawRat(Option_t *option = "");
  void DrawDif(Option_t *option = "");
  void DrawRelDif(Option_t *option = "");

  ClassDef(TH1Comp,1)// TH1Comp
};

#endif
