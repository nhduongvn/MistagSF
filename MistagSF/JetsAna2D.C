#define JetTree_cxx
#include "JetsAna2D.h" 
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "PS.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include <iostream>
#include <algorithm> // From Mauro

#include <fstream>
#include <TString.h>

#define minPtBin  20.
#define maxPtBin  1000.
#define NPtBin  ( (int) ((maxPtBin-minPtBin)/10.) )

#define minEtaBin  0.
#define maxEtaBin  2.5
#define NEtaBin  ( (int) ((maxEtaBin-minEtaBin)*10.) )

#define NDelimiters 7
#define MAXPU 150

int ReadJSon(unsigned int N_MaxRun, unsigned int N_MaxLumiRangePerRun, unsigned int &RunSize, unsigned int *JS_Run, unsigned int *JS_Run_NLumi, unsigned int *JS_LumiMin, unsigned int *JS_LumiMax, TString fileName)
{
  RunSize = 0;
  unsigned int iLumi = 0;
  int Depth = -1;
  int SubDepth0 = -1;
  int SubDepth2 = -1;
  if( fileName == "one" )
  {
    return 0;
  }
  fileName = "JSonFiles/"+fileName;

  ifstream in;
  in.open(fileName); 
  string line;

  int minI = 1000;

  string Delimiters[NDelimiters] = {"{","}","[","]",":",",","\""};
  int iDelimiterMin = -1; 

  //float tmp;

  if( in.is_open() )
  {
    while ( getline(in, line) )
    {
      while(line.size() > 0)
      {
        minI = 1000;
        int tmpI = 1000;

        for( unsigned int i = 0; i < NDelimiters; ++i )
        {
          tmpI = line.find_first_of(Delimiters[i]);
          if( tmpI != (int) string::npos && tmpI < minI )
          {
            minI = tmpI;
            iDelimiterMin = i;
          }
        }
        if( Depth == 2 )
        {
          if( iLumi >= N_MaxLumiRangePerRun ) return -2;
          if( Delimiters[iDelimiterMin] == "," )
          {
            ++SubDepth2;
            TString tmpLine = line.substr(0,minI );
            JS_LumiMin[RunSize*N_MaxLumiRangePerRun+iLumi] = tmpLine.Atoi();
          }
          if( Delimiters[iDelimiterMin] == "]" )
          {
            ++SubDepth2;
            TString tmpLine = line.substr(0,minI );
            JS_LumiMax[RunSize*N_MaxLumiRangePerRun+iLumi] = tmpLine.Atoi();
            cout << JS_Run[RunSize] <<  "|||" << JS_LumiMin[RunSize*N_MaxLumiRangePerRun+iLumi] << " ; " << JS_LumiMax[RunSize*N_MaxLumiRangePerRun+iLumi] << endl;
            iLumi++;
          }

        } else
        {
          SubDepth2 = -1;
        }

        if( Delimiters[iDelimiterMin] == "{" ) ++Depth;
        if( Delimiters[iDelimiterMin] == "[" ) ++Depth;
        if( Delimiters[iDelimiterMin] == "}" ) --Depth;
        if( Delimiters[iDelimiterMin] == "]" ) --Depth;


        if( Depth == 0 )
        {
          if( Delimiters[iDelimiterMin] == "\"" ) ++SubDepth0;
          if( Delimiters[iDelimiterMin] == ":" ) ++SubDepth0;

          if( SubDepth0 == 1 ) // New Run
          {
            if( RunSize < N_MaxRun ) JS_Run_NLumi[RunSize] = iLumi;
            //if( RunSize >= 0 && RunSize < N_MaxRun ) JS_Run_NLumi[RunSize] = iLumi;
            RunSize++;
            if( RunSize >= N_MaxRun ) return -3;
            iLumi = 0;
            TString tmpRun = line.substr(0,minI);
            JS_Run[RunSize] = tmpRun.Atoi();
          }
        } else
        {
          SubDepth0 = -1;
        }
        
        line.erase(0,minI+1);

        
      }
    }
  }
  else
  {
    cout << "WARNING: Fail to open file " << fileName << endl;
    RunSize = 0;
    return -1;
  }

  in.close();
  cout << "===============" << endl;
  cout << RunSize << endl;
  return RunSize;
}

  // float WeightNS[18*52];
  // float WeightNS[N_ptbin_nseltrack*N_SelTracks];
void ReadNSelTracksWeight(float *WeightNS, unsigned int N_SelTracks, unsigned int N_ptbin_nseltrack,
                          TString fileName, TString *St_fhseltrackPtMin, TString *St_fhseltrackPtMax)
{
  if( fileName == "one" )
  {
    for( unsigned int i = 0; i < N_ptbin_nseltrack; ++i )
    for( unsigned int j = 0; j < N_SelTracks; ++j )
    {
       WeightNS[i*N_SelTracks+j] = 1;
    }
    return;
  }
  fileName = "Weights/"+fileName;

  for( unsigned int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
  {
    ifstream in;
    in.open(fileName+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt]+".txt"); // change name of file

    float tmp;
    unsigned int nlines=0;

    if( in.is_open() )
    {
      while (1) {
        if( nlines >= N_SelTracks )
        {
          in >> tmp; // not used
          if( in.good() )
            cout << "WARNING: file " << fileName << " contains more than " << N_SelTracks << " lines ! " << endl;
          break;
        }
  
        in >> tmp;
        if( !in.good() ) break;
        WeightNS[iPt*N_SelTracks+nlines] = tmp;
        nlines++;
     }
    }
    else cout << "WARNING: Fail to open file " << fileName << endl;

    for( unsigned int i = nlines; i < N_SelTracks; ++i )
    {
       cout << "Warning, file " << fileName<< St_fhseltrackPtMin[iPt]<< "to"<< St_fhseltrackPtMax[iPt]<< ".txt" << " does not contain enough lines ! " << endl;
       cout << "Setting remaining weights to 1 " << endl;
       WeightNS[iPt*N_SelTracks+i] = 1;
    }
    in.close();
  }
}

int ReadPUWeight(float *WeightPU, unsigned int maxPU, TString fileName = "one")
{
  cout << "WARNING WARNING WARNING " << endl;
  cout << "WARNING WARNING WARNING " << endl;
  cout << "WARNING WARNING WARNING " << endl;
  cout << "WARNING WARNING WARNING " << endl;
  cout << "I limit pileup weight to nPV = 80 !!! " << endl;
  int TMPMAXPU = 80;
  if( fileName.Contains("one") )
  {
    for( unsigned int i = 0; i < maxPU; ++i )
    {
       WeightPU[i] = 1;
    }
    return 1;
  }
  fileName = "Weights/"+fileName;
  ifstream in;
  in.open(fileName); // change name of file

  float tmp;
  int nlines=0;

  if( in.is_open() )
  {
    while (1) {
      if( nlines >= (int) maxPU )
      {
        in >> tmp; // not used
        if( in.good() )
          cout << "WARNING: file " << fileName << " contains more than " << maxPU << " lines ! " << endl;
        break;
      }

      in >> tmp; // not used
      if( !in.good() ) break;
      if( tmp - nlines > 0.1 || tmp - nlines < -0.1 )
      {
        cout << "Warning: Index seems not correct, we have tmp = " << tmp << " while nlines = " << nlines << endl;
      }
  
      in >> tmp; // not used
      if( !in.good() ) break;
      WeightPU[nlines] = tmp;
  
      nlines++;
     }
     for( int i = TMPMAXPU; i < nlines; ++i )
     {
       WeightPU[i] = WeightPU[TMPMAXPU];
     }
     return nlines;
  }
  else
  {
    cout << "WARNING: Fail to open file " << fileName << endl;

    for( unsigned int i = nlines; i < maxPU; ++i )
    {
       WeightPU[i] = 0;
    }
    return 0;
  }
}

void ReadpthatWeight(float *pthatMin, float *pthatMax, float *pthatweight, unsigned int maxLines, TString fileName = "one")
{
  if( fileName.Contains("one") )
  {
    for( unsigned int i = 0; i < maxLines; ++i )
    {
       pthatMin[i] = 0;
       pthatMax[i] = 999999;
       pthatweight[i] = 1;
    }
    return;
  }

  fileName = "Weights/"+fileName;
  ifstream in;
  in.open(fileName); // change name of file


  float tmp;
  int nlines=0;

  if( in.is_open() )
  {
    while (1) {
      if( nlines >= (int) maxLines )
      {
        in >> tmp; // not used
        if( in.good() )
          cout << "WARNING: file " << fileName << " contains more than " << maxLines << " lines ! " << endl;
        break;
      }

      in >> tmp; // not used
      if( !in.good() ) break;
      pthatMin[nlines] = tmp;
  
      in >> tmp; // not used
      if( !in.good() ) break;
      pthatMax[nlines] = tmp;
  
      in >> tmp; // not used
      if( !in.good() ) break;
      pthatweight[nlines] = tmp;
  
      nlines++;
     }
  }
  else cout << "WARNING: Fail to open file " << fileName << endl;

  for( unsigned int i = nlines; i < maxLines; ++i )
  {
     pthatMin[i] = 99999;
     pthatMax[i] = 999999;
     pthatweight[i] = 0;
  }
}

void SetSebStyle()
{
 gStyle->SetTitleFillColor(42);
 gStyle->SetTitleFont(1);
 gStyle->SetStatColor(29);
 gStyle->SetCanvasColor(25);   
 gStyle->SetOptStat(1111111);
 gStyle->SetHistFillColor(5);
}

void SaveInFile(TH1* ahisto, TFile* afile)
{
 if (!ahisto) { std::cout << "!! no histo !!" << std::endl; return ;}
 TDirectory* current = gDirectory ;
 afile->cd();
 ahisto->Write();
 current->cd();
}

void JetTree::Loop(int minPV, int maxPV, TString TagName, float aPtMin, float aPtMax, 
                int aIntCut,
                float minCutJetPtMax, float maxCutJetPtMax,
		TString afilename, TString weightPU_file,
                TString weightPthat_file, TString JSONFile, bool truePU, bool WeightTracks, TString TrigType)
{
cout << " Re PUFile : " << weightPU_file << endl;
cout << " Re PthatFile : " << weightPthat_file << endl;

int nLinesPUFile = 0;

TString TagLevel[3] = {"L", "M", "T" };
///////////////////////////////////////////////////////////////////
//                  T I T L E      C A R D S                     //

//$$
  float  TagCut[3]       = {0,0,0};        // tag cut to be determined from TagName and TagLevel
  float  TagCut_CvsL[3] = {0, 0, 0} ;
  float  TagCut_CvsB[3] = {0, 0, 0} ;
  float  PtMin        = aPtMin;   // pt jet min
  float  PtMax        = aPtMax;   // pt jet max
  //float  FreeCut      = aFreeCut; 
  int    IntCut       = aIntCut; 
  TString  filename   = afilename;

// Operating points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X 
  if( TagName == "JP" )
  {	// tagger changed according to Luca talk at https://indico.cern.ch/event/593959/.
    TagCut[0] = .262; //.245;
    TagCut[1] = .5344; //.515;
    TagCut[2] = .7914; //.760;
  }
  else if( TagName == "CSVv2" )
  {	// tagger changed according to Luca mail of 26/01/2018 10h25
    TagCut[0] = 0.5803;	//.460;
    TagCut[1] = 0.8838;	//.800;
    TagCut[2] = 0.9693;	//.935;
  }
  else if( TagName == "DeepCSV" )
  {	// tagger set on 27/07/2021 according to Luca talk at 
        // https://indico.cern.ch/event/1053254/contributions/4425737/attachments/2271561/3857890/BTagPerf_210623_UL16WPUpdate.pdf)
    TagCut[0] = 0.1918;
    TagCut[1] = 0.5847;
    TagCut[2] = 0.8767;
  }
  else if( TagName == "DeepFlavour" ) // "DeepFlavour" = "DeepJet"
  {	// tagger set on 27/07/2021 according to Luca talk at 
        // https://indico.cern.ch/event/1053254/contributions/4425737/attachments/2271561/3857890/BTagPerf_210623_UL16WPUpdate.pdf)
    TagCut[0] = 0.0480;
    TagCut[1] = 0.2489;
    TagCut[2] = 0.6377;
  }
  else if( TagName == "cMVAv2" )
  {	// tagger changed according to Luca talk at https://indico.cern.ch/event/593959/.
    TagCut[0] = -0.5884;	//-0.715;
    TagCut[1] =  0.4432;	//0.185;
    TagCut[2] =  0.9432;	//0.875;
  }
  else if ( TagName == "CvsB" || TagName == "CvsL") //CTag is 2D cut of CvsB and CvsL
  {
    //From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    TagCut_CvsL[0] = -0.53 ;
    TagCut_CvsL[1] = 0.07 ;
    TagCut_CvsL[2] = 0.87 ;
    TagCut_CvsB[0] = -0.26 ; 
    TagCut_CvsB[1] = -0.10 ; 
    TagCut_CvsB[2] = -0.30 ; 
  }
  else if ( TagName == "DeepCSVCvsB" || TagName == "DeepCSVCvsL") //CTag is 2D cut of CvsB and CvsL
  {	// tagger set on 29/10/2020 according to Seth talk at 
        // https://indico.cern.ch/event/967689/contributions/4083041/attachments/2130779/3592013/CTagger_WP_UltraLegacy2018_DeepCSV_DeepJet.pdf
    TagCut_CvsL[0] = 0.064 ;
    TagCut_CvsL[1] = 0.153 ;
    TagCut_CvsL[2] = 0.405 ;
    TagCut_CvsB[0] = 0.313 ; 
    TagCut_CvsB[1] = 0.363 ; 
    TagCut_CvsB[2] = 0.288 ; 
  }
  else if ( TagName == "DeepFlavourCvsB" || TagName == "DeepFlavourCvsL") //CTag is 2D cut of CvsB and CvsL
  {	// tagger set on 29/10/2020 according to Seth talk at 
        // https://indico.cern.ch/event/967689/contributions/4083041/attachments/2130779/3592013/CTagger_WP_UltraLegacy2018_DeepCSV_DeepJet.pdf
    TagCut_CvsL[0] = 0.038 ;
    TagCut_CvsL[1] = 0.099 ;
    TagCut_CvsL[2] = 0.282 ;
    TagCut_CvsB[0] = 0.246 ; 
    TagCut_CvsB[1] = 0.325 ; 
    TagCut_CvsB[2] = 0.267 ; 
  }
  else
  {
    cout << "Undefined TagName for " << TagName << endl;
    return;
  }
//$$
  int NtrackMin = 1;         // Taggability
//
//  int NtrackMin = int(FreeCut);
//
//  float VetoPos = 4.;
//  float VetoPos = FreeCut;
  //float VetoPos = 1000.;
//$$

  //int NpuMin = 0;
//$$
//$$  NpuMin = int(FreeCut);
//$$

//$$
  //bool GluonSplit = false;
  //bool ReweightNSelTracks = false; // To activate the reweighting according to the nSelTracks.
  //TString weightNS_file = "NSEL_ak5Jets/nselWeight_pt_"; // Set to "one" to set all weights to 1.
  TString weightNS_file = "one"; // Set to "one" to set all weights to 1.

  bool BTagMu  = false;
  //bool MuonJet = false;  // need to have MuOp    = false
  //bool MuOp    = false;  // need to have MuonJet = false and Away = false
  //bool Away    = false;
  //bool EtaBin  = false;

  int Year = 2011;
  if( TrigType == "2011" ) Year = 2011;
  else if( TrigType == "2012" || TrigType == "2015" || TrigType == "2015_Cert" || TrigType == "2015_LooserSel"
        || TrigType == "2016_BX_4X2" || TrigType == "2016_BX_4X3" || TrigType == "2016_BX_4X4" || TrigType == "2016_BX_4X5"
        || TrigType == "2016_BX_4X6" || TrigType == "2016_BX_Cont" || TrigType == "2016_All"
        || TrigType == "2016_flat1"	// To provide SF for the flat period 1 (Run in 274300 to 276097) defined in mail of 13/7/2016
        || TrigType == "2016_flat2"	// To provide SF for the flat period 2 (Run in 275784 to 276097) defined in mail of 13/7/2016
        || TrigType == "2016_76fb_L1T"  // To provide SF for the whole 7.6fb-1 period with L1T cert as defined in mail of 13/7/2016
         ) Year = 2012;
  else if( TrigType == "Moriond_v1_2016") Year = 2016;
  else if( TrigType == "Run2016UL") Year = 2016;
  else if( TrigType.Contains("Run2016") ) Year = 2016;
  else if( TrigType.Contains("Run2017") ) Year = 2016;
  //else if( TrigType.Contains("Run2017")&&TrigType.Contains("-17Nov2017-v1") ) Year = 2016;
  else if( TrigType.Contains("Run2018") ) Year = 2016;
  else
  {
    cout << "This TrigType is not foreseen ( " << TrigType << ")" << endl;
  }

  //bool AlaLuca = false;
  //bool AlaPablo = false;

  //bool GSplitMin = false;
  //bool GSplitMax = false;

  //bool BFragMin = false;
  //bool BFragMax = false;

  //bool GenTrkMin = false;
  //bool GenTrkMax = false;

  //bool CFrag = false;
  //bool V0weight = false;

  //bool dRJetCut = false;

  bool TopAna = false;
//$$

///////////////////////////////////////////////////////////////////
 
//$$
  TH1D::SetDefaultSumw2(kTRUE);
//$$

//**********************************
// Data
//**********************************
  TH1D* hBeforeTrig_MaxJetPt       = new TH1D("hBeforeTrig_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig30_MaxJetPt       = new TH1D("hAfterTrig30_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig60_MaxJetPt       = new TH1D("hAfterTrig60_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig80_MaxJetPt       = new TH1D("hAfterTrig80_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig110_MaxJetPt       = new TH1D("hAfterTrig110_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig150_MaxJetPt       = new TH1D("hAfterTrig150_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig190_MaxJetPt       = new TH1D("hAfterTrig190_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig240_MaxJetPt       = new TH1D("hAfterTrig240_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hAfterTrig300_MaxJetPt       = new TH1D("hAfterTrig300_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hAfterTrig370_MaxJetPt       = new TH1D("hAfterTrig370_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);

  TH1D* hData_All_nPV         = new TH1D("hData_All_nPV","nb. of PV",MAXPU,0.5,MAXPU+0.5);
  TH1D* hData_All_NJets       = new TH1D("hData_All_NJets","nb. of jets",16,-0.5,15.5);
  TH1D* hData_All_tracks      = new TH1D("hData_All_tracks","nb. of tracks",20,0.5,20.5);
  TH1D* hData_All_JetPV       = new TH1D("hData_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hData_All_JetPt       = new TH2D("hData_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_All_JetEta      = new TH1D("hData_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_All_JetRun      = new TH1D("hData_All_JetRun","Run number",400,190700.5,191100.5);
  TH1D* hData_nPV             = new TH1D("hData_nPV","nb. of PV",MAXPU,0.5,MAXPU+0.5);
  TH1D* hData_NJets           = new TH1D("hData_NJets","nb. of jets",16,-0.5,15.5);
  TH1D* hData_NJet30          = new TH1D("hData_NJet30","nb. jets pt>30",16,-0.5,15.5);
  TH1D* hData_JetPV           = new TH1D("hData_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hData_JetPt           = new TH2D("hData_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_JetEta          = new TH1D("hData_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_JetRun          = new TH1D("hData_JetRun","Run number",400,190700.5,191100.5);
  TH2D* hData_Tagger           = new TH2D("hData_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hData_PTagger           = new TH2D("hData_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hData_MI2_Tagger           = new TH2D("hData_MI2_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hData_MI2_PTagger           = new TH2D("hData_MI2_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_JetPt_LE5pv    = new TH1D("hData_JetPt_LE5pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hData_JetEta_LE5pv   = new TH1D("hData_JetEta_LE5pv","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_JetPt_GE6pv    = new TH1D("hData_JetPt_GE6pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hData_JetEta_GE6pv   = new TH1D("hData_JetEta_GE6pv","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hData_1JetPt           = new TH1D("hData_1JetPt","pt(jet)",30,20.,320.);
  TH1D* hData_2JetPt           = new TH1D("hData_2JetPt","pt(jet)",30,20.,320.);
  TH1D* hData_3JetPt           = new TH1D("hData_3JetPt","pt(jet)",30,20.,320.);
  TH1D* hData_4JetPt           = new TH1D("hData_4JetPt","pt(jet)",30,20.,320.);
  TH1D* hData_PosTag           = new TH1D("hData_PosTag","Tag(+)",50,0.,25.);
  TH1D* hData_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hData_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hData_PosTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hData_PosTag_JetRun[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hData_PosTag_JetPV[i]    = new TH1D("hData_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hData_PosTag_JetPt[i]    = new TH2D("hData_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hData_PosTag_JetEta[i]   = new TH1D("hData_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hData_PosTag_JetRun[i]   = new TH1D("hData_PosTag"+TagLevel[i]+"_JetRun","Run number",70,190400,191100);
  }

  TH1D* hData_NegTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hData_NegTag_JetPt[3] = {NULL, NULL, NULL};    
  TH1D* hData_NegTag_JetEta[3] = {NULL, NULL, NULL};   
  TH1D* hData_NegTag_JetRun[3] = {NULL, NULL, NULL};   
  TH1D* hData_NegTag_JetPt_LE5pv[3] = {NULL, NULL, NULL};  
  TH1D* hData_NegTag_JetEta_LE5pv[3] = {NULL, NULL, NULL}; 
  TH1D* hData_NegTag_JetPt_GE6pv[3] = {NULL, NULL, NULL};  
  TH1D* hData_NegTag_JetEta_GE6pv[3] = {NULL, NULL, NULL}; 
  TH1D* hData_NegTag_1JetPt[3] = {NULL, NULL, NULL};    
  TH1D* hData_NegTag_2JetPt[3] = {NULL, NULL, NULL};    
  TH1D* hData_NegTag_3JetPt[3] = {NULL, NULL, NULL};    
  TH1D* hData_NegTag_4JetPt[3] = {NULL, NULL, NULL};    

  for(unsigned int i = 0; i < 3; ++i )
  {
    hData_NegTag_JetPV[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hData_NegTag_JetPt[i] = new TH2D("hData_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hData_NegTag_JetEta[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hData_NegTag_JetRun[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetRun","Run number",70,190400,191100);

    hData_NegTag_JetPt_LE5pv[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetPt_LE5pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
    hData_NegTag_JetEta_LE5pv[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetEta_LE5pv","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hData_NegTag_JetPt_GE6pv[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetPt_GE6pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
    hData_NegTag_JetEta_GE6pv[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_JetEta_GE6pv","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);


    hData_NegTag_1JetPt[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_1JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_2JetPt[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_2JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_3JetPt[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_3JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_4JetPt[i] = new TH1D("hData_NegTag"+TagLevel[i]+"_4JetPt","pt(jet)",30,20.,320.);
  }

  //TH1D* hData_Trigger     = new TH1D("hData_Trigger","Trigger",100,-0.5,99.5);
  //TH1D* hData_Trig1000    = new TH1D("hData_Trig1000","Trig1000",1000,-0.5,999.5);
  //TH1D* hData_Trig1_JetPt = new TH1D("hData_Trig1_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hData_Trig2_JetPt = new TH1D("hData_Trig2_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hData_Trig3_JetPt = new TH1D("hData_Trig3_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hData_Trig4_JetPt = new TH1D("hData_Trig4_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hData_Trig5_JetPt = new TH1D("hData_Trig5_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  //TH1D* hData_Trig6_JetPt = new TH1D("hData_Trig6_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);

//**********************************
// All flavours in Monte Carlo
//**********************************
  TH1D* hAllFlav_pthatAll    = new TH1D("hAllFlav_pthatAll","pthat",80,0.,800.);
  TH1D* hAllFlav_pthat       = new TH1D("hAllFlav_pthat","pthat",80,0.,800.);
  TH1D* hAllFlav_pthatTrig   = new TH1D("hAllFlav_pthatTrig","pthat",80,0.,800.);
  TH1D* hAllFlav_All_nPU     = new TH1D("hAllFlav_All_nPU","nb. of PU",60,0.,60.);
  TH1D* hAllFlav_nPU         = new TH1D("hAllFlav_nPU","nb. of PU",60,0.,60.);
  TH1D* hAllFlav_All_Flavour = new TH1D("hAllFlav_All_Flavour","Flavour",22,-0.5,21.5);
  TH1D* hAllFlav_Flavour     = new TH1D("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
  TH1D* hAllFlav_Tagger_Bwd  = new TH1D("hAllFlav_Tagger_Bwd","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Cwd  = new TH1D("hAllFlav_Tagger_Cwd","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Tau  = new TH1D("hAllFlav_Tagger_Tau","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Gam  = new TH1D("hAllFlav_Tagger_Gam","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_K0s  = new TH1D("hAllFlav_Tagger_K0s","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Lam  = new TH1D("hAllFlav_Tagger_Lam","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Int  = new TH1D("hAllFlav_Tagger_Int","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Fak  = new TH1D("hAllFlav_Tagger_Fak","Tagger",100,-25.,25.);
  TH1D* hAllFlav_Tagger_Oth  = new TH1D("hAllFlav_Tagger_Oth","Tagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Bwd  = new TH1D("hAllFlav_PTagger_Bwd","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Cwd  = new TH1D("hAllFlav_PTagger_Cwd","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Tau  = new TH1D("hAllFlav_PTagger_Tau","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Gam  = new TH1D("hAllFlav_PTagger_Gam","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_K0s  = new TH1D("hAllFlav_PTagger_K0s","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Lam  = new TH1D("hAllFlav_PTagger_Lam","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Int  = new TH1D("hAllFlav_PTagger_Int","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Fak  = new TH1D("hAllFlav_PTagger_Fak","PTagger",100,-25.,25.);
  TH1D* hAllFlav_PTagger_Oth  = new TH1D("hAllFlav_PTagger_Oth","PTagger",100,-25.,25.);

  TH1D* hAllFlav_Gam_JetPV  	   = new TH1D("hAllFlav_Gam_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hAllFlav_Gam_JetPt	   = new TH2D("hAllFlav_Gam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_Gam_JetEta	   = new TH1D("hAllFlav_Gam_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_K0s_JetPV  	   = new TH1D("hAllFlav_K0s_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hAllFlav_K0s_JetPt	   = new TH2D("hAllFlav_K0s_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_K0s_JetEta	   = new TH1D("hAllFlav_K0s_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);

  TH1D* hAllFlav_Gam_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_Gam_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Gam_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_K0s_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_K0s_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_K0s_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hAllFlav_Gam_NegTag_JetPV[i]   = new TH1D("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_Gam_NegTag_JetPt[i]   = new TH2D("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Gam_NegTag_JetEta[i]  = new TH1D("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_K0s_NegTag_JetPV[i]   = new TH1D("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_K0s_NegTag_JetPt[i]   = new TH2D("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_K0s_NegTag_JetEta[i]  = new TH1D("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }
//
  TH1D* hAllFlav_K0s_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_K0s_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_K0s_PosTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Lam_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_Lam_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Lam_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hAllFlav_K0s_PosTag_JetPV[i]   = new TH1D("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_K0s_PosTag_JetPt[i]   = new TH2D("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_K0s_PosTag_JetEta[i]  = new TH1D("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Lam_PosTag_JetPV[i]   = new TH1D("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_Lam_PosTag_JetPt[i]   = new TH2D("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Lam_PosTag_JetEta[i]  = new TH1D("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }
//
  TH1D* hAllFlav_Lam_JetPV  	   = new TH1D("hAllFlav_Lam_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hAllFlav_Lam_JetPt	   = new TH2D("hAllFlav_Lam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_Lam_JetEta	   = new TH1D("hAllFlav_Lam_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_Fak_JetPV  	   = new TH1D("hAllFlav_Fak_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hAllFlav_Fak_JetPt	   = new TH2D("hAllFlav_Fak_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hAllFlav_Fak_JetEta	   = new TH1D("hAllFlav_Fak_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);

  TH1D* hAllFlav_Lam_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_Lam_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Lam_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Fak_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hAllFlav_Fak_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hAllFlav_Fak_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hAllFlav_Lam_NegTag_JetPV[i]   = new TH1D("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_Lam_NegTag_JetPt[i]   = new TH2D("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Lam_NegTag_JetEta[i]  = new TH1D("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Fak_NegTag_JetPV[i]   = new TH1D("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hAllFlav_Fak_NegTag_JetPt[i]   = new TH2D("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Fak_NegTag_JetEta[i]  = new TH1D("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

//**********************************
// udsg-jets
//**********************************
  TH1D* hLightFlav_All_JetPV      = new TH1D("hLightFlav_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_All_JetPt      = new TH2D("hLightFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_All_JetEta     = new TH1D("hLightFlav_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_JetPU          = new TH1D("hLightFlav_JetPU","#PU",MAXPU,0.5,MAXPU+0.5);
  TH1D* hLightFlav_JetPV          = new TH1D("hLightFlav_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_JetPt          = new TH2D("hLightFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_JetEta         = new TH1D("hLightFlav_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hLightFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hLightFlav_NegTag_JetPV[i]   = new TH1D("hLightFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_NegTag_JetPt[i]   = new TH2D("hLightFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_NegTag_JetEta[i]  = new TH1D("hLightFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

  TH2D* hLightFlav_JetPt_0pu      = new TH2D("hLightFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hLightFlav_JetPt_GE8pu    = new TH2D("hLightFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hLightFlav_Tagger         = new TH2D("hLightFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hLightFlav_MI2_Tagger           = new TH2D("hLightFlav_MI2_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hLightFlav_MI2_PTagger           = new TH2D("hLightFlav_MI2_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);

  TH1D* hLightFlav_Tagger_Gam     = new TH1D("hLightFlav_Tagger_Gam","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_K0s     = new TH1D("hLightFlav_Tagger_K0s","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Lam     = new TH1D("hLightFlav_Tagger_Lam","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Bwd     = new TH1D("hLightFlav_Tagger_Bwd","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Cwd     = new TH1D("hLightFlav_Tagger_Cwd","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Tau     = new TH1D("hLightFlav_Tagger_Tau","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Int     = new TH1D("hLightFlav_Tagger_Int","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Fak     = new TH1D("hLightFlav_Tagger_Fak","Tagger",100,-25.,25.);
  TH1D* hLightFlav_Tagger_Oth     = new TH1D("hLightFlav_Tagger_Oth","Tagger",100,-25.,25.);
  TH2D* hLightFlav_PTagger         = new TH2D("hLightFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_PTagger_Gam     = new TH1D("hLightFlav_PTagger_Gam","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_K0s     = new TH1D("hLightFlav_PTagger_K0s","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Lam     = new TH1D("hLightFlav_PTagger_Lam","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Bwd     = new TH1D("hLightFlav_PTagger_Bwd","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Cwd     = new TH1D("hLightFlav_PTagger_Cwd","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Tau     = new TH1D("hLightFlav_PTagger_Tau","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Int     = new TH1D("hLightFlav_PTagger_Int","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Fak     = new TH1D("hLightFlav_PTagger_Fak","PTagger",100,-25.,25.);
  TH1D* hLightFlav_PTagger_Oth     = new TH1D("hLightFlav_PTagger_Oth","PTagger",100,-25.,25.);
  TH1D* hLightFlav_1JetPt         = new TH1D("hLightFlav_1JetPt","pt(jet)",30,20.,320.);
  TH1D* hLightFlav_2JetPt         = new TH1D("hLightFlav_2JetPt","pt(jet)",30,20.,320.);
  TH1D* hLightFlav_3JetPt         = new TH1D("hLightFlav_3JetPt","pt(jet)",30,20.,320.);
  TH1D* hLightFlav_4JetPt         = new TH1D("hLightFlav_4JetPt","pt(jet)",30,20.,320.);

  TH1D* hLightFlav_PosTag_JetPU[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_PosTagger      = new TH1D("hLightFlav_PosTagger","Tagger",1000,0.,25.);
  TH1D* hLightFlav_PosTagger_0pu  = new TH1D("hLightFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1D* hLightFlav_PosTagger_GE8pu= new TH1D("hLightFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1D* hLightFlav_PosTag_1JetPt[3] = {NULL, NULL, NULL};  
  TH1D* hLightFlav_PosTag_2JetPt[3] = {NULL, NULL, NULL};  
  TH1D* hLightFlav_PosTag_3JetPt[3] = {NULL, NULL, NULL};  
  TH1D* hLightFlav_PosTag_4JetPt[3] = {NULL, NULL, NULL};  
  TH1D* hLightFlav_BCT_PosTag_JetPV[3] = {NULL, NULL, NULL};   
  TH2D* hLightFlav_BCT_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_BCT_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Gam_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_Gam_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Gam_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_K0s_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_K0s_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_K0s_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Lam_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_Lam_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Lam_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Fak_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_Fak_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Fak_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Fak_PosTag_pthat[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Fak_PosTag_nPU[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Oth_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2D* hLightFlav_Oth_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Oth_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Oth_PosTag_pthat[3] = {NULL, NULL, NULL};
  TH1D* hLightFlav_Oth_PosTag_nPU[3] = {NULL, NULL, NULL};

  for( unsigned int i = 0; i < 3; ++i )
  {
    hLightFlav_PosTag_JetPU[i]   = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_JetPU","#PU",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_PosTag_JetPV[i]   = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_PosTag_JetPt[i]   = new TH2D("hLightFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_PosTag_JetEta[i]  = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_PosTag_1JetPt[i]  = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_1JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_2JetPt[i]  = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_2JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_3JetPt[i]  = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_3JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_4JetPt[i]  = new TH1D("hLightFlav_PosTag"+TagLevel[i]+"_4JetPt","pt(jet)",30,20.,320.);
    hLightFlav_BCT_PosTag_JetPV[i]   = new TH1D("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_BCT_PosTag_JetPt[i]   = new TH2D("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_BCT_PosTag_JetEta[i]  = new TH1D("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Gam_PosTag_JetPV[i]   = new TH1D("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Gam_PosTag_JetPt[i]   = new TH2D("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Gam_PosTag_JetEta[i]  = new TH1D("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_K0s_PosTag_JetPV[i]   = new TH1D("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_K0s_PosTag_JetPt[i]   = new TH2D("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_K0s_PosTag_JetEta[i]  = new TH1D("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Lam_PosTag_JetPV[i]   = new TH1D("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Lam_PosTag_JetPt[i]   = new TH2D("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Lam_PosTag_JetEta[i]  = new TH1D("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_PosTag_JetPV[i]   = new TH1D("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Fak_PosTag_JetPt[i]   = new TH2D("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_PosTag_JetEta[i]  = new TH1D("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_PosTag_pthat[i]  = new TH1D("hLightFlav_Fak_PosTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Fak_PosTag_nPU[i]  = new TH1D("hLightFlav_Fak_PosTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
    hLightFlav_Oth_PosTag_JetPV[i]   = new TH1D("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Oth_PosTag_JetPt[i]   = new TH2D("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_PosTag_JetEta[i]  = new TH1D("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_PosTag_pthat[i]  = new TH1D("hLightFlav_Oth_PosTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Oth_PosTag_nPU[i]  = new TH1D("hLightFlav_Oth_PosTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
  }


  TH1D* hLightFlav_Gam_JetPV  	     = new TH1D("hLightFlav_Gam_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_Gam_JetPt	     = new TH2D("hLightFlav_Gam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Gam_JetEta	     = new TH1D("hLightFlav_Gam_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_K0s_JetPV  	     = new TH1D("hLightFlav_K0s_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_K0s_JetPt	     = new TH2D("hLightFlav_K0s_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_K0s_JetEta	     = new TH1D("hLightFlav_K0s_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Lam_JetPV  	     = new TH1D("hLightFlav_Lam_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_Lam_JetPt	     = new TH2D("hLightFlav_Lam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Lam_JetEta	     = new TH1D("hLightFlav_Lam_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Fak_JetPV  	     = new TH1D("hLightFlav_Fak_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_Fak_JetPt	     = new TH2D("hLightFlav_Fak_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Fak_JetEta	     = new TH1D("hLightFlav_Fak_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Oth_JetPV  	     = new TH1D("hLightFlav_Oth_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hLightFlav_Oth_JetPt	     = new TH2D("hLightFlav_Oth_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hLightFlav_Oth_JetEta	     = new TH1D("hLightFlav_Oth_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);

  TH1D* hLightFlav_Fak_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hLightFlav_Fak_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_Fak_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_Oth_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hLightFlav_Oth_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_Oth_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_Oth_NegTag_pthat[3] = {NULL, NULL, NULL };
  TH1D* hLightFlav_Oth_NegTag_nPU[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hLightFlav_Fak_NegTag_JetPV[i]   = new TH1D("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Fak_NegTag_JetPt[i]   = new TH2D("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_NegTag_JetEta[i]  = new TH1D("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_NegTag_JetPV[i]   = new TH1D("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hLightFlav_Oth_NegTag_JetPt[i]   = new TH2D("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_NegTag_JetEta[i]  = new TH1D("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_NegTag_pthat[i]  = new TH1D("hLightFlav_Oth_NegTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Oth_NegTag_nPU[i]  = new TH1D("hLightFlav_Oth_NegTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
  }

//**********************************
// gluon-jets
//**********************************
  TH1D* hGluonFlav_All_JetPV      = new TH1D("hGluonFlav_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hGluonFlav_All_JetPt      = new TH2D("hGluonFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hGluonFlav_All_JetEta     = new TH1D("hGluonFlav_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hGluonFlav_JetPV          = new TH1D("hGluonFlav_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hGluonFlav_JetPt          = new TH2D("hGluonFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hGluonFlav_JetEta         = new TH1D("hGluonFlav_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hGluonFlav_Tagger         = new TH2D("hGluonFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hGluonFlav_PTagger         = new TH2D("hGluonFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hGluonFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hGluonFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hGluonFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hGluonFlav_PosTag = new TH1D("hGluonFlav_PosTag","Tag(+)",50,0.,25.);
  TH1D* hGluonFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hGluonFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hGluonFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hGluonFlav_NegTag_JetPV[i]   = new TH1D("hGluonFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hGluonFlav_NegTag_JetPt[i]   = new TH2D("hGluonFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hGluonFlav_NegTag_JetEta[i]  = new TH1D("hGluonFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);

    hGluonFlav_PosTag_JetPV[i]   = new TH1D("hGluonFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hGluonFlav_PosTag_JetPt[i]   = new TH2D("hGluonFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hGluonFlav_PosTag_JetEta[i]  = new TH1D("hGluonFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

//**********************************
// uds-jets
//**********************************
  TH1D* hUDSFlav_All_JetPV      = new TH1D("hUDSFlav_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hUDSFlav_All_JetPt      = new TH2D("hUDSFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hUDSFlav_All_JetEta     = new TH1D("hUDSFlav_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hUDSFlav_JetPV          = new TH1D("hUDSFlav_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hUDSFlav_JetPt          = new TH2D("hUDSFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hUDSFlav_JetEta         = new TH1D("hUDSFlav_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hUDSFlav_Tagger         = new TH2D("hUDSFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hUDSFlav_PTagger         = new TH2D("hUDSFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hUDSFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hUDSFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hUDSFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hUDSFlav_PosTag = new TH1D("hUDSFlav_PosTag","Tag(+)",50,0.,25.);
  TH1D* hUDSFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hUDSFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hUDSFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hUDSFlav_NegTag_JetPV[i]   = new TH1D("hUDSFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hUDSFlav_NegTag_JetPt[i]   = new TH2D("hUDSFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hUDSFlav_NegTag_JetEta[i]  = new TH1D("hUDSFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hUDSFlav_PosTag_JetPV[i]   = new TH1D("hUDSFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hUDSFlav_PosTag_JetPt[i]   = new TH2D("hUDSFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hUDSFlav_PosTag_JetEta[i]  = new TH1D("hUDSFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

//**********************************
// c-jets
//**********************************
  TH1D* hCFlav_All_JetPV      = new TH1D("hCFlav_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hCFlav_All_JetPt      = new TH2D("hCFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hCFlav_All_JetEta     = new TH1D("hCFlav_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hCFlav_JetPV          = new TH1D("hCFlav_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hCFlav_JetPt          = new TH2D("hCFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hCFlav_JetEta         = new TH1D("hCFlav_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hCFlav_Tagger         = new TH2D("hCFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hCFlav_PTagger         = new TH2D("hCFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hCFlav_MI2_Tagger           = new TH2D("hCFlav_MI2_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hCFlav_MI2_PTagger           = new TH2D("hCFlav_MI2_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  //TH1D* hCFlav_deltaR         = new TH1D("hCFlav_deltaR","#Delta R gluon split",50,0.,1.);
  TH2D* hCFlav_JetPt_0pu      = new TH2D("hCFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hCFlav_JetPt_GE8pu    = new TH2D("hCFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hCFlav_PosTagger      = new TH1D("hCFlav_PosTagger","Tagger",1000,0.,25.);
  TH1D* hCFlav_PosTagger_0pu  = new TH1D("hCFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1D* hCFlav_PosTagger_GE8pu= new TH1D("hCFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1D* hCFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hCFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hCFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1D* hCFlav_PosTag = new TH1D("hCFlav_PosTag","Tag(+)",50,0.,25.);
  TH1D* hCFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hCFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hCFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hCFlav_NegTag_JetPV[i]   = new TH1D("hCFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hCFlav_NegTag_JetPt[i]   = new TH2D("hCFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hCFlav_NegTag_JetEta[i]  = new TH1D("hCFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);

    hCFlav_PosTag_JetPV[i]   = new TH1D("hCFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hCFlav_PosTag_JetPt[i]   = new TH2D("hCFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hCFlav_PosTag_JetEta[i]  = new TH1D("hCFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

//**********************************
// b-jets
//**********************************
  TH1D* hBFlav_All_JetPV      = new TH1D("hBFlav_All_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hBFlav_All_JetPt      = new TH2D("hBFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hBFlav_All_JetEta     = new TH1D("hBFlav_All_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hBFlav_JetPU          = new TH1D("hBFlav_JetPU","#PU",MAXPU,0.5,MAXPU+0.5);
  TH1D* hBFlav_JetPV          = new TH1D("hBFlav_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
  TH2D* hBFlav_JetPt          = new TH2D("hBFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hBFlav_JetPt_etaLT12  = new TH1D("hBFlav_JetPt_etaLT12","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hBFlav_JetPt_etaGT12  = new TH1D("hBFlav_JetPt_etaGT12","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1D* hBFlav_JetEta         = new TH1D("hBFlav_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hBFlav_PosTagger      = new TH1D("hBFlav_PosTagger","Tagger",1000,0.,25.);
  TH2D* hBFlav_Tagger         = new TH2D("hBFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hBFlav_PTagger         = new TH2D("hBFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hBFlav_MI2_Tagger           = new TH2D("hBFlav_MI2_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hBFlav_MI2_PTagger           = new TH2D("hBFlav_MI2_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  //TH1D* hBFlav_deltaR         = new TH1D("hBFlav_deltaR","#Delta R gluon split",50,0.,1.);
  TH2D* hBFlav_JetPt_0pu      = new TH2D("hBFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hBFlav_JetPt_GE8pu    = new TH2D("hBFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1D* hBFlav_PosTagger_0pu  = new TH1D("hBFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1D* hBFlav_PosTagger_GE8pu= new TH1D("hBFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1D* hBFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hBFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hBFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };

  TH1D* hBFlav_PosTag= new TH1D("hBFlav_PosTag","Tag(+)",50,0.,25.);
  TH1D* hBFlav_PosTag_JetPU[3] = {NULL, NULL, NULL };
  TH1D* hBFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2D* hBFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1D* hBFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( unsigned int i = 0; i < 3; ++i )
  {
    hBFlav_NegTag_JetPV[i]   = new TH1D("hBFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hBFlav_NegTag_JetPt[i]   = new TH2D("hBFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hBFlav_NegTag_JetEta[i]  = new TH1D("hBFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
    hBFlav_PosTag_JetPU[i]   = new TH1D("hBFlav_PosTag"+TagLevel[i]+"_JetPU","#PU",MAXPU,0.5,MAXPU+0.5);
    hBFlav_PosTag_JetPV[i]   = new TH1D("hBFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",MAXPU,0.5,MAXPU+0.5);
    hBFlav_PosTag_JetPt[i]   = new TH2D("hBFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hBFlav_PosTag_JetEta[i]  = new TH1D("hBFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",NEtaBin, minEtaBin, maxEtaBin);
  }

//**********************************
// No flavour
//**********************************
  TH2D* hNoFlav_Tagger        = new TH2D("hNoFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hNoFlav_PTagger        = new TH2D("hNoFlav_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);

  TH2D* hNoFlav_MI2_Tagger           = new TH2D("hNoFlav_MI2_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH2D* hNoFlav_MI2_PTagger           = new TH2D("hNoFlav_MI2_PTagger","PTagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);

//--------- Start use of vector of histos
  //-1° Eta binning
  unsigned int NhEta = 3; 
  float fhEtaMax[3] = { 0.8, 1.6, 2.5 };
  float fhEtaMin[3];
  TString St_hEtaMin[3];
  TString St_hEtaMax[3];
  for(unsigned int i = 0; i < NhEta; ++i)
  {
    if( i == 0 )
    {
      St_hEtaMin[i] = Form("%1.1f",0.0);
      St_hEtaMin[i].Replace(1,1,"_");
      fhEtaMin[i] = 0.;
    }
    else
    {
      St_hEtaMin[i] = Form("%1.1f",fhEtaMax[i-1]);
      St_hEtaMin[i].Replace(1,1,"_");
      fhEtaMin[i] = fhEtaMax[i-1];
    }
    St_hEtaMax[i] = Form("%1.1f",fhEtaMax[i]);
    St_hEtaMax[i].Replace(1,1,"_");
  } 
  //-2° Pt binning
  /* For december 2013 comparison with ttH
  int NhPt = 6;
  // float fhPtMax[6] = { 30, 40, 60,100,160, 1000}; // For december 2013 comparison with ttH
  float fhPtMin[6];
  */
  // For december 2015 comparison with ttH
  unsigned int NhPt = 4;
  float fhPtMax[4] = { 30, 40, 60,1000}; // For december 2013 comparison with ttH
  float fhPtMin[4];
  TString St_hPtMin[4];
  TString St_hPtMax[4];
  for(unsigned int i = 0; i < NhPt; ++i)
  {
    if( i == 0 )
    {
      St_hPtMin[i] = Form("%4.0f",0.0);
      St_hPtMin[i].ReplaceAll(" ","0");
      fhPtMin[i] = 20;
    }
    else
    {
      St_hPtMin[i] = Form("%4.0f",fhPtMax[i-1]);
      St_hPtMin[i].ReplaceAll(" ","0");
      fhPtMin[i] = fhPtMax[i-1];
    }
    St_hPtMax[i] = Form("%4.0f",fhPtMax[i]);
    St_hPtMax[i].ReplaceAll(" ","0");
  } 
  //-3° Histograms Pt binning
  TH1D *hData_JetPtBinned_nPV[4];
  //-4° Histograms for Eta and Pt binning
  TH1D *h_AllFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_UDSFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_CFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_BFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_GluonFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_K0s_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_K0s_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Lam_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Lam_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Fak_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Fak_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Gam_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Gam_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Oth_IntTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Oth_IntTag[3][4]; //[NhEta][NhPt];

  TH1D *h_AllFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_UDSFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_CFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_BFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_GluonFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_K0s_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_K0s_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Lam_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Lam_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Fak_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Fak_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Gam_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Gam_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Oth_NegTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Oth_NegTag[3][4]; //[NhEta][NhPt];

  TH1D *h_AllFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_UDSFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_CFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_BFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_GluonFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_K0s_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_K0s_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Lam_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Lam_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Fak_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Fak_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Gam_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Gam_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_AllFlavour_Oth_PosTag[3][4]; //[NhEta][NhPt];
  TH1D *h_LightFlavour_Oth_PosTag[3][4]; //[NhEta][NhPt];
/*
double csvbins[] = { -10.0, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };
TH1D *h = new TH1D("h","test",17,csvbins);
*/
  for( unsigned int iPt = 0; iPt < NhPt; ++iPt )
  {
    TString hName = "hData_nPV_PtMax_";
    hName += St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
    hData_JetPtBinned_nPV[iPt] = new TH1D(hName,hName,MAXPU,0.5,MAXPU+0.5);
  }

  unsigned int Ncsvbins = 17;
  double csvbins[18] = { -1.0, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };
  for( unsigned int i = 0; i < Ncsvbins+1; ++i ) csvbins[i]*=25.;

  for( unsigned int iEta = 0; iEta < NhEta; ++iEta )
  {
    for( unsigned int iPt = 0; iPt < NhPt; ++iPt )
    {
      TString hName = "IntTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "IntTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "IntTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_IntTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      //-------
      hName = "NegTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "NegTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "NegTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_NegTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      //-------
      hName = "PosTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "PosTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_K0sPart_Eta"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);

      hName = "PosTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_PosTag[iEta][iPt] = new TH1D(hName,hName,Ncsvbins,csvbins);
    }
  }

  // Pt Binning for nselTrack reweighting
  unsigned int N_ptbin_nseltrack = 18;
  float fhseltrackPtMin[18] = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000};
  float fhseltrackPtMax[18];
  TString St_fhseltrackPtMin[18];
  TString St_fhseltrackPtMax[18];
  for(unsigned int i = 0; i < N_ptbin_nseltrack; ++i)
  {
    if( i == (N_ptbin_nseltrack-1) )
    {
      fhseltrackPtMax[i] = 10000;
    }
    else
    {
      fhseltrackPtMax[i] = fhseltrackPtMin[i+1];
    }
    St_fhseltrackPtMin[i] = Form("%4.0f",fhseltrackPtMin[i]);
    St_fhseltrackPtMin[i].ReplaceAll(" ","0");

    St_fhseltrackPtMax[i] = Form("%4.0f",fhseltrackPtMax[i]);
    St_fhseltrackPtMax[i].ReplaceAll(" ","0");
  } 
  TH1D *hData_tracks[18]; //[N_ptbin_nseltrack];
  for( unsigned int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
  {
    TString hName = "hData_tracks_Pt_"+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt];
    hData_tracks[iPt] = new TH1D(hName,hName,51, -0.5, 50.5);
  }

  TH1D* hData_All_nseltracks      = new TH1D("hData_All_nseltracks","nb. of tracks",51,-0.5,50.5); // No reweighting
  TH1D* hData_nseltracks      = new TH1D("hData_nseltracks","nb. of tracks",51,-0.5,50.5);     // PU and pthat reweighting
  TH1D* hData_W_nseltracks      = new TH1D("hData_W_nseltracks","nb. of tracks",51,-0.5,50.5);     // selTracks, PU and pthat reweighting
 //PVH 

///////////////////////////////////////////////////////////////////

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Total Entries : " << nentries << std::endl;
  
  Long64_t nbytes = 0, nb = 0;

  bool TagPos[3] = {false, false, false}, TagNeg[3] = {false, false, false}; //, Veto = false;  
  float varpos, varneg, varPos;
  int numjet = -1, njets = -1;
  int ntagjet = 0;
  int allevents = 0;
  
  int e1 = int(1e+1);
  int e2 = int(1e+2);
  int e3 = int(1e+3);
  int e4 = int(1e+4);
  int e5 = int(1e+5);
  int e6 = int(1e+6);
  int e7 = int(1e+7);
  int e8 = int(1e+8);

//@@
  //int itest = 0;
//@@

// Weight in pthat: QCD Fall11 prod
  float ww0 = 1.;  // pthat weight
  float ww = 1.;   // pthat and PU weight

// Read the weights related to the trigger prescales for data:
  LTANA::PS ps;
  ps.init(1);
// Get the trigger name from IntCut (to be used ig getPS...)
  TString PS_ToUse = "";
  if ( IntCut ==  40 )
  {
    PS_ToUse = "PFJet40";
    weightPU_file += "40.out";
  }else if ( IntCut ==  60 )
  {
    PS_ToUse = "PFJet60";
    weightPU_file += "60.out";
  }else if ( IntCut ==  80 )
  {
    PS_ToUse = "PFJet80";
    weightPU_file += "80.out";
  }else if ( IntCut ==  140 )
  {
    PS_ToUse = "PFJet140";
    weightPU_file += "140.out";
  }else if ( IntCut ==  200 )
  {
    PS_ToUse = "PFJet200";
    weightPU_file += "200.out";
  }else if ( IntCut ==  260 )
  {
    PS_ToUse = "PFJet260";
    weightPU_file += "260.out";
  }else if ( IntCut ==  320 )
  {
    PS_ToUse = "PFJet320";
    weightPU_file += "320.out";
  }else if ( IntCut ==  400 )
  {
    PS_ToUse = "PFJet400";
    weightPU_file += "400.out";
  }else if ( IntCut ==  450 )
  {
    PS_ToUse = "PFJet450";
    weightPU_file += "450.out";
  }else if ( IntCut ==  500 )
  {
    PS_ToUse = "PFJet500";
    weightPU_file += "500.out";
  }else if ( IntCut ==  550 )
  {
    PS_ToUse = "PFJet550";
    weightPU_file += "550.out";
  }else if ( IntCut ==  0 )
  {
    PS_ToUse = "PFJet400";
    weightPU_file += "40.out";
  }else
  {
    cout << "Could not find proper prescale file for IntCut = " << IntCut << endl;
  }

// Read pthat weight from file:
  float pthatMin[200], pthatMax[200], pthatweight[200];
  unsigned int maxLines = 200;
  ReadpthatWeight(pthatMin, pthatMax, pthatweight, maxLines, weightPthat_file);

// Read PU weight from file:
  float WeightPU[MAXPU]; //, WeightPUmin[MAXPU], WeightPUmax[MAXPU];
  for (int k = 0; k < MAXPU; k++) WeightPU[k] = 1.;

  unsigned int maxPU = MAXPU; // Must be coherent (less or equal dim of vector) with definition above !
  //if( weightPU_file != "inJetAna" )
  if( !weightPU_file.Contains("inJetAna") )
  {
    nLinesPUFile = ReadPUWeight(WeightPU, maxPU, weightPU_file); // This is the number of lines red.
    if ( nLinesPUFile == 0 ) return;  // If not a success, stop immediately.
    //if ( ReadPUWeight(WeightPU, maxPU, weightPU_file) == 0 ) return;  // If not a success, stop immediately.
  }
  else
  {
    // MC pu reweighting
    if ( Year == 2011 )
    {
      if ( BTagMu )
      {
// WeightPUmin is computed from 2011A alone
// WeightPUmax is computed from 2011A+B but without runs 165088-167913
        WeightPU[0] = 0.42814;
        WeightPU[1] = 0.852011;
        WeightPU[2] = 1.23782;
        WeightPU[3] = 1.50857;
        WeightPU[4] = 1.67743;
        WeightPU[5] = 1.71335;
        WeightPU[6] = 1.63117;
        WeightPU[7] = 1.4899;
        WeightPU[8] = 1.31817;
        WeightPU[9] = 1.13156;
        WeightPU[10] = 0.959846;
        WeightPU[11] = 0.793424;
        WeightPU[12] = 0.658323;
        WeightPU[13] = 0.536642;
        WeightPU[14] = 0.442476;
        WeightPU[15] = 0.364995;
        WeightPU[16] = 0.301626;
        WeightPU[17] = 0.245887;
        WeightPU[18] = 0.199188;
        WeightPU[19] = 0.158801;
        WeightPU[20] = 0.125015;
        WeightPU[21] = 0.0965054;
        WeightPU[22] = 0.0743175;
        WeightPU[23] = 0.0552261;
        WeightPU[24] = 0.0410448;
        WeightPU[25] = 0.0300012;
        WeightPU[26] = 0.0216022;
        WeightPU[27] = 0.01554;
        WeightPU[28] = 0.0109019;
        WeightPU[29] = 0.0077661;
        WeightPU[30] = 0.00527189;
        WeightPU[31] = 0.00367071;
        WeightPU[32] = 0.00247471;
        WeightPU[33] = 0.0016687;
        WeightPU[34] = 0.00223246;
        //WeightPUmin[0] = 0.508648;
        //WeightPUmin[1] = 1.00906;
        //WeightPUmin[2] = 1.45877;
        //WeightPUmin[3] = 1.76461;
        //WeightPUmin[4] = 1.94019;
        //WeightPUmin[5] = 1.9487;
        //WeightPUmin[6] = 1.80937;
        //WeightPUmin[7] = 1.59264;
        //WeightPUmin[8] = 1.33497;
        //WeightPUmin[9] = 1.06071;
        //WeightPUmin[10] = 0.80788;
        //WeightPUmin[11] = 0.577805;
        //WeightPUmin[12] = 0.397835;
        //WeightPUmin[13] = 0.257752;
        //WeightPUmin[14] = 0.162172;
        //WeightPUmin[15] = 0.0985342;
        //WeightPUmin[16] = 0.0582882;
        //WeightPUmin[17] = 0.0332803;
        //WeightPUmin[18] = 0.0185821;
        //WeightPUmin[19] = 0.0100938;
        //WeightPUmin[20] = 0.0053696;
        //WeightPUmin[21] = 0.00278412;
        //WeightPUmin[22] = 0.00143351;
        //WeightPUmin[23] = 0.000709667;
        //WeightPUmin[24] = 0.000350299;
        //WeightPUmin[25] = 0.00016959;
        //WeightPUmin[26] = 8.06722e-05;
        //WeightPUmin[27] = 3.82429e-05;
        //WeightPUmin[28] = 1.76357e-05;
        //WeightPUmin[29] = 8.23732e-06;
        //WeightPUmin[30] = 3.65711e-06;
        //WeightPUmin[31] = 1.66109e-06;
        //WeightPUmin[32] = 7.28646e-07;
        //WeightPUmin[33] = 3.18864e-07;
        //WeightPUmin[34] = 2.06626e-07;
        //WeightPUmax[0] = 0.321623;
        //WeightPUmax[1] = 0.633854;
        //WeightPUmax[2] = 0.929572;
        //WeightPUmax[3] = 1.1604;
        //WeightPUmax[4] = 1.34017;
        //WeightPUmax[5] = 1.44133;
        //WeightPUmax[6] = 1.46356;
        //WeightPUmax[7] = 1.44183;
        //WeightPUmax[8] = 1.3872;
        //WeightPUmax[9] = 1.30012;
        //WeightPUmax[10] = 1.20261;
        //WeightPUmax[11] = 1.07688;
        //WeightPUmax[12] = 0.957174;
        //WeightPUmax[13] = 0.824407;
        //WeightPUmax[14] = 0.708206;
        //WeightPUmax[15] = 0.601179;
        //WeightPUmax[16] = 0.506317;
        //WeightPUmax[17] = 0.417744;
        //WeightPUmax[18] = 0.340911;
        //WeightPUmax[19] = 0.272991;
        //WeightPUmax[20] = 0.215471;
        //WeightPUmax[21] = 0.166586;
        //WeightPUmax[22] = 0.128399;
        //WeightPUmax[23] = 0.0954633;
        //WeightPUmax[24] = 0.0709704;
        //WeightPUmax[25] = 0.0518837;
        //WeightPUmax[26] = 0.0373623;
        //WeightPUmax[27] = 0.0268788;
        //WeightPUmax[28] = 0.0188572;
        //WeightPUmax[29] = 0.0134334;
        //WeightPUmax[30] = 0.00911912;
        //WeightPUmax[31] = 0.00634949;
        //WeightPUmax[32] = 0.0042807;
        //WeightPUmax[33] = 0.00288649;
        //WeightPUmax[34] = 0.00107161;
        nLinesPUFile = 34+1; // Must be consistent with max defined weight +1
      }
      else if ( IntCut == 300 )
      {
// WeightPUmin is computed from 2011A+B but without runs 178098-180252
// WeightPUmax is computed from 2011A+B but without runs 165088-167913
        WeightPU[0] = 0.0236327;
        WeightPU[1] = 0.186895;
        WeightPU[2] = 0.411244;
        WeightPU[3] = 0.715173;
        WeightPU[4] = 0.976423;
        WeightPU[5] = 1.17936;
        WeightPU[6] = 1.29761;
        WeightPU[7] = 1.36739;
        WeightPU[8] = 1.39617;
        WeightPU[9] = 1.40456;
        WeightPU[10] = 1.41265;
        WeightPU[11] = 1.41307;
        WeightPU[12] = 1.45094;
        WeightPU[13] = 1.46775;
        WeightPU[14] = 1.50175;
        WeightPU[15] = 1.52329;
        WeightPU[16] = 1.5588;
        WeightPU[17] = 1.58469;
        WeightPU[18] = 1.59661;
        WeightPU[19] = 1.57766;
        WeightPU[20] = 1.56339;
        WeightPU[21] = 1.52986;
        WeightPU[22] = 1.50984;
        WeightPU[23] = 1.41909;
        WeightPU[24] = 1.38554;
        WeightPU[25] = 1.34946;
        WeightPU[26] = 1.28489;
        WeightPU[27] = 1.15493;
        WeightPU[28] = 1.03597;
        WeightPU[29] = 0.984955;
        WeightPU[30] = 0.876273;
        WeightPU[31] = 0.772276;
        WeightPU[32] = 0.720553;
        WeightPU[33] = 0.707666;
        WeightPU[34] = 0.378344;
        //WeightPUmin[0] = 0.293038;
        //WeightPUmin[1] = 0.58781;
        //WeightPUmin[2] = 0.86736;
        //WeightPUmin[3] = 1.08442;
        //WeightPUmin[4] = 1.2539;
        //WeightPUmin[5] = 1.35408;
        //WeightPUmin[6] = 1.38846;
        //WeightPUmin[7] = 1.39185;
        //WeightPUmin[8] = 1.37383;
        //WeightPUmin[9] = 1.33004;
        //WeightPUmin[10] = 1.27531;
        //WeightPUmin[11] = 1.18245;
        //WeightPUmin[12] = 1.08192;
        //WeightPUmin[13] = 0.950148;
        //WeightPUmin[14] = 0.822697;
        //WeightPUmin[15] = 0.695634;
        //WeightPUmin[16] = 0.577288;
        //WeightPUmin[17] = 0.465037;
        //WeightPUmin[18] = 0.367816;
        //WeightPUmin[19] = 0.283842;
        //WeightPUmin[20] = 0.214971;
        //WeightPUmin[21] = 0.158955;
        //WeightPUmin[22] = 0.116885;
        //WeightPUmin[23] = 0.0827463;
        //WeightPUmin[24] = 0.0584831;
        //WeightPUmin[25] = 0.0405946;
        //WeightPUmin[26] = 0.0277256;
        //WeightPUmin[27] = 0.0188999;
        //WeightPUmin[28] = 0.0125536;
        //WeightPUmin[29] = 0.00846052;
        //WeightPUmin[30] = 0.00543007;
        //WeightPUmin[31] = 0.00357258;
        //WeightPUmin[32] = 0.00227472;
        //WeightPUmin[33] = 0.00144797;
        //WeightPUmin[34] = 0.000481535;
        //WeightPUmax[0] = 0.138424;
        //WeightPUmax[1] = 0.28206;
        //WeightPUmax[2] = 0.435051;
        //WeightPUmax[3] = 0.582968;
        //WeightPUmax[4] = 0.739521;
        //WeightPUmax[5] = 0.892853;
        //WeightPUmax[6] = 1.03605;
        //WeightPUmax[7] = 1.18049;
        //WeightPUmax[8] = 1.3209;
        //WeightPUmax[9] = 1.43841;
        //WeightPUmax[10] = 1.53504;
        //WeightPUmax[11] = 1.56652;
        //WeightPUmax[12] = 1.56158;
        //WeightPUmax[13] = 1.48143;
        //WeightPUmax[14] = 1.37646;
        //WeightPUmax[15] = 1.24272;
        //WeightPUmax[16] = 1.09721;
        //WeightPUmax[17] = 0.937959;
        //WeightPUmax[18] = 0.785909;
        //WeightPUmax[19] = 0.641764;
        //WeightPUmax[20] = 0.513977;
        //WeightPUmax[21] = 0.401756;
        //WeightPUmax[22] = 0.312279;
        //WeightPUmax[23] = 0.233716;
        //WeightPUmax[24] = 0.174679;
        //WeightPUmax[25] = 0.128265;
        //WeightPUmax[26] = 0.0927113;
        //WeightPUmax[27] = 0.0669142;
        //WeightPUmax[28] = 0.0470794;
        //WeightPUmax[29] = 0.0336246;
        //WeightPUmax[30] = 0.0228792;
        //WeightPUmax[31] = 0.0159644;
        //WeightPUmax[32] = 0.0107841;
        //WeightPUmax[33] = 0.00728486;
        //WeightPUmax[34] = 0.00271264;
        nLinesPUFile = 34+1; // Must be consistent with max defined weight +1
      }
      else
      {
        WeightPU[0] = 0.0428918;
        WeightPU[1] = 0.335376;
        WeightPU[2] = 0.721533;
        WeightPU[3] = 1.20833;
        WeightPU[4] = 1.55725;
        WeightPU[5] = 1.7342;
        WeightPU[6] = 1.71538;
        WeightPU[7] = 1.58616;
        WeightPU[8] = 1.39386;
        WeightPU[9] = 1.19426;
        WeightPU[10] = 1.02425;
        WeightPU[11] = 0.88516;
        WeightPU[12] = 0.802558;
        WeightPU[13] = 0.735522;
        WeightPU[14] = 0.698884;
        WeightPU[15] = 0.67222;
        WeightPU[16] = 0.662791;
        WeightPU[17] = 0.656711;
        WeightPU[18] = 0.649984;
        WeightPU[19] = 0.634296;
        WeightPU[20] = 0.622947;
        WeightPU[21] = 0.605549;
        WeightPU[22] = 0.594589;
        WeightPU[23] = 0.556591;
        WeightPU[24] = 0.541617;
        WeightPU[25] = 0.526013;
        WeightPU[26] = 0.499605;
        WeightPU[27] = 0.448078;
        WeightPU[28] = 0.40112;
        WeightPU[29] = 0.380669;
        WeightPU[30] = 0.338094;
        WeightPU[31] = 0.297503;
        WeightPU[32] = 0.277176;
        WeightPU[33] = 0.271852;
        WeightPU[34] = 0.145013;
        nLinesPUFile = 34+1; // Must be consistent with max defined weight +1
      }
    }
  }

// Read NSelTrack weight from files:

  unsigned int N_SelTracks = 52;
  float WeightNS[18*52];
  ReadNSelTracksWeight(WeightNS, N_SelTracks, N_ptbin_nseltrack, weightNS_file, St_fhseltrackPtMin, St_fhseltrackPtMax);

// Read JSon file:

  unsigned int N_MaxRun = 100;
  unsigned int N_MaxLumiRangePerRun = 1000;
  unsigned int RunSize;
  unsigned int JS_Run[100];
  unsigned int JS_Run_NLumi[100];
  unsigned int JS_LumiMin[100*1000];
  unsigned int JS_LumiMax[100*1000];

  int JSON_RunN = ReadJSon(N_MaxRun, N_MaxLumiRangePerRun, RunSize, JS_Run, JS_Run_NLumi, JS_LumiMin, JS_LumiMax, JSONFile);
  if ( JSON_RunN < 0 )
  {
    cout << "Failed to read JSon file with return value " << JSON_RunN << endl;
    return;
  } 
  
// 

////////////////
// Event loop //
////////////////

  bool testJSONSelect = true;
  Long64_t NTreeEvents = 0;

  int N_repEvts = 0;
  for (Long64_t jentry=0; jentry<nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    // Count events in the chain
    NTreeEvents++;
//----------------------------------------------------------------
//Start selection (remove bad runs, bad PtVspthat...
// Bad run in promptreco data (temporary before reprocessing)



    if ( Run > 1 )	// Selections for data files
    {
      // Select using JSon file
      if ( JSON_RunN > 0 )
      {
        bool notInRange = true;
        for( unsigned int i = 0; i <= RunSize && notInRange ; ++i)
        {
          for( unsigned int j = 0; j < JS_Run_NLumi[i]; ++j )
          {
            if( testJSONSelect )
            {
              cout << Run << " : " << LumiBlock << " =?= " << JS_Run[i] << " : " << JS_LumiMin[i*1000+j] << " success = " << ( Run == (int) JS_Run[i] && ( LumiBlock >= (int) JS_LumiMin[i*1000+j] && LumiBlock <= (int) JS_LumiMax[i*1000+j] ) ) << endl;
              testJSONSelect = false;
            }
            if( Run == (int) JS_Run[i] && ( LumiBlock >= (int) JS_LumiMin[i*1000+j] && LumiBlock <= (int) JS_LumiMax[i*1000+j] ) )
            {
              notInRange = false;
              break;
            }
          }
        }
        if( notInRange ) continue; 		// Go to next event
      }

      if( TrigType == "2016_BX_4X2" )
      {
        if( !( Run < 273440 || (Run > 274000 && Run < 274100) ) ) continue;
      }
      else if ( TrigType == "2016_BX_4X3" )
      {
        if( !( Run >= 273440 && Run <= 273550 ) ) continue;
      }
      else if ( TrigType == "2016_BX_4X4" )
      {
        if( !( (Run > 273550 && Run <= 274000) || (Run > 274100 && Run < 274150 ) ) ) continue;
      }
      else if ( TrigType == "2016_BX_4X5" )
      {
        if( !( Run > 274150 && Run <= 274170 ) ) continue;
      }
      else if ( TrigType == "2016_BX_4X6" )
      {
        if( !( Run > 274170 && Run <= 274300 ) ) continue;
      }
      else if ( TrigType == "2016_BX_Cont" )
      {
        if( !( Run >= 274300 ) ) continue;
      }
      else if ( TrigType == "2016_flat1" )
      {
        if( !( Run >= 274300 && Run <= 276097 ) ) continue;
      }
      else if ( TrigType == "2016_flat2" )
      {
        if( !( Run >= 275784 && Run <= 276097 ) ) continue;
      }
      else if ( TrigType == "2016_76fb_L1T" )
      {
        if( !( Run >= 273158 && Run <= 276097 ) ) continue;
      }
      else if ( TrigType == "Run2017Elt304671-17Nov2017-v1" )	// Runs from E period, below 304671
      {
        if( Run >= 304671 ) continue;
      }
      else if ( TrigType == "Run2017Ele304671-17Nov2017-v1" )	// Runs from E period, below or eq 304671
      {
        if( Run > 304671 ) continue;
      }
      else if ( TrigType == "Run2017Eeq304671-17Nov2017-v1" ) // Runs from E period, above 304671
      {
        if( Run != 304671 ) continue;
      }
      else if ( TrigType == "Run2017EGT304671-17Nov2017-v1" ) // Runs from E period, above 304671
      {
        if( Run <= 304671 ) continue;
      }
      // Filter laser events which show up with large nJet or large Jet_pt[0]
      if( nJet > 0 )
      {
        if( nJet > 16 || Jet_pt[0] > 1600 ) continue;
      }
    } 
    else if( Run < 0 )
    {
      bool Skip_StrangeEvts = false;
      if( Skip_StrangeEvts )
      {
        if( pthat <  30. ) { if( Jet_pt[0] > 130. ) continue; }		// Go to next event
        else if( pthat <  50. ) { if( Jet_pt[0] > 130. ) continue; }	// Go to next event
        else if( pthat <  80. ) { if( Jet_pt[0] > 180. ) continue; }	// Go to next event
        else if( pthat <  120. ) { if( Jet_pt[0] > 250. ) continue; }	// Go to next event
        else if( pthat <  170. ) { if( Jet_pt[0] > 310. ) continue; }	// Go to next event
        else if( pthat <  300. ) { if( Jet_pt[0] > 460. ) continue; }	// Go to next event

        // Add special trick to remove some kind of repetitive event... Or the result of a bug:
        //   some events with more than 3 jets with values around 382, 267 and 135
        //                    and corresponding eta values around -0.6, -2.12, -1.27
        //                    should be skipped
        if( nJet > 0 && pthat < 200 )
        {
          if ( Jet_pt[0] > 1000 )
          {
            N_repEvts++;
            cout << "Warning: one more strange event type: " << N_repEvts << endl;
            cout << "Event = " << Evt << " \tnJet = " << nJet << endl;
            cout << "Jet0 with pt = " << Jet_pt[0] << " \t and eta = " << Jet_eta[0] << endl;
            continue;		// Go to next event
          }
        }
        if( nJet >= 3 && pthat < 200 )
        {
          bool BadJet0 = false;
          bool BadJet1 = false;
          bool BadJet2 = false;
          Double_t BadJet0_pt, BadJet0_eta;
          Double_t BadJet1_pt, BadJet1_eta;
          Double_t BadJet2_pt, BadJet2_eta;
          for( unsigned int iJet = 0; (int) iJet < nJet; ++iJet )
          {
            if( Jet_pt[iJet] > 360 && Jet_pt[iJet] < 400 && Jet_eta[iJet] > -0.70 && Jet_eta[iJet] < -0.60 )
            {
              BadJet0_pt = Jet_pt[iJet];
              BadJet0_eta = Jet_eta[iJet];
              BadJet0 = true;
              cout << "Warning0" << endl;
            }
            if( Jet_pt[iJet] > 250 && Jet_pt[iJet] < 305 && Jet_eta[iJet] > -2.20 && Jet_eta[iJet] < -2.10 )
            {
              BadJet1_pt = Jet_pt[iJet];
              BadJet1_eta = Jet_eta[iJet];
              BadJet1 = true;
              cout << "Warning1" << endl;
            }
            if( Jet_pt[iJet] > 125 && Jet_pt[iJet] < 145 && Jet_eta[iJet] > -1.40 && Jet_eta[iJet] < -1.20 )
            {
              BadJet2_pt = Jet_pt[iJet];
              BadJet2_eta = Jet_eta[iJet];
              BadJet2 = true;
              cout << "Warning2" << endl;
            }
          }
          if( BadJet0 && BadJet1 && BadJet2 )
          {
            N_repEvts++;
            cout << "Warning: one more repetitive event type: " << N_repEvts << endl;
            cout << "Event = " << Evt << " \tnJet = " << nJet << endl;
            cout << "Jet0 with pt = " << BadJet0_pt << " \t and eta = " << BadJet0_eta << endl;
            cout << "Jet1 with pt = " << BadJet1_pt << " \t and eta = " << BadJet1_eta << endl;
            cout << "Jet2 with pt = " << BadJet2_pt << " \t and eta = " << BadJet2_eta << endl;
            continue;		// Go to next event
          }
        }
      }
    }
//
//



    if( nPV < minPV || nPV > maxPV) continue; // Run only on some nPV events
    if( nJet < 1 ) continue;
// Filter event if maxpt < given cut || > given cut
    if( Jet_pt[0] < minCutJetPtMax || Jet_pt[0] > maxCutJetPtMax ) continue; // Go to next event

//End of bad runs, bad PtVspthat selection 
//----------------------------------------------------------------

// pthat reweighting procedure
//$$
    int npv = nPV, npu = nPU; 
    if ( truePU ) npu = nPUtrue; 
    if ( npv >= MAXPU ) npv = MAXPU;
    if ( npu >= (int) MAXPU-1 ) npu = MAXPU-1;

    // Filter against monsters in MC
    //$$
    bool MONSTER = false;
    if ( pthat > 0 && pthat < 120. ) {
      for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++) {
        if ( Jet_pt[ijet] > pthat * 7. ) MONSTER = true;
      }
    }
    if ( MONSTER ) continue;


    ww = ww0 = 0.;
    if ( Run >= 0 ) 
    {
      ww0 = 1;		// No pthat reweighting for Data
      //ww = 1;		// Default value
      //if ( IntCut ==  40 ) ww = ps.getPSWeight((std::string) PS_ToUse,(int) Run,(int) LumiBlock);    // No reweighting for Data
// Temporarily comment the prescale reweight for 2018 and set it to 1 instead of getPSWeight(... //
      ww = 1; // See line above. ps.getPSWeight((std::string) PS_ToUse,(int) Run,(int) LumiBlock);    // No reweighting for Data
    } else {
      bool inPthatBin = false;
      for( unsigned int i = 0; i < maxLines; ++i )
      {
        if( pthat > pthatMin[i] && pthat < pthatMax[i] )
        {
          ww0 = pthatweight[i];
          inPthatBin = true;
          break;
        }
      }
      if( !inPthatBin ) continue;
      if( weightPU_file.Contains("_nPV") )      // Then weight according to nPV
      {
        //cout << "Here is Weightfilename: " << weightPU_file << endl;
        //cout << "Yes it does contain nPV" << endl;
        if( nPV > 0 && nPV < nLinesPUFile ) ww = ww0 * WeightPU[nPV];
        else if( nPV >= nLinesPUFile ) ww =  ww0 * WeightPU[nLinesPUFile-1]; // If weight not defined, use last defined one
        else ww = 0;
      } else {					// Weight according to npu
        //cout << "Here is Weightfilename: " << weightPU_file << endl;
        //cout << "NO it doesn't contain nPV" << endl;
        ww = ww0 * WeightPU[npu]; 
      }
    }
//$$
    hAllFlav_pthatAll->Fill( pthat, 1. );
    hAllFlav_pthat->Fill( pthat , ww0 );

 
//----------------------------------------------------------------
// Trigger selection
//$$

 if( nJet > 0 ) hBeforeTrig_MaxJetPt->Fill(Jet_pt[0],ww0);

 bool Jet30  = false, Jet60  = false, Jet150 = false, Jet190 = false, Jet240 = false;
 bool Jet40  = false, Jet80  = false, Jet140 = false;
 bool Jet200 = false, Jet260 = false, Jet320 = false, Jet400 = false, Jet450 = false, Jet500 = false;
 bool Jet20  = false, Jet70  = false, Jet110 = false, Jet300 = false;

 int njet20 = 0, njet40 = 0, njet60 = 0, njet70 = 0, njet110 = 0, njet300 = 0;
 int njet80 = 0, njet140 = 0, njet200 = 0, njet260 = 0, njet320 = 0;
 int njet400 = 0, njet450 = 0, njet500 = 0;

 int triggerIdx = 0, bitIdx = 0;

 if ( !BTagMu ) {
   if ( Year == 2011 ) {
     triggerIdx = 1;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet30  = true;
       if( nJet > 0 ) hAfterTrig30_MaxJetPt->Fill(Jet_pt[0],ww0);
     }

     triggerIdx = 4;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet60  = true;
       if( nJet > 0 ) hAfterTrig60_MaxJetPt->Fill(Jet_pt[0],ww0);
     }

     triggerIdx = 6;
     bitIdx = int(triggerIdx/32);
     if( Run > 0 )
     {
       if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
       {
         Jet80  = true;
         if( nJet > 0 ) hAfterTrig80_MaxJetPt->Fill(Jet_pt[0],ww0);
       }
     } else {
       if( nJet > 0 && Jet_pt[0] >= 90 ) // Used to simulate the inexisting MC trigger. Also used for Data for compatibility
       {
         Jet80  = true;
         if( nJet > 0 ) hAfterTrig80_MaxJetPt->Fill(Jet_pt[0],ww0);
       }
     }

     triggerIdx =  9;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet110 = true;
       if( nJet > 0 ) hAfterTrig110_MaxJetPt->Fill(Jet_pt[0],ww0);
     }

     triggerIdx = 11;
     bitIdx = int(triggerIdx/32);
     if( Run > 0 )
     {
       if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
       {
         Jet150 = true;
         if( nJet > 0 ) hAfterTrig150_MaxJetPt->Fill(Jet_pt[0],ww0);
       }
     } else {
       if( nJet > 0 && Jet_pt[0] >= 70 ) // Used to simulate the inexisting MC trigger. Also used for Data for compatibility
       {
         Jet150 = true;
         if( nJet > 0 ) hAfterTrig150_MaxJetPt->Fill(Jet_pt[0],ww0);
       }
     }

     triggerIdx = 14;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet190 = true;
       if( nJet > 0 ) hAfterTrig190_MaxJetPt->Fill(Jet_pt[0],ww0);
     }

     triggerIdx = 16;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet240 = true;
       if( nJet > 0 ) hAfterTrig240_MaxJetPt->Fill(Jet_pt[0],ww0);
     }

     triggerIdx = 18;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
     {
       Jet300 = true;
       if( nJet > 0 ) hAfterTrig300_MaxJetPt->Fill(Jet_pt[0],ww0);
     }
   }
   if ( Year == 2012 ) {
     triggerIdx = 2;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;

     triggerIdx = 7;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;

     triggerIdx = 12;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;

     triggerIdx = 15;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;
 
     triggerIdx = 17;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet260 = true;

     triggerIdx = 19;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet320 = true;

     triggerIdx = 21;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet400 = true;

     triggerIdx = 87;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60 = true;

     triggerIdx = 88;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

     triggerIdx = 89;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet500 = true;

     for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++) {
       float residual = 1.;
       float ptjet = Jet_pt[ijet] * residual;
       //float etajet = fabs(Jet_eta[ijet]);

       //hData_Trigger_JetPt->Fill( ptjet , ww );
       //if ( Jet40 )  hData_Trig40_JetPt->Fill(  ptjet , ww );
       //if ( Jet80 )  hData_Trig80_JetPt->Fill(  ptjet , ww );
       //if ( Jet140 ) hData_Trig140_JetPt->Fill( ptjet , ww );
       //if ( Jet200 ) hData_Trig200_JetPt->Fill( ptjet , ww );
       //if ( Jet260 ) hData_Trig260_JetPt->Fill( ptjet , ww );
       //if ( Jet320 ) hData_Trig320_JetPt->Fill( ptjet , ww );

       if ( ptjet >  50. ) njet40++;
       if ( ptjet >  70. ) njet60++;
       if ( ptjet > 100. ) njet80++;
       if ( ptjet > 160. ) njet140++;
       if ( ptjet > 220. ) njet200++;
       if ( ptjet > 300. ) njet260++;
       if ( ptjet > 360. ) njet320++;
       if ( ptjet > 420. ) njet400++;
       if ( ptjet > 470. ) njet450++;
       if ( ptjet > 520. ) njet500++;
     }
     if ( njet40  < 1 ) Jet40  = false;
     if ( njet60  < 1 ) Jet60  = false;
     if ( njet80  < 1 ) Jet80  = false;
     if ( njet140 < 1 ) Jet140 = false;
     if ( njet200 < 1 ) Jet200 = false;
     if ( njet260 < 1 ) Jet260 = false;
     if ( njet320 < 1 ) Jet320 = false;
     if ( njet400 < 1 ) Jet400 = false;
     if ( njet450 < 1 ) Jet450 = false;
     if ( njet500 < 1 ) Jet500 = false;
   }
   if ( Year == 2016 ) {
     triggerIdx = 0;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;

     triggerIdx = 1;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60  = true;

     triggerIdx = 2;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;

     triggerIdx = 3;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;

     triggerIdx = 4;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;

     triggerIdx = 5;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet260 = true;

     triggerIdx = 6;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet320 = true;

     triggerIdx = 7;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet400 = true;

     triggerIdx = 8;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

     triggerIdx = 9;
     bitIdx = int(triggerIdx/32);
     if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet500 = true;

     for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++) {
       float residual = 1.;
       float ptjet = Jet_pt[ijet] * residual;
       //float etajet = fabs(Jet_eta[ijet]);

       if ( ptjet >  50. ) njet40++;
       if ( ptjet >  70. ) njet60++;
       if ( ptjet > 100. ) njet80++;
       if ( ptjet > 160. ) njet140++;
       if ( ptjet > 220. ) njet200++;
       if ( ptjet > 300. ) njet260++;
       if ( ptjet > 360. ) njet320++;
       if ( ptjet > 420. ) njet400++;
       if ( ptjet > 470. ) njet450++;
       if ( ptjet > 520. ) njet500++;
     }
     if ( njet40  < 1 ) Jet40  = false;
     if ( njet60  < 1 ) Jet60  = false;
     if ( njet80  < 1 ) Jet80  = false;
     if ( njet140 < 1 ) Jet140 = false;
     if ( njet200 < 1 ) Jet200 = false;
     if ( njet260 < 1 ) Jet260 = false;
     if ( njet320 < 1 ) Jet320 = false;
     if ( njet400 < 1 ) Jet400 = false;
     if ( njet450 < 1 ) Jet450 = false;
     if ( njet500 < 1 ) Jet500 = false;
   }
 }

 else if ( BTagMu ) {

   triggerIdx = 35;
   bitIdx = int(triggerIdx/32);
   if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet20  = true;

   triggerIdx = 41;
   bitIdx = int(triggerIdx/32);
   if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;

   triggerIdx = 44;
   bitIdx = int(triggerIdx/32);
   if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet70  = true;

   triggerIdx = 47;
   bitIdx = int(triggerIdx/32);
   if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet110 = true;

   triggerIdx = 50;
   bitIdx = int(triggerIdx/32);
   if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet300 = true;

   for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++) {
     float residual = 1.;
     float ptjet = Jet_pt[ijet] * residual;
     //float etajet = fabs(Jet_eta[ijet]);

     //hData_Trigger_JetPt->Fill( ptjet , ww );
     //if ( Jet20 )  hData_Trig20_JetPt->Fill(  ptjet , ww );
     //if ( Jet40 )  hData_Trig40_JetPt->Fill(  ptjet , ww );
     //if ( Jet70 )  hData_Trig70_JetPt->Fill(  ptjet , ww );
     //if ( Jet110 ) hData_Trig110_JetPt->Fill( ptjet , ww );
     //if ( Jet300 ) hData_Trig300_JetPt->Fill( ptjet , ww );

     if ( ptjet > 20. )  njet20++;
     if ( ptjet > 50. )  njet40++;
     if ( ptjet > 80. )  njet70++;
     if ( ptjet > 120. ) njet110++;
     if ( ptjet > 320. ) njet300++;
   }
   if ( Year == 2012 ) {   
     if ( njet20  < 1 ) Jet20  = false;
     if ( njet40  < 1 ) Jet40  = false;
     if ( njet70  < 1 ) Jet70  = false;
     if ( njet110 < 1 ) Jet110 = false;
     if ( njet300 < 1 ) Jet300 = false;
   }
 }
      
// Trigger selection
 if ( !BTagMu && !TopAna ) {
   if ( IntCut ==  30 && !Jet30 )  continue;
   if ( IntCut ==  40 && !Jet40 )  continue;
   if ( IntCut ==  60 && !Jet60 )  continue;
   if ( IntCut ==  80 && !Jet80 )  continue;
   //if ( IntCut == 80 && Jet_pt[0] < minCutJetPtMax ) continue;		// Go to next event
   if ( IntCut == 110 && !Jet110 ) continue;
   if ( IntCut == 140 && !Jet140 ) continue;
   if ( IntCut == 150 && !Jet150 ) continue;
   //if ( IntCut == 150 && Jet_pt[0] < minCutJetPtMax ) continue;		// Go to next event
   if ( IntCut == 190 && !Jet190 ) continue;
   if ( IntCut == 200 && !Jet200 ) continue;
   if ( IntCut == 240 && !Jet240 ) continue;
   if ( IntCut == 260 && !Jet260 ) continue;
   if ( IntCut == 300 && !Jet300 ) continue;
   if ( IntCut == 320 && !Jet320 ) continue;
   if ( IntCut == 400 && !Jet400 ) continue;
   if ( IntCut == 450 && !Jet450 ) continue;
   if ( IntCut == 500 && !Jet500 ) continue;
   if ( IntCut == 1 && Year == 2011
        && !Jet30 && !Jet60 && !Jet80 && !Jet110 && !Jet150 
                  && !Jet190 && !Jet240 && !Jet300 ) continue;
   if ( IntCut == 1 && Year == 2012
        && !Jet40 && !Jet60 && !Jet80 && !Jet140
                  && !Jet200 && !Jet260 && !Jet320 && !Jet400 && !Jet450 && !Jet500 ) continue;
 }
 else if ( BTagMu ) {
   if ( IntCut == 20  && !Jet20 )  continue;
   if ( IntCut == 40  && !Jet40 )  continue;
   if ( IntCut == 70  && !Jet70 )  continue;
   if ( IntCut == 110 && !Jet110 ) continue;
   if ( IntCut == 300 && !Jet300 ) continue;
   if ( IntCut == 1 && Year == 2011
        && !Jet20 && !Jet40 && !Jet70 && !Jet110 ) continue;
   if ( IntCut == 1 && Year == 2012
        && !Jet20 && !Jet40 && !Jet70 && !Jet110 && !Jet300 ) continue;
 }
      



// End of trigger selection
//----------------------------------------------------------------
//

//----------------------------------------------------------------
// Eta and ptmax selection
    int njet30 = 0;
    float maxJetPt;
    maxJetPt = 0;
    for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++)
    {
      float residual = 1.;
      if ( Run > 0 ) residual = Jet_residual[ijet];
      residual = 1.; // Comment if above has to be applied
      float ptjet = Jet_pt[ijet] * residual;
      //float etajet = fabs(Jet_eta[ijet]);
      // Only used info on jet between EtaMin and EtaMax to make the selection... A bit strange ?
      if ( ptjet > 30. ) njet30++;

      if( ptjet > maxJetPt ) maxJetPt = ptjet;
    }
// End of eta and ptmax selection
//----------------------------------------------------------------


    hAllFlav_pthatTrig->Fill( pthat, 1. );
      
    njets = nJet;

    hData_All_nPV->Fill( npv , ww0 );		// Only pthat reweighted. Weights are "1" if data
    for( unsigned int iPt = 0; iPt < NhPt; ++iPt )
    {
      if( maxJetPt >= fhPtMin[iPt] && maxJetPt < fhPtMax[iPt] )
      {
        hData_JetPtBinned_nPV[iPt]->Fill( npv , ww );
      }
    }
    hData_nPV->Fill( npv , ww );		// pthat and PU reweighted or prescale if data
    hAllFlav_All_nPU->Fill( npu , ww0 );
    hAllFlav_nPU->Fill( nPU , ww );

    if ( npv >= MAXPU ) npv = MAXPU;
    if ( npu >= MAXPU ) npu = MAXPU;

//------------------------------------------------
// Loop on jets
    float tmp_ww = ww;
    //size_t begin=0;		// From Mauro https://gist.github.com/mverzett/886b35bc03303638c6c4795541a67116
    for (unsigned int ijet = 0; ijet < (unsigned int) nJet; ijet++) {

      if ( ijet == 0 ) {
        allevents++;
        if ( allevents%10000 == 0 ) std::cout << "events : " << allevents << std::endl ;

        hData_All_NJets->Fill( njets , ww );
        hData_NJets->Fill( numjet , ww );
        hData_NJet30->Fill( njet30 , ww );
   
        numjet = 0;
        ntagjet = 0;
      }

//    int ntracks = Jet_nLastTrack[ijet] - Jet_nFirstTrack[ijet];
      int ntracks = Jet_ntracks[ijet];
      int nseltracks = Jet_nseltracks[ijet];
      float residual = 1.;
      if ( Run > 0 ) residual = Jet_residual[ijet];
      residual = 1.; // Comment if above has to be applied
      float ptjet = Jet_pt[ijet] * residual;
      float etajet = fabs(Jet_eta[ijet]);
      //float phijet = Jet_phi[ijet];
      int flavour = abs( Jet_flavour[ijet] );
      if( Run > 0 )
      {
        flavour = -1;
      }
      else
      {      
        if ( (flavour != 4) && (flavour != 5)  && (flavour != 21) ) flavour = 1; // This is used if Gluon info is not stored
        // the flavour inforation is not always complete. Anyway, c and b flavour are stored as 4 and 5.
        // if gluon flavor is kept, it is number 21
        // light flavour should contain all the rest (u,d, s)  so the line above should always be correct
        
      }


  // float WeightNS[18*52];
  // float WeightNS[N_ptbin_nseltrack*N_SelTracks];
      if( WeightTracks )
      {
        for( unsigned int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
        {
          if( ptjet >= fhseltrackPtMin[iPt] && ptjet < fhseltrackPtMax[iPt] )
          {
            ww = tmp_ww * WeightNS[iPt*N_SelTracks+nseltracks];
            break;
          }
        }
      }
      hData_All_nseltracks->Fill(nseltracks, 1.);   // No reweighting
      hData_nseltracks->Fill(nseltracks, tmp_ww);   // pthat and PU reweighting
      hData_W_nseltracks->Fill(nseltracks, ww); // pthat, PU and nseltracks reweighting

//  std::cout << " jet " << ijet << " pt eta " << ptjet << " " << etajet 
//       << " nJet numjet " << njets << " " << numjet << std::endl;

//*********************************
// Jet selection

//$$
      if ( Run < 0 && Jet_genpt[ijet] < 8 ) continue;			// Go to next jet
      if ( !(ptjet >= PtMin && ptjet < PtMax) ) continue;		// Go to next jet
//$$
      if ( ptjet > maxPtBin) ptjet = maxPtBin -0.01;
//$$

      numjet++;

      //if ( trig1 == 1 || trig1 == 3 || trig1 == 5 || trig1 >= 7 ) hData_Trig1_JetPt->Fill( ptjet , ww );
      //if ( trig1 == 2 || trig1 == 3 || trig1 >= 6 ) hData_Trig2_JetPt->Fill( ptjet , ww );
      //if ( trig1 >= 4 ) hData_Trig3_JetPt->Fill( ptjet , ww );
      //if ( trig2 == 1 || trig2 == 3 || trig2 == 5 || trig2 >= 7 ) hData_Trig4_JetPt->Fill( ptjet , ww );
      //if ( trig2 == 2 || trig2 == 3 || trig2 >= 6 ) hData_Trig5_JetPt->Fill( ptjet , ww );
      //if ( trig2 >= 4 ) hData_Trig6_JetPt->Fill( ptjet , ww );

//*********************************

      for( unsigned int i = 0; i < 3; ++i )
      {
        TagPos[i] = false;
        TagNeg[i] = false;
      }
      //Veto = false;

      varpos = -1000.;
      varneg = -1000.;
      varPos = -1000.;

// Track history
      int cat = 0, catP = 0, catN = 0;
      int cat1 = 0, cat1P = 0, cat1N = 0;
      int cat2 = 0, cat2P = 0, cat2N = 0;

      bool GamP = false, GamN = false, GamPN = false;
      bool K0sP = false, K0sN = false, K0sPN = false;
      bool LamP = false, LamN = false, LamPN = false;
      bool FakP = false, FakN = false, FakPN = false;
      bool OthP = false, OthN = false, OthPN = false;
   
//***************************
// Tagger

// Combined CSVv2
     if ( TagName == "CSVv2" ) {
        if ( Jet_CombIVF[ijet] > 0. )   varPos = 25.*Jet_CombIVF[ijet];
        if ( Jet_CombIVF_N[ijet] > 0. ) varneg = 25.*Jet_CombIVF_N[ijet];
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_CombIVF[ijet] > TagCut[i] )   TagPos[i] = true;
          if ( Jet_CombIVF_N[ijet] > TagCut[i] ) TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
      }
// DeepCSV
      else if ( TagName == "DeepCSV" ) {
// First test of DeepCSV Use Branch "DeepCSVBDisc" as proposed in Ivan Marchesini mail from 5/12/2016
        if ( Jet_DeepCSVBDisc[ijet] > 0. )   varPos = 25.*Jet_DeepCSVBDisc[ijet];
        //if ( Jet_DeepCSVBDiscP[ijet] > 0. ) varpos = 25.*Jet_DeepCSVBDiscP[ijet];
        if ( Jet_DeepCSVBDiscN[ijet] > 0. ) varneg = 25.*Jet_DeepCSVBDiscN[ijet];
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepCSVBDisc[ijet] > TagCut[i] )   TagPos[i] = true;
          if ( Jet_DeepCSVBDiscN[ijet] > TagCut[i] )  TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
//
      }
      else if ( TagName == "DeepFlavour" ) {
// First test of DeepFlavour Use Branch "DeepFlavourBDisc" as proposed in Ivan Marchesini mail from 5/12/2016
        if ( Jet_DeepFlavourBDisc[ijet] > 0. )   varPos = 25.*Jet_DeepFlavourBDisc[ijet];
        //if ( Jet_DeepFlavourBDiscP[ijet] > 0. ) varpos = 25.*Jet_DeepFlavourBDiscP[ijet];
        if ( Jet_DeepFlavourBDiscN[ijet] > 0. ) varneg = 25.*Jet_DeepFlavourBDiscN[ijet];
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepFlavourBDisc[ijet] > TagCut[i] )   TagPos[i] = true;
          if ( Jet_DeepFlavourBDiscN[ijet] > TagCut[i] )  TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
//
      }
// Combined cMVAv2
/*
      else if ( TagName == "cMVAv2" ) {
        if ( Jet_cMVAv2[ijet] > -1. )   varPos = 12.5*(1+Jet_cMVAv2[ijet]);		// Hopefully this is correct rescaling and test ?!
        //if ( Jet_cMVAv2P[ijet] > -1. ) varpos = 12.5*(1+Jet_cMVAv2P[ijet]);
        if ( Jet_cMVAv2N[ijet] > -1. ) varneg = 12.5*(1+Jet_cMVAv2N[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_cMVAv2[ijet] > TagCut[i] )   TagPos[i] = true;
          if ( Jet_cMVAv2N[ijet] > TagCut[i] ) TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
      }
*/

// CTag_CvsB
/*
      else if (TagName == "CvsB") {
        if ( CTag_Jet_CvsB[ijet] > -1. )  varPos = 12.5*(1+CTag_Jet_CvsB[ijet]);	//CvsB in -1, 1
        if ( CTag_Jet_CvsBP[ijet] > -1. ) varpos = 12.5*(1+CTag_Jet_CvsBP[ijet]);
        if ( CTag_Jet_CvsBN[ijet] > -1. ) varneg = 12.5*(1+CTag_Jet_CvsBN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( CTag_Jet_CvsB[ijet] > TagCut_CvsB[i] && CTag_Jet_CvsL[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( CTag_Jet_CvsBN[ijet] > TagCut_CvsB[i] && CTag_Jet_CvsLN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];

      } 
*/

// CTag_CvsL
/*
      else if (TagName == "CvsL") {
        if ( CTag_Jet_CvsL[ijet] > -1. )  varPos = 12.5*(1+CTag_Jet_CvsL[ijet]);		// Hopefully this is correct rescaling and test ?!
        if ( CTag_Jet_CvsLP[ijet] > -1. ) varpos = 12.5*(1+CTag_Jet_CvsLP[ijet]);
        if ( CTag_Jet_CvsLN[ijet] > -1. ) varneg = 12.5*(1+CTag_Jet_CvsLN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( CTag_Jet_CvsB[ijet] > TagCut_CvsB[i] && CTag_Jet_CvsL[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( CTag_Jet_CvsBN[ijet] > TagCut_CvsB[i] && CTag_Jet_CvsLN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true; 
        }
        cat = Jet_histJet[ijet];

      } 
*/
      else if ( TagName == "DeepCSVCvsB" )
      {
        if ( Jet_DeepCSVCvsBDisc[ijet] > 0. )  varPos = 25*(Jet_DeepCSVCvsBDisc[ijet]);	//CvsB in -1, 1
        //if ( Jet_DeepCSVCvsBDiscP[ijet] > 0. ) varpos = 25*(Jet_DeepCSVCvsBDiscP[ijet]);
        if ( Jet_DeepCSVCvsBDiscN[ijet] > 0. ) varneg = 25*(Jet_DeepCSVCvsBDiscN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepCSVCvsBDisc[ijet] > TagCut_CvsB[i] && Jet_DeepCSVCvsLDisc[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( Jet_DeepCSVCvsBDiscN[ijet] > TagCut_CvsB[i] && Jet_DeepCSVCvsLDiscN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
      }
      else if ( TagName == "DeepCSVCvsL" )
      {
        if ( Jet_DeepCSVCvsLDisc[ijet] > 0. )  varPos = 25*(Jet_DeepCSVCvsLDisc[ijet]);	// Hopefully this is correct rescaling and test ?!
        //if ( Jet_DeepCSVCvsLDiscP[ijet] > 0. ) varpos = 25*(Jet_DeepCSVCvsLDiscP[ijet]);
        if ( Jet_DeepCSVCvsLDiscN[ijet] > 0. ) varneg = 25*(Jet_DeepCSVCvsLDiscN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepCSVCvsBDisc[ijet] > TagCut_CvsB[i] && Jet_DeepCSVCvsLDisc[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( Jet_DeepCSVCvsBDiscN[ijet] > TagCut_CvsB[i] && Jet_DeepCSVCvsLDiscN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true; 
        }
        cat = Jet_histJet[ijet];
      }
      else if ( TagName == "DeepFlavourCvsB" )
      {
        if ( Jet_DeepFlavourCvsBDisc[ijet] > 0. )  varPos = 25*(Jet_DeepFlavourCvsBDisc[ijet]);	//CvsB in -1, 1
        //if ( Jet_DeepFlavourCvsBDiscP[ijet] > 0. ) varpos = 25*(Jet_DeepFlavourCvsBDiscP[ijet]);
        if ( Jet_DeepFlavourCvsBDiscN[ijet] > 0. ) varneg = 25*(Jet_DeepFlavourCvsBDiscN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepFlavourCvsBDisc[ijet] > TagCut_CvsB[i] && Jet_DeepFlavourCvsLDisc[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( Jet_DeepFlavourCvsBDiscN[ijet] > TagCut_CvsB[i] && Jet_DeepFlavourCvsLDiscN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true;
        }
        cat = Jet_histJet[ijet];
      }
      else if ( TagName == "DeepFlavourCvsL" )
      {
        if ( Jet_DeepFlavourCvsLDisc[ijet] > 0. )  varPos = 25*(Jet_DeepFlavourCvsLDisc[ijet]);	// Hopefully this is correct rescaling and test ?!
        //if ( Jet_DeepFlavourCvsLDiscP[ijet] > 0. ) varpos = 25*(Jet_DeepFlavourCvsLDiscP[ijet]);
        if ( Jet_DeepFlavourCvsLDiscN[ijet] > 0. ) varneg = 25*(Jet_DeepFlavourCvsLDiscN[ijet]);
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( Jet_DeepFlavourCvsBDisc[ijet] > TagCut_CvsB[i] && Jet_DeepFlavourCvsLDisc[ijet] > TagCut_CvsL[i])   TagPos[i] = true;
          if ( Jet_DeepFlavourCvsBDiscN[ijet] > TagCut_CvsB[i] && Jet_DeepFlavourCvsLDiscN[ijet] > TagCut_CvsL[i]) TagNeg[i] = true; 
        }
        cat = Jet_histJet[ijet];
      }

   //   cout << "\n CvsL: " << CTag_Jet_CvsL[ijet] << "  " << CTag_Jet_CvsLP[ijet] << "  " << CTag_Jet_CvsLN[ijet] ;
   //   cout << "\n CvsB: " << CTag_Jet_CvsB[ijet] << "  " << CTag_Jet_CvsBP[ijet] << "  " << CTag_Jet_CvsBN[ijet] ;


//$$
      if ( varPos > 24.9 ) varPos = 24.9; 
      if ( varpos > 24.9 ) varpos = 24.9; 
      if ( varneg > 24.9 ) varneg = 24.9; 

      if ( varPos < 0.1 ) varPos = 25*-0.5; 
      if ( varpos < 0.1 ) varpos = 25*-0.5; 
      if ( varneg < 0.1 ) varneg = 25*-0.5; 
//$$

// Track History information
      if      ( cat%e1 == 1		    || cat%e1 == 3)	   catP = 1; // BWeakDecay
      else if ((cat%e2 >= e1 && cat%e2 < 2*e1) || cat%e2 >= 3*e1 ) catP = 2; // CWeakDecay
      //else if ((cat%e3 >= e2 && cat%e3 < 2*e2) || cat%e3 >= 3*e2 ) catP = 3; // TauDecay
      else if ((cat%e4 >= e3 && cat%e4 < 2*e3) || cat%e4 >= 3*e3 ) catP = 4; // Conversion
      else if ((cat%e5 >= e4 && cat%e5 < 2*e4) || cat%e5 >= 3*e4 ) catP = 5; // KsDecay
      else if ((cat%e6 >= e5 && cat%e6 < 2*e5) || cat%e6 >= 3*e5 ) catP = 6; // LambdaDecay
      else if ((cat%e7 >= e6 && cat%e7 < 2*e6) || cat%e7 >= 3*e6 ) catP = 7; // Interaction
      else if ((cat%e8 >= e7 && cat%e8 < 2*e7) || cat%e8 >= 3*e7 ) catP = 8; // Fake
      if      ( cat%e1 >= 2 )    catN = 1; // BWeakDecay
      else if ( cat%e2 >= 2*e1 ) catN = 2; // CWeakDecay
      //else if ( cat%e3 >= 2*e2 ) catN = 3; // TauDecay
      else if ( cat%e4 >= 2*e3 ) catN = 4; // Conversion
      else if ( cat%e5 >= 2*e4 ) catN = 5; // KsDecay
      else if ( cat%e6 >= 2*e5 ) catN = 6; // LambdaDecay
      else if ( cat%e7 >= 2*e6 ) catN = 7; // Interaction
      else if ( cat%e8 >= 2*e7 ) catN = 8; // Fake

      if ( TagName == "TCHE" || TagName == "TCHP" ) {
        cat1 = Jet_hist1[ijet];
        if	     ( cat1%e1 == 1			|| cat1%e1 == 3)     cat1P = 1;
        else if ((cat1%e2 >= e1 && cat1%e2 < 2*e1) || cat1%e2 >= 3*e1 ) cat1P = 2;
        else if ((cat1%e3 >= e2 && cat1%e3 < 2*e2) || cat1%e3 >= 3*e2 ) cat1P = 3;
        else if ((cat1%e4 >= e3 && cat1%e4 < 2*e3) || cat1%e4 >= 3*e3 ) cat1P = 4;
        else if ((cat1%e5 >= e4 && cat1%e5 < 2*e4) || cat1%e5 >= 3*e4 ) cat1P = 5;
        else if ((cat1%e6 >= e5 && cat1%e6 < 2*e5) || cat1%e6 >= 3*e5 ) cat1P = 6;
        else if ((cat1%e7 >= e6 && cat1%e7 < 2*e6) || cat1%e7 >= 3*e6 ) cat1P = 7;
        else if ((cat1%e8 >= e7 && cat1%e8 < 2*e7) || cat1%e8 >= 3*e7 ) cat1P = 8;
        else if ((cat1    >= e8 && cat1    < 2*e8) || cat1    >= 3*e8 ) cat1P = 9;
        if	     ( cat1%e1 >= 2 )    cat1N = 1;
        else if ( cat1%e2 >= 2*e1 ) cat1N = 2;
        else if ( cat1%e3 >= 2*e2 ) cat1N = 3;
        else if ( cat1%e4 >= 2*e3 ) cat1N = 4;
        else if ( cat1%e5 >= 2*e4 ) cat1N = 5;
        else if ( cat1%e6 >= 2*e5 ) cat1N = 6;
        else if ( cat1%e7 >= 2*e6 ) cat1N = 7;
        else if ( cat1%e8 >= 2*e7 ) cat1N = 8;
        else if ( cat1    >= 2*e8 ) cat1N = 9;
        if ( cat1P > 0 && ( cat1P < catP || catP == 0 ) ) catP = cat1P; 
        if ( cat1N > 0 && ( cat1N < catN || catN == 0 ) ) catN = cat1N; 
      }

      if ( TagName == "TCHP" || TagName == "SSVHP" ) {
        if ( TagName == "TCHP" ) cat2 = Jet_hist2[ijet];
        if ( TagName == "SSVHP" ) cat2 = Jet_histSvx[ijet];
        if	     ( cat2%e1 == 1			|| cat2%e1 == 3)     cat2P = 1;
        else if ((cat2%e2 >= e1 && cat2%e2 < 2*e1) || cat2%e2 >= 3*e1 ) cat2P = 2;
        else if ((cat2%e3 >= e2 && cat2%e3 < 2*e2) || cat2%e3 >= 3*e2 ) cat2P = 3;
        else if ((cat2%e4 >= e3 && cat2%e4 < 2*e3) || cat2%e4 >= 3*e3 ) cat2P = 4;
        else if ((cat2%e5 >= e4 && cat2%e5 < 2*e4) || cat2%e5 >= 3*e4 ) cat2P = 5;
        else if ((cat2%e6 >= e5 && cat2%e6 < 2*e5) || cat2%e6 >= 3*e5 ) cat2P = 6;
        else if ((cat2%e7 >= e6 && cat2%e7 < 2*e6) || cat2%e7 >= 3*e6 ) cat2P = 7;
        else if ((cat2%e8 >= e7 && cat2%e8 < 2*e7) || cat2%e8 >= 3*e7 ) cat2P = 8;
        else if ((cat2    >= e8 && cat2    < 2*e8) || cat2    >= 3*e8 ) cat2P = 9;
        if	     ( cat2%e1 >= 2 )    cat2N = 1;
        else if ( cat2%e2 >= 2*e1 ) cat2N = 2;
        else if ( cat2%e3 >= 2*e2 ) cat2N = 3;
        else if ( cat2%e4 >= 2*e3 ) cat2N = 4;
        else if ( cat2%e5 >= 2*e4 ) cat2N = 5;
        else if ( cat2%e6 >= 2*e5 ) cat2N = 6;
        else if ( cat2%e7 >= 2*e6 ) cat2N = 7;
        else if ( cat2%e8 >= 2*e7 ) cat2N = 8;
        else if ( cat2    >= 2*e8 ) cat2N = 9;
        if ( TagName == "TCHP" ) {
          if ( cat2P > 0 && ( cat2P < catP || catP == 0 ) ) catP = cat2P; 
          if ( cat2N > 0 && ( cat2N < catN || catN == 0 ) ) catN = cat2N;
        } 
        //if ( TagName == "SSVHP" ) {
          //if ( Jet_Svx[ijet] > 1.  )  catP = cat2P; 
          //if ( Jet_SvxN[ijet] < -1. ) catN = cat2N;
        //} 
      }

      if      ( catP == 4 || catP == 7 ) GamP = true;
      else if ( catP == 5 )              K0sP = true;
      else if ( catP == 6 )              LamP = true;
      else if ( catP == 8 )              FakP = true;
      if      ( catP == 0 || catP == 9 ) OthP = true;
   
      if      ( catN == 4 || catN == 7 ) GamN = true;
      else if ( catN == 5 )              K0sN = true;
      else if ( catN == 6 )              LamN = true;
      else if ( catN == 8 )              FakN = true;
      if      ( catN == 0 || catN == 9 ) OthN = true;

      if      ( catP == 4 || catN == 4 || catP == 7 || catN == 7 ) GamPN = true;
      else if ( catP == 5 || catN == 5 )                           K0sPN = true;
      else if ( catP == 6 || catN == 6 )                           LamPN = true;
      else if ( catP == 8 || catN == 8 )                           FakPN = true;
      if      ( catP == 0 || catN == 0 || catP == 9 || catN == 9 ) OthPN = true;

// veto on positive tag
//    if ( Jet_Ip1P[ijet] > VetoPos ) Veto = true;
   
//*********************************

      hData_All_tracks->Fill( ntracks , 1. );
      hData_All_JetPV->Fill( npv , 1. );
      hData_All_JetPt->Fill( ptjet, etajet , 1. );
      hData_All_JetEta->Fill( etajet , 1. );
      hData_All_JetRun->Fill( Run , 1. );

      hAllFlav_All_Flavour->Fill( flavour , 1. );

      if ( flavour == 1 || flavour == 21 ) {
        hLightFlav_All_JetPV->Fill( npv , 1. );
        hLightFlav_All_JetPt->Fill( ptjet, etajet , 1. );
        hLightFlav_All_JetEta->Fill( etajet , 1. );
      }

      if (flavour == 21) {				// In principle this will now be empty !
        hGluonFlav_All_JetPV->Fill( npv , 1. );
        hGluonFlav_All_JetPt->Fill( ptjet, etajet , 1. );
        hGluonFlav_All_JetEta->Fill( etajet , 1. );
      }
      else if (flavour == 1) {				// In principle, it will contain all udsg jets !
        hUDSFlav_All_JetPV->Fill( npv , 1. );
        hUDSFlav_All_JetPt->Fill( ptjet, etajet , 1. );
        hUDSFlav_All_JetEta->Fill( etajet , 1. );
      }
      else if (flavour == 4) {				// C jets
        hCFlav_All_JetPV->Fill( npv , 1. );
        hCFlav_All_JetPt->Fill( ptjet, etajet , 1. );
        hCFlav_All_JetEta->Fill( etajet , 1. );
      }
      else if (flavour == 5) {				// B jets
        hBFlav_All_JetPV->Fill( npv , 1. );
        hBFlav_All_JetPt->Fill( ptjet, etajet , 1. );
        hBFlav_All_JetEta->Fill( etajet , 1. );
      }

      //PVH nseltracks
      for( unsigned int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
      {
        //TString hName = "hData_tracks_Pt_"+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt];
        if( ptjet >= fhseltrackPtMin[iPt] && ptjet < fhseltrackPtMax[iPt] ) hData_tracks[iPt]->Fill(nseltracks,tmp_ww);
      }
   
//*********************************
// Taggability
//$$
      if ( ntracks < NtrackMin ) continue;		// Go to next jet
//$$
      ntagjet++;

      hData_JetPV->Fill( npv , ww ) ;
      hData_JetPt->Fill( ptjet, etajet , ww );
      hData_JetEta->Fill( etajet , ww );
      hData_JetRun->Fill( Run , ww );
      if ( npv <= 5 ) {
        hData_JetPt_LE5pv->Fill( ptjet , ww );
        hData_JetEta_LE5pv->Fill( etajet , ww );
      } 
      else {
        hData_JetPt_GE6pv->Fill( ptjet , ww );
        hData_JetEta_GE6pv->Fill( etajet , ww );
      }
      if ( njet30 == 1 ) hData_1JetPt->Fill( ptjet , ww );
      if ( njet30 == 2 ) hData_2JetPt->Fill( ptjet , ww );
      if ( njet30 == 3 ) hData_3JetPt->Fill( ptjet , ww );
      if ( njet30 >= 4 ) hData_4JetPt->Fill( ptjet , ww );

      hAllFlav_Flavour->Fill( flavour , ww );
      if ( GamPN ) {
        hAllFlav_Gam_JetPV->Fill( npv , ww ) ;
        hAllFlav_Gam_JetPt->Fill( ptjet, etajet , ww );
        hAllFlav_Gam_JetEta->Fill( etajet , ww );
      }
      if ( K0sPN ) {
        hAllFlav_K0s_JetPV->Fill( npv , ww ) ;
        hAllFlav_K0s_JetPt->Fill( ptjet, etajet , ww );
        hAllFlav_K0s_JetEta->Fill( etajet , ww );
      }
      if ( LamPN ) {
        hAllFlav_Lam_JetPV->Fill( npv , ww ) ;
        hAllFlav_Lam_JetPt->Fill( ptjet, etajet , ww );
        hAllFlav_Lam_JetEta->Fill( etajet , ww );
      }
      if ( FakPN ) {
        hAllFlav_Fak_JetPV->Fill( npv , ww ) ;
        hAllFlav_Fak_JetPt->Fill( ptjet, etajet , ww );
        hAllFlav_Fak_JetEta->Fill( etajet , ww );
      }
   
      if ( flavour == 1 || flavour == 21 ) {
        hLightFlav_JetPU->Fill( npu , ww0 );
        hLightFlav_JetPV->Fill( npv , ww );
        hLightFlav_JetPt->Fill( ptjet, etajet , ww );
        hLightFlav_JetEta->Fill( etajet , ww );
        if ( nPU < 10 ) hLightFlav_JetPt_0pu->Fill( ptjet , ww );
        if ( nPU >= 20) hLightFlav_JetPt_GE8pu->Fill( ptjet , ww );
        if ( GamPN ) {
          hLightFlav_Gam_JetPV->Fill( npv , ww ) ;
          hLightFlav_Gam_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_Gam_JetEta->Fill( etajet , ww );
        }
        if ( K0sPN ) {
          hLightFlav_K0s_JetPV->Fill( npv , ww ) ;
          hLightFlav_K0s_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_K0s_JetEta->Fill( etajet , ww );
        }
        if ( LamPN ) {
          hLightFlav_Lam_JetPV->Fill( npv , ww ) ;
          hLightFlav_Lam_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_Lam_JetEta->Fill( etajet , ww );
        }
        if ( FakPN ) {
          hLightFlav_Fak_JetPV->Fill( npv , ww ) ;
          hLightFlav_Fak_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_Fak_JetEta->Fill( etajet , ww );
        }
        if ( OthPN ) {
          hLightFlav_Oth_JetPV->Fill( npv , ww ) ;
          hLightFlav_Oth_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_Oth_JetEta->Fill( etajet , ww );
        }
        if ( njet30 == 1 ) hLightFlav_1JetPt->Fill( ptjet , ww );
        if ( njet30 == 2 ) hLightFlav_2JetPt->Fill( ptjet , ww );
        if ( njet30 == 3 ) hLightFlav_3JetPt->Fill( ptjet , ww );
        if ( njet30 >= 4 ) hLightFlav_4JetPt->Fill( ptjet , ww );
      }

      if (flavour == 21) {
        hGluonFlav_JetPV->Fill( npv , ww );
        hGluonFlav_JetPt->Fill( ptjet, etajet , ww );
        hGluonFlav_JetEta->Fill( etajet , ww );
      }
      else if (flavour == 1) {
        hUDSFlav_JetPV->Fill( npv , ww );
        hUDSFlav_JetPt->Fill( ptjet, etajet , ww );
        hUDSFlav_JetEta->Fill( etajet , ww );
      }
      else if (flavour == 4) {
        hCFlav_JetPV->Fill( npv , ww );
        hCFlav_JetPt->Fill( ptjet, etajet , ww );
        hCFlav_JetEta->Fill( etajet , ww );
        if ( nPU < 10 ) hCFlav_JetPt_0pu->Fill( ptjet , ww );
        if ( nPU >= 20) hCFlav_JetPt_GE8pu->Fill( ptjet , ww );
      }
      else if (flavour == 5) {
        hBFlav_JetPU->Fill( npu , ww0 );
        hBFlav_JetPV->Fill( npv , ww );
        hBFlav_JetPt->Fill( ptjet, etajet , ww );
        if ( etajet < 1.2 ) hBFlav_JetPt_etaLT12->Fill( ptjet , ww );
        else                hBFlav_JetPt_etaGT12->Fill( ptjet , ww );
        hBFlav_JetEta->Fill( etajet , ww );
        if ( nPU < 10 ) hBFlav_JetPt_0pu->Fill( ptjet , ww );
        if ( nPU >= 20) hBFlav_JetPt_GE8pu->Fill( ptjet , ww );
      }
   
//*********************************
// Tagging

      for( unsigned int i = 0; i < 3; ++i )
      {
        if ( TagNeg[i] ) {
          hData_NegTag_JetPV[i]->Fill( npv , ww );
          hData_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
          hData_NegTag_JetEta[i]->Fill( etajet , ww );
          hData_NegTag_JetRun[i]->Fill( Run , ww );
          if ( npv <= 5 ) {
            hData_NegTag_JetPt_LE5pv[i]->Fill( ptjet , ww );
            hData_NegTag_JetEta_LE5pv[i]->Fill( etajet , ww );
          } 
          else {
            hData_NegTag_JetPt_GE6pv[i]->Fill( ptjet , ww );
            hData_NegTag_JetEta_GE6pv[i]->Fill( etajet , ww );
          }
          if ( njet30 == 1 ) hData_NegTag_1JetPt[i]->Fill( ptjet , ww );
          if ( njet30 == 2 ) hData_NegTag_2JetPt[i]->Fill( ptjet , ww );
          if ( njet30 == 3 ) hData_NegTag_3JetPt[i]->Fill( ptjet , ww );
          if ( njet30 >= 4 ) hData_NegTag_4JetPt[i]->Fill( ptjet , ww );
        } else break;
      }
      for( unsigned int i = 0; i < 3; ++i )
      {
        if ( TagPos[i] ) {
          hData_PosTag_JetPV[i]->Fill( npv , ww );
          hData_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
          hData_PosTag_JetEta[i]->Fill( etajet , ww );
          hData_PosTag_JetRun[i]->Fill( Run , ww );
          if ( K0sPN ) {
            hAllFlav_K0s_PosTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_K0s_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_K0s_PosTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( LamPN ) {
            hAllFlav_Lam_PosTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Lam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Lam_PosTag_JetEta[i]->Fill( etajet , ww );
          }
        } else break;
      }
      if ( varneg > 0 )
      {
        hData_Tagger->Fill(-varneg , etajet , ww );
        hData_PTagger->Fill( -varneg , etajet , ww );
      }
      if ( varpos > 0 )
      {
        hData_Tagger->Fill( varpos , etajet , ww );
        hData_MI2_Tagger->Fill( varpos , etajet , ww );
      }
      if ( varPos > 0 )
      {
        hData_PosTag->Fill( varPos , ww );
        hData_PTagger->Fill( varPos , etajet , ww );
        hData_MI2_PTagger->Fill( varPos , etajet , ww );
      }

      if( TagName == "DeepCSV" )
      {
/*
        if ( varpos > 0 )
        {
          hData_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hData_MI2_PTagger->Fill( varPos , etajet , ww );
        }
*/
      }
      if( TagName == "DeepFlavour" )
      {
      
/*
        if ( varpos > 0 )
        {
          hData_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hData_MI2_PTagger->Fill( varPos , etajet , ww );
        }
*/
      }
      for( unsigned int iEta = 0; iEta < NhEta; ++iEta )
      {
        for( unsigned int iPt = 0; iPt < NhPt; ++iPt )
        {
          if( etajet >= fhEtaMin[iEta] && etajet < fhEtaMax[iEta] && ptjet >= fhPtMin[iPt] && ptjet < fhPtMax[iPt] )
          {
            for( unsigned int iTagBin = 0; iTagBin <= (unsigned int) h_AllFlavour_All_IntTag[iEta][iPt]->GetNbinsX()+1; ++iTagBin )
            {
              Double_t xTagBin = h_AllFlavour_All_IntTag[iEta][iPt]->GetBinCenter(iTagBin);

              h_AllFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( K0sN || K0sP ) h_AllFlavour_K0s_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( LamN || LamP ) h_AllFlavour_Lam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( FakN || FakP ) h_AllFlavour_Fak_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( GamN || GamP ) h_AllFlavour_Gam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( OthN || OthP ) h_AllFlavour_Oth_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              if( flavour == 1 || flavour == 21 ) // UDS and Gluon
              {
                h_LightFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( K0sN || K0sP ) h_LightFlavour_K0s_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( LamN || LamP ) h_LightFlavour_Lam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( FakN || FakP ) h_LightFlavour_Fak_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( GamN || GamP ) h_LightFlavour_Gam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( OthN || OthP ) h_LightFlavour_Oth_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              }
              if( flavour == 1 ) // UDS only
              {
                h_UDSFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              }
              if( flavour == 21 ) // Gluon only
              {
                h_GluonFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              }
              if( flavour == 4 ) // C only
              {
                h_CFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              }
              if( flavour == 5 ) // B only
              {
                h_BFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
              }
            }
            if( varneg > 0 )
            {
              h_AllFlavour_All_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( K0sN ) h_AllFlavour_K0s_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( LamN ) h_AllFlavour_Lam_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( FakN ) h_AllFlavour_Fak_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( GamN ) h_AllFlavour_Gam_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( OthN ) h_AllFlavour_Oth_NegTag[iEta][iPt]->Fill( varneg, ww );
              if( flavour == 1 || flavour == 21 )
              {
                h_LightFlavour_All_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( K0sN ) h_LightFlavour_K0s_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( LamN ) h_LightFlavour_Lam_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( FakN ) h_LightFlavour_Fak_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( GamN ) h_LightFlavour_Gam_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( OthN ) h_LightFlavour_Oth_NegTag[iEta][iPt]->Fill( varneg, ww );
              }
              if( flavour == 1 ) // UDS only
              {
                h_UDSFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
              }
              if( flavour == 21 ) // Gluon only
              {
                h_GluonFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
              }
              if( flavour == 4 ) // C only
              {
                h_CFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
              }
              if( flavour == 5 ) // B only
              {
                h_BFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
              }
            }
            if( 1 ) // Used to be: varPos > 0 )
            {
              h_AllFlavour_All_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( K0sP ) h_AllFlavour_K0s_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( LamP ) h_AllFlavour_Lam_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( FakP ) h_AllFlavour_Fak_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( GamP ) h_AllFlavour_Gam_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( OthP ) h_AllFlavour_Oth_PosTag[iEta][iPt]->Fill( varPos, ww );
              if( flavour == 1 || flavour == 21 )
              {
                h_LightFlavour_All_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( K0sP ) h_LightFlavour_K0s_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( LamP ) h_LightFlavour_Lam_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( FakP ) h_LightFlavour_Fak_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( GamP ) h_LightFlavour_Gam_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( OthP ) h_LightFlavour_Oth_PosTag[iEta][iPt]->Fill( varPos, ww );
              }
              if( flavour == 1 ) // UDS only
              {
                h_UDSFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
              }
              if( flavour == 21 ) // Gluon only
              {
                h_GluonFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
              }
              if( flavour == 4 ) // C only
              {
                h_CFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
              }
              if( flavour == 5 ) // B only
              {
                h_BFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
              }
            }
          }
        }
      }
   
      if ( varneg > 0 ) {
        if      ( catN == 1 ) hAllFlav_Tagger_Bwd->Fill(-varneg , ww );
        else if ( catN == 2 ) hAllFlav_Tagger_Cwd->Fill(-varneg , ww );
        else if ( catN == 3 ) hAllFlav_Tagger_Tau->Fill(-varneg , ww );
        else if ( catN == 4 ) hAllFlav_Tagger_Gam->Fill(-varneg , ww );
        else if ( catN == 5 ) hAllFlav_Tagger_K0s->Fill(-varneg , ww );
        else if ( catN == 6 ) hAllFlav_Tagger_Lam->Fill(-varneg , ww );
        else if ( catN == 7 ) hAllFlav_Tagger_Int->Fill(-varneg , ww );
        else if ( catN == 8 ) hAllFlav_Tagger_Fak->Fill(-varneg , ww );
        else		   hAllFlav_Tagger_Oth->Fill(-varneg , ww );

        if      ( catP == 1 ) hAllFlav_PTagger_Bwd->Fill( -varneg , ww );
        else if ( catP == 2 ) hAllFlav_PTagger_Cwd->Fill( -varneg , ww );
        else if ( catP == 3 ) hAllFlav_PTagger_Tau->Fill( -varneg , ww );
        else if ( catP == 4 ) hAllFlav_PTagger_Gam->Fill( -varneg , ww );
        else if ( catP == 5 ) hAllFlav_PTagger_K0s->Fill( -varneg , ww );
        else if ( catP == 6 ) hAllFlav_PTagger_Lam->Fill( -varneg , ww );
        else if ( catP == 7 ) hAllFlav_PTagger_Int->Fill( -varneg , ww );
        else if ( catP == 8 ) hAllFlav_PTagger_Fak->Fill( -varneg , ww );
        else		   hAllFlav_PTagger_Oth->Fill( -varneg , ww );
      }
      if ( varpos > 0 ) {
        if      ( catP == 1 ) hAllFlav_Tagger_Bwd->Fill( varpos , ww );
        else if ( catP == 2 ) hAllFlav_Tagger_Cwd->Fill( varpos , ww );
        else if ( catP == 3 ) hAllFlav_Tagger_Tau->Fill( varpos , ww );
        else if ( catP == 4 ) hAllFlav_Tagger_Gam->Fill( varpos , ww );
        else if ( catP == 5 ) hAllFlav_Tagger_K0s->Fill( varpos , ww );
        else if ( catP == 6 ) hAllFlav_Tagger_Lam->Fill( varpos , ww );
        else if ( catP == 7 ) hAllFlav_Tagger_Int->Fill( varpos , ww );
        else if ( catP == 8 ) hAllFlav_Tagger_Fak->Fill( varpos , ww );
        else		   hAllFlav_Tagger_Oth->Fill( varpos , ww );
      }
      if ( varPos > 0 ) {
        if      ( catP == 1 ) hAllFlav_PTagger_Bwd->Fill( varPos , ww );
        else if ( catP == 2 ) hAllFlav_PTagger_Cwd->Fill( varPos , ww );
        else if ( catP == 3 ) hAllFlav_PTagger_Tau->Fill( varPos , ww );
        else if ( catP == 4 ) hAllFlav_PTagger_Gam->Fill( varPos , ww );
        else if ( catP == 5 ) hAllFlav_PTagger_K0s->Fill( varPos , ww );
        else if ( catP == 6 ) hAllFlav_PTagger_Lam->Fill( varPos , ww );
        else if ( catP == 7 ) hAllFlav_PTagger_Int->Fill( varPos , ww );
        else if ( catP == 8 ) hAllFlav_PTagger_Fak->Fill( varPos , ww );
        else		   hAllFlav_PTagger_Oth->Fill( varPos , ww );
      }
   
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
      for( unsigned int i = 0; i < 3; ++i )
      {
        if ( TagNeg[i] ) {
          hLightFlav_NegTag_JetPV[i]->Fill( npv , ww );
          hLightFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
          hLightFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          if ( FakN ) {
   	    hLightFlav_Fak_NegTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Fak_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Fak_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( OthN ) {
   	    hLightFlav_Oth_NegTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Oth_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Oth_NegTag_JetEta[i]->Fill( etajet , ww );
   	    hLightFlav_Oth_NegTag_nPU[i]->Fill( npu , ww );
   	    hLightFlav_Oth_NegTag_pthat[i]->Fill( pthat , ww );
          }
        } else break;
      }
      for( unsigned int i = 0; i < 3; ++i )
      {
        if ( TagPos[i] ) {
          hLightFlav_PosTag_JetPU[i]->Fill( npu , ww0 );
          hLightFlav_PosTag_JetPV[i]->Fill( npv , ww );
          hLightFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
          hLightFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          if ( catP >= 1 && catP <= 3 ) {
   	    hLightFlav_BCT_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_BCT_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_BCT_PosTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( GamP ) {
   	    hLightFlav_Gam_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Gam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Gam_PosTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( K0sP ) {
   	    hLightFlav_K0s_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_K0s_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_K0s_PosTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( LamP ) {
   	    hLightFlav_Lam_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Lam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Lam_PosTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( FakP ) {
   	    hLightFlav_Fak_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Fak_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Fak_PosTag_JetEta[i]->Fill( etajet , ww );
   	    hLightFlav_Fak_PosTag_nPU[i]->Fill( npu , ww );
   	    hLightFlav_Fak_PosTag_pthat[i]->Fill( pthat , ww );
          }
          if ( OthP ) {
   	    hLightFlav_Oth_PosTag_JetPV[i]->Fill( npv , ww );
   	    hLightFlav_Oth_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	    hLightFlav_Oth_PosTag_JetEta[i]->Fill( etajet , ww );
   	    hLightFlav_Oth_PosTag_nPU[i]->Fill( npu , ww );
   	    hLightFlav_Oth_PosTag_pthat[i]->Fill( pthat , ww );
          }
          if ( njet30 == 1 ) hLightFlav_PosTag_1JetPt[i]->Fill( ptjet , ww );
          if ( njet30 == 2 ) hLightFlav_PosTag_2JetPt[i]->Fill( ptjet , ww );
          if ( njet30 == 3 ) hLightFlav_PosTag_3JetPt[i]->Fill( ptjet , ww );
          if ( njet30 >= 4 ) hLightFlav_PosTag_4JetPt[i]->Fill( ptjet , ww );
        } else break;
      }
        if ( varneg > 0 )
        {
          hLightFlav_Tagger->Fill(-varneg , etajet , ww );
          hLightFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 )
        {
          hLightFlav_Tagger->Fill( varpos , etajet , ww );
          hLightFlav_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hLightFlav_PTagger->Fill( varPos , etajet , ww );
          hLightFlav_MI2_PTagger->Fill( varPos , etajet , ww );
        }
        if ( varPos > 0 ) hLightFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hLightFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hLightFlav_PosTagger_GE8pu->Fill( varPos , ww );
        if ( varneg > 0 ) {
          if      ( catN == 1 ) hLightFlav_Tagger_Bwd->Fill(-varneg , ww );
          else if ( catN == 2 ) hLightFlav_Tagger_Cwd->Fill(-varneg , ww );
          else if ( catN == 3 ) hLightFlav_Tagger_Tau->Fill(-varneg , ww );
          else if ( catN == 4 ) hLightFlav_Tagger_Gam->Fill(-varneg , ww );
          else if ( catN == 5 ) hLightFlav_Tagger_K0s->Fill(-varneg , ww );
          else if ( catN == 6 ) hLightFlav_Tagger_Lam->Fill(-varneg , ww );
          else if ( catN == 7 ) hLightFlav_Tagger_Int->Fill(-varneg , ww );
          else if ( catN == 8 ) hLightFlav_Tagger_Fak->Fill(-varneg , ww );
          else		     hLightFlav_Tagger_Oth->Fill(-varneg , ww );

          if      ( catP == 1 ) hLightFlav_PTagger_Bwd->Fill( -varneg , ww );
          else if ( catP == 2 ) hLightFlav_PTagger_Cwd->Fill( -varneg , ww );
          else if ( catP == 3 ) hLightFlav_PTagger_Tau->Fill( -varneg , ww );
          else if ( catP == 4 ) hLightFlav_PTagger_Gam->Fill( -varneg , ww );
          else if ( catP == 5 ) hLightFlav_PTagger_K0s->Fill( -varneg , ww );
          else if ( catP == 6 ) hLightFlav_PTagger_Lam->Fill( -varneg , ww );
          else if ( catP == 7 ) hLightFlav_PTagger_Int->Fill( -varneg , ww );
          else if ( catP == 8 ) hLightFlav_PTagger_Fak->Fill( -varneg , ww );
          else		     hLightFlav_PTagger_Oth->Fill( -varneg , ww );
        }
        if ( varpos > 0 ) {
          if      ( catP == 1 ) hLightFlav_Tagger_Bwd->Fill( varpos , ww );
          else if ( catP == 2 ) hLightFlav_Tagger_Cwd->Fill( varpos , ww );
          else if ( catP == 3 ) hLightFlav_Tagger_Tau->Fill( varpos , ww );
          else if ( catP == 4 ) hLightFlav_Tagger_Gam->Fill( varpos , ww );
          else if ( catP == 5 ) hLightFlav_Tagger_K0s->Fill( varpos , ww );
          else if ( catP == 6 ) hLightFlav_Tagger_Lam->Fill( varpos , ww );
          else if ( catP == 7 ) hLightFlav_Tagger_Int->Fill( varpos , ww );
          else if ( catP == 8 ) hLightFlav_Tagger_Fak->Fill( varpos , ww );
          else		     hLightFlav_Tagger_Oth->Fill( varpos , ww );
        }
        if ( varPos > 0 ) {
          if      ( catP == 1 ) hLightFlav_PTagger_Bwd->Fill( varPos , ww );
          else if ( catP == 2 ) hLightFlav_PTagger_Cwd->Fill( varPos , ww );
          else if ( catP == 3 ) hLightFlav_PTagger_Tau->Fill( varPos , ww );
          else if ( catP == 4 ) hLightFlav_PTagger_Gam->Fill( varPos , ww );
          else if ( catP == 5 ) hLightFlav_PTagger_K0s->Fill( varPos , ww );
          else if ( catP == 6 ) hLightFlav_PTagger_Lam->Fill( varPos , ww );
          else if ( catP == 7 ) hLightFlav_PTagger_Int->Fill( varPos , ww );
          else if ( catP == 8 ) hLightFlav_PTagger_Fak->Fill( varPos , ww );
          else		     hLightFlav_PTagger_Oth->Fill( varPos , ww );
        }
      } // light q+gluon

      if ( flavour == 21 ) { // gluon jets
        if ( varpos > 0 ) hGluonFlav_PosTag->Fill( varpos , ww );
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hGluonFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hGluonFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hGluonFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 )
        {
          hGluonFlav_Tagger->Fill(-varneg , etajet , ww );
          hGluonFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 ) hGluonFlav_Tagger->Fill( varpos , etajet , ww );
        if ( varPos > 0 ) hGluonFlav_PTagger->Fill( varPos , etajet , ww );
      } // gluon jets
   
      else if ( flavour == 1 ) { // uds jets
        if ( varpos > 0 ) hUDSFlav_PosTag->Fill( varpos , ww );
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hUDSFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hUDSFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hUDSFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 )
        {
          hUDSFlav_Tagger->Fill(-varneg , etajet , ww );
          hUDSFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 ) hUDSFlav_Tagger->Fill( varpos , etajet , ww );
        if ( varPos > 0 ) hUDSFlav_PTagger->Fill( varPos , etajet , ww );
      } // uds jets
   
      else if ( flavour == 4 ) { // c jets
        if ( varpos > 0 ) hCFlav_PosTag->Fill( varpos , ww );
        if ( varPos > 0 ) hCFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hCFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hCFlav_PosTagger_GE8pu->Fill( varPos , ww );
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hCFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hCFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hCFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 )
        {
          hCFlav_Tagger->Fill(-varneg , etajet , ww );
          hCFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 )
        {
          hCFlav_Tagger->Fill( varpos , etajet , ww );
          hCFlav_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hCFlav_PTagger->Fill( varPos , etajet , ww );
          hCFlav_MI2_PTagger->Fill( varPos , etajet , ww );
        }
      } // c jets
   
      else if ( flavour == 5 ) { // b jets
        if ( varpos > 0 ) hBFlav_PosTag->Fill( varpos , ww );
        if ( varPos > 0 ) hBFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hBFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hBFlav_PosTagger_GE8pu->Fill( varPos , ww );
        for( unsigned int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hBFlav_PosTag_JetPU[i]->Fill( npu , ww0 );
            hBFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hBFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hBFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 )
        {
          hBFlav_Tagger->Fill(-varneg , etajet , ww );
          hBFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 )
        {
          hBFlav_Tagger->Fill( varpos , etajet , ww );
          hBFlav_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hBFlav_PTagger->Fill( varPos , etajet , ww );
          hBFlav_MI2_PTagger->Fill( varPos , etajet , ww );
        }
      } // b jets
      else { // no flavour
        if ( varneg > 0 )
        {
          hNoFlav_Tagger->Fill(-varneg , etajet , ww );
          hNoFlav_PTagger->Fill( -varneg , etajet , ww );
        }
        if ( varpos > 0 )
        {
          hNoFlav_Tagger->Fill( varpos , etajet , ww );
          hNoFlav_MI2_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 )
        {
          hNoFlav_PTagger->Fill( varPos , etajet , ww );
          hNoFlav_MI2_PTagger->Fill( varPos , etajet , ww );
        }
      } // no flavour

//*********************************
// Negative Tag 
//$$
      for( unsigned int i = 0; i < 3; ++i )
      {
        if ( TagNeg[i] ) {

          if ( GamN ) {
            hAllFlav_Gam_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Gam_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Gam_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( K0sN ) {
            hAllFlav_K0s_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_K0s_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_K0s_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( LamN ) {
            hAllFlav_Lam_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Lam_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Lam_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( FakN ) {
            hAllFlav_Fak_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Fak_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Fak_NegTag_JetEta[i]->Fill( etajet , ww );
          }
       
          if (flavour == 21) {
            hGluonFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hGluonFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hGluonFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 1) {
            hUDSFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hUDSFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hUDSFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 4) {
            hCFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hCFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hCFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 5) {
            hBFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hBFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hBFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
        } else break;
      }

    } // end loop on jet
  
  } // end loop on events 

  std::cout << "total number of events in the files : " << NTreeEvents << std::endl ;   
  std::cout << "total events passing trigger cuts : " << allevents << std::endl ;   

//################################################

// Output Postscript

//   TCanvas* c = new TCanvas("c");
  hData_All_NJets -> Draw(); 
//   c->Print("output.ps(");
  
//################################################
  HistogramManager h ;
  
  h.WriteAllHistogramsInFile(filename.Data(),"recreate");
  h.DeleteAllHistograms();	// In order to avoid memory leaks
//################################################
}   
