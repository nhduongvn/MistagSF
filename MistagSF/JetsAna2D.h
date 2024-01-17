//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 17:42:51 2020 by ROOT version 5.34/18
// from TTree ttree/ttree
// found on file: ../../../NTUPLES/BTagAnalyzer_2017UL/MC/Pt_15to30/JetTree_mc_FatJets_Subjets_1.root
//////////////////////////////////////////////////////////

#ifndef JetTree_h
#define JetTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define D_nBitTrigger 100
#define D_nPU 100
#define D_nPV 100
#define D_ncQuarks 1000
#define D_nbQuarks 1000
#define D_nBHadrons 1000
#define D_nDHadrons 1000
#define D_nDaughters 1000
#define D_nGenlep 1000
#define D_nGenquark 1000
#define D_nGenPruned 1000
#define D_nGenV0 1000
#define D_nJet 1000
#define D_nSVTagVar 100
#define D_nPFElectron 1000
#define D_nPFMuon 1000
#define D_nTrkInc 1000
#define D_nTrkTagVarCSV 1000
#define D_nTrkCTagVar 1000
#define D_nTrkEtaRelTagVarCSV 1000
#define D_nTrkEtaRelCTagVar 1000
#define D_nLeptons 1000


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nBitTrigger;
   Long64_t        BitTrigger[D_nBitTrigger];
   UInt_t          Run;
   //Int_t           Run;
   ULong64_t       Evt;
   UInt_t          LumiBlock;
   //Int_t          LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           BX;
   Int_t           nPV;
   Int_t           nPVGood;
   Float_t         PVz;
   Float_t         PVez;
   Float_t         GenPVz;
   Float_t         nPUtrue;
   Int_t           nPU;
   Float_t         PV_x[D_nPV];
   Float_t         PV_y[D_nPV];
   Float_t         PV_chi2[D_nPV];
   Float_t         PV_ndf[D_nPV];
   Int_t           nJet;
   Float_t         Jet_pt[D_nJet];
   Float_t         Jet_uncorrpt[D_nJet];
   Double_t        Jet_genpt[D_nJet];
   Float_t         Jet_eta[D_nJet];
   Float_t         Jet_phi[D_nJet];
   Long64_t        Jet_flavour[D_nJet];
   Long64_t        Jet_flavourCleaned[D_nJet];
   //Int_t           Jet_flavour[D_nJet];
   //Int_t           Jet_flavourCleaned[D_nJet];
   Float_t         Jet_DeepFlavourBDisc[D_nJet];
   Float_t         Jet_DeepFlavourBDiscN[D_nJet];
   Float_t         Jet_PNetBDisc[D_nJet];
   Float_t         Jet_PNetBDiscN[D_nJet];
   Float_t         Jet_ParTBDisc[D_nJet];
   Float_t         Jet_ParTBDiscN[D_nJet];
   Int_t           Jet_nFirstSE[D_nJet];
   Int_t           Jet_nLastSE[D_nJet];
   Int_t           Jet_nFirstSM[D_nJet];
   Int_t           Jet_nLastSM[D_nJet];
   Int_t           Jet_nFirstTrack[D_nJet];
   Int_t           Jet_nLastTrack[D_nJet];
   Bool_t          Jet_tightID[D_nJet];
   Bool_t          Jet_tightlepvetoID[D_nJet];
   Float_t         Jet_vetomap[D_nJet];

   // List of branches
   TBranch        *b_nBitTrigger;   //!
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPVGood;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVez;   //!
   TBranch        *b_GenPVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ndf;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_uncorrpt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_nseltracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_flavourCleaned;   //!
   TBranch        *b_Jet_DeepFlavourBDisc;   //!
   TBranch        *b_Jet_DeepFlavourBDiscN;   //!
   TBranch        *b_Jet_PNetBDisc;   //!
   TBranch        *b_Jet_PNetBDiscN;   //!
   TBranch        *b_Jet_ParTBDisc;   //!
   TBranch        *b_Jet_ParTBDiscN;   //!
   TBranch        *b_Jet_nFirstSE;   //!
   TBranch        *b_Jet_nLastSE;   //!
   TBranch        *b_Jet_nFirstSM;   //!
   TBranch        *b_Jet_nLastSM;   //!
   TBranch        *b_Jet_tightlepvetoID;
   TBranch        *b_Jet_tightID;
   TBranch        *b_Jet_vetomap;

   JetTree(TTree *tree=0);
   virtual ~JetTree();
   //virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void   Loop(TString sampleType, int minPV, int maxPV, TString TagName, TString TagCutSetName,
               Float_t  aPtMin = 20.,
               Float_t  aPtMax = 1000.,
               //Float_t  aFreeCut = 0.,
               Int_t    aIntCut = 0,
               float    minCutJetPtMax = -1,
               float maxCutJetPtMax = 999999,
               TString  afilename = "output.root",
               TString weightPU_file = "",
               TString weightPthat_file = "",
               TString JSONFile = "",
               bool truePU = true,
               bool WeightTracks = true,
               TString TrigType = "2011");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
