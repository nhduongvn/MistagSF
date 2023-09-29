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
   Int_t           BitTrigger[D_nBitTrigger];
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           BX;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         PVez;
   Float_t         GenPVz;
   Float_t         nPUtrue;
   Int_t           nPU;
   Int_t           PU_bunch[D_nPU];
   Float_t         PV_x[D_nPV];
   Float_t         PV_y[D_nPV];
   Float_t         PV_chi2[D_nPV];
   Float_t         PV_ndf[D_nPV];
   Int_t           nJet;
   Float_t         Jet_pt[D_nJet];
   Float_t         Jet_uncorrpt[D_nJet];
   Float_t         Jet_genpt[D_nJet];
   Float_t         Jet_residual[D_nJet];
   Float_t         Jet_jes[D_nJet];
   Float_t         Jet_eta[D_nJet];
   Float_t         Jet_phi[D_nJet];
   Float_t         Jet_mass[D_nJet];
   Int_t           Jet_ntracks[D_nJet];
   Int_t           Jet_nseltracks[D_nJet];
   Int_t           Jet_flavour[D_nJet];
   Int_t           Jet_flavourCleaned[D_nJet];
   Float_t         Jet_ProbaN[D_nJet];
   Float_t         Jet_Proba[D_nJet];
   Float_t         Jet_BprobN[D_nJet];
   Float_t         Jet_Bprob[D_nJet];
   Float_t         Jet_CombIVF[D_nJet];
   Float_t         Jet_CombIVF_P[D_nJet];
   Float_t         Jet_CombIVF_N[D_nJet];
   Float_t         Jet_SoftMu[D_nJet];
   Float_t         Jet_SoftEl[D_nJet];
   Float_t         Jet_cMVAv2[D_nJet];
   Float_t         Jet_cMVAv2N[D_nJet];
   Float_t         Jet_cMVAv2P[D_nJet];
   Int_t           Jet_hist1[D_nJet];
   Int_t           Jet_hist2[D_nJet];
   Int_t           Jet_hist3[D_nJet];
   Int_t           Jet_histJet[D_nJet];
   Int_t           Jet_histSvx[D_nJet];
   Int_t           Jet_SV_multi[D_nJet];
   Float_t         Jet_DeepFlavourBDisc[D_nJet];
   Float_t         Jet_DeepFlavourCvsLDisc[D_nJet];
   Float_t         Jet_DeepFlavourCvsBDisc[D_nJet];
   Float_t         Jet_DeepFlavourBDiscN[D_nJet];
   Float_t         Jet_DeepFlavourCvsLDiscN[D_nJet];
   Float_t         Jet_DeepFlavourCvsBDiscN[D_nJet];
   Float_t         Jet_DeepFlavourPrunedBDisc[D_nJet];
   Float_t         Jet_DeepFlavourPrunedCvsLDisc[D_nJet];
   Float_t         Jet_DeepFlavourPrunedCvsBDisc[D_nJet];
   Float_t         Jet_DeepFlavourPrunedBDiscN[D_nJet];
   Float_t         Jet_DeepFlavourPrunedCvsLDiscN[D_nJet];
   Float_t         Jet_DeepFlavourPrunedCvsBDiscN[D_nJet];
   Float_t         Jet_DeepCSVBDisc[D_nJet];
   Float_t         Jet_DeepCSVBDiscN[D_nJet];
   Float_t         Jet_DeepCSVBDiscP[D_nJet];
   Float_t         Jet_DeepCSVCvsLDisc[D_nJet];
   Float_t         Jet_DeepCSVCvsLDiscN[D_nJet];
   Float_t         Jet_DeepCSVCvsLDiscP[D_nJet];
   Float_t         Jet_DeepCSVCvsBDisc[D_nJet];
   Float_t         Jet_DeepCSVCvsBDiscN[D_nJet];
   Float_t         Jet_DeepCSVCvsBDiscP[D_nJet];
   Int_t           Jet_nFirstSE[D_nJet];
   Int_t           Jet_nLastSE[D_nJet];
   Int_t           Jet_nFirstSM[D_nJet];
   Int_t           Jet_nLastSM[D_nJet];
   Int_t           Jet_nFirstTrack[D_nJet];
   Int_t           Jet_nLastTrack[D_nJet];
   Float_t         CTag_Jet_CvsB[D_nJet];
   Float_t         CTag_Jet_CvsBN[D_nJet];
   Float_t         CTag_Jet_CvsBP[D_nJet];
   Float_t         CTag_Jet_CvsL[D_nJet];
   Float_t         CTag_Jet_CvsLN[D_nJet];
   Float_t         CTag_Jet_CvsLP[D_nJet];
   Float_t         Jet_ptPruned[D_nJet];
   Float_t         Jet_etaPruned[D_nJet];
   Float_t         Jet_phiPruned[D_nJet];
   Float_t         Jet_massPruned[D_nJet];

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
   TBranch        *b_PVz;   //!
   TBranch        *b_PVez;   //!
   TBranch        *b_GenPVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_PU_bunch;   //!
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
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_nseltracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_flavourCleaned;   //!
   TBranch        *b_Jet_ProbaN;   //!
   TBranch        *b_Jet_Proba;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_Bprob;   //!
   TBranch        *b_Jet_CombIVF;   //!
   TBranch        *b_Jet_CombIVF_P;   //!
   TBranch        *b_Jet_CombIVF_N;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_SoftEl;   //!
   TBranch        *b_Jet_cMVAv2;   //!
   TBranch        *b_Jet_cMVAv2N;   //!
   TBranch        *b_Jet_cMVAv2P;   //!
   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_SV_multi;   //!
   TBranch        *b_Jet_DeepFlavourBDisc;   //!
   TBranch        *b_Jet_DeepFlavourCvsLDisc;   //!
   TBranch        *b_Jet_DeepFlavourCvsBDisc;   //!
   TBranch        *b_Jet_DeepFlavourBDiscN;   //!
   TBranch        *b_Jet_DeepFlavourCvsLDiscN;   //!
   TBranch        *b_Jet_DeepFlavourCvsBDiscN;   //!
   TBranch        *b_Jet_DeepFlavourPrunedBDisc;   //!
   TBranch        *b_Jet_DeepFlavourPrunedCvsLDisc;   //!
   TBranch        *b_Jet_DeepFlavourPrunedCvsBDisc;   //!
   TBranch        *b_Jet_DeepFlavourPrunedBDiscN;   //!
   TBranch        *b_Jet_DeepFlavourPrunedCvsLDiscN;   //!
   TBranch        *b_Jet_DeepFlavourPrunedCvsBDiscN;   //!
   TBranch        *b_Jet_DeepCSVBDisc;   //!
   TBranch        *b_Jet_DeepCSVBDiscN;   //!
   TBranch        *b_Jet_DeepCSVBDiscP;   //!
   TBranch        *b_Jet_DeepCSVCvsLDisc;   //!
   TBranch        *b_Jet_DeepCSVCvsLDiscN;   //!
   TBranch        *b_Jet_DeepCSVCvsLDiscP;   //!
   TBranch        *b_Jet_DeepCSVCvsBDisc;   //!
   TBranch        *b_Jet_DeepCSVCvsBDiscN;   //!
   TBranch        *b_Jet_DeepCSVCvsBDiscP;   //!
   TBranch        *b_Jet_nFirstSE;   //!
   TBranch        *b_Jet_nLastSE;   //!
   TBranch        *b_Jet_nFirstSM;   //!
   TBranch        *b_Jet_nLastSM;   //!
   TBranch        *b_Jet_nFirstTrack;   //!
   TBranch        *b_Jet_nLastTrack;   //!
   TBranch        *b_CTag_Jet_CvsB;   //!
   TBranch        *b_CTag_Jet_CvsBN;   //!
   TBranch        *b_CTag_Jet_CvsBP;   //!
   TBranch        *b_CTag_Jet_CvsL;   //!
   TBranch        *b_CTag_Jet_CvsLN;   //!
   TBranch        *b_CTag_Jet_CvsLP;   //!
   TBranch        *b_Jet_ptPruned;   //!
   TBranch        *b_Jet_etaPruned;   //!
   TBranch        *b_Jet_phiPruned;   //!
   TBranch        *b_Jet_massPruned;   //!

   JetTree(TTree *tree=0);
   virtual ~JetTree();
   //virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void   Loop(int minPV, int maxPV, TString TagName,
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

#ifdef JetTree_cxx
JetTree::JetTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../../NTUPLES/BTagAnalyzer_2017UL/MC/Pt_15to30/JetTree_mc_FatJets_Subjets_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../../NTUPLES/BTagAnalyzer_2017UL/MC/Pt_15to30/JetTree_mc_FatJets_Subjets_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../../../NTUPLES/BTagAnalyzer_2017UL/MC/Pt_15to30/JetTree_mc_FatJets_Subjets_1.root:/btagana");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

JetTree::~JetTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nBitTrigger", &nBitTrigger, &b_nBitTrigger);
   fChain->SetBranchAddress("BitTrigger", BitTrigger, &b_BitTrigger);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("BX", &BX, &b_BX);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVez", &PVez, &b_PVez);
   fChain->SetBranchAddress("GenPVz", &GenPVz, &b_GenPVz);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("PU_bunch", PU_bunch, &b_PU_bunch);
   fChain->SetBranchAddress("PV_x", PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_chi2", PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ndf", PV_ndf, &b_PV_ndf);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_uncorrpt", Jet_uncorrpt, &b_Jet_uncorrpt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_nseltracks", Jet_nseltracks, &b_Jet_nseltracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   fChain->SetBranchAddress("Jet_flavourCleaned", Jet_flavourCleaned, &b_Jet_flavourCleaned);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_Proba", Jet_Proba, &b_Jet_Proba);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_Bprob", Jet_Bprob, &b_Jet_Bprob);
   fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   fChain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P, &b_Jet_CombIVF_P);
   fChain->SetBranchAddress("Jet_CombIVF_N", Jet_CombIVF_N, &b_Jet_CombIVF_N);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl, &b_Jet_SoftEl);
   fChain->SetBranchAddress("Jet_cMVAv2", Jet_cMVAv2, &b_Jet_cMVAv2);
   fChain->SetBranchAddress("Jet_cMVAv2N", Jet_cMVAv2N, &b_Jet_cMVAv2N);
   fChain->SetBranchAddress("Jet_cMVAv2P", Jet_cMVAv2P, &b_Jet_cMVAv2P);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_SV_multi", Jet_SV_multi, &b_Jet_SV_multi);
   fChain->SetBranchAddress("Jet_DeepFlavourBDisc", Jet_DeepFlavourBDisc, &b_Jet_DeepFlavourBDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourCvsLDisc", Jet_DeepFlavourCvsLDisc, &b_Jet_DeepFlavourCvsLDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourCvsBDisc", Jet_DeepFlavourCvsBDisc, &b_Jet_DeepFlavourCvsBDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourBDiscN", Jet_DeepFlavourBDiscN, &b_Jet_DeepFlavourBDiscN);
   fChain->SetBranchAddress("Jet_DeepFlavourCvsLDiscN", Jet_DeepFlavourCvsLDiscN, &b_Jet_DeepFlavourCvsLDiscN);
   fChain->SetBranchAddress("Jet_DeepFlavourCvsBDiscN", Jet_DeepFlavourCvsBDiscN, &b_Jet_DeepFlavourCvsBDiscN);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedBDisc", Jet_DeepFlavourPrunedBDisc, &b_Jet_DeepFlavourPrunedBDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedCvsLDisc", Jet_DeepFlavourPrunedCvsLDisc, &b_Jet_DeepFlavourPrunedCvsLDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedCvsBDisc", Jet_DeepFlavourPrunedCvsBDisc, &b_Jet_DeepFlavourPrunedCvsBDisc);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedBDiscN", Jet_DeepFlavourPrunedBDiscN, &b_Jet_DeepFlavourPrunedBDiscN);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedCvsLDiscN", Jet_DeepFlavourPrunedCvsLDiscN, &b_Jet_DeepFlavourPrunedCvsLDiscN);
   fChain->SetBranchAddress("Jet_DeepFlavourPrunedCvsBDiscN", Jet_DeepFlavourPrunedCvsBDiscN, &b_Jet_DeepFlavourPrunedCvsBDiscN);
   fChain->SetBranchAddress("Jet_DeepCSVBDisc", Jet_DeepCSVBDisc, &b_Jet_DeepCSVBDisc);
   fChain->SetBranchAddress("Jet_DeepCSVBDiscN", Jet_DeepCSVBDiscN, &b_Jet_DeepCSVBDiscN);
   fChain->SetBranchAddress("Jet_DeepCSVBDiscP", Jet_DeepCSVBDiscP, &b_Jet_DeepCSVBDiscP);
   fChain->SetBranchAddress("Jet_DeepCSVCvsLDisc", Jet_DeepCSVCvsLDisc, &b_Jet_DeepCSVCvsLDisc);
   fChain->SetBranchAddress("Jet_DeepCSVCvsLDiscN", Jet_DeepCSVCvsLDiscN, &b_Jet_DeepCSVCvsLDiscN);
   fChain->SetBranchAddress("Jet_DeepCSVCvsLDiscP", Jet_DeepCSVCvsLDiscP, &b_Jet_DeepCSVCvsLDiscP);
   fChain->SetBranchAddress("Jet_DeepCSVCvsBDisc", Jet_DeepCSVCvsBDisc, &b_Jet_DeepCSVCvsBDisc);
   fChain->SetBranchAddress("Jet_DeepCSVCvsBDiscN", Jet_DeepCSVCvsBDiscN, &b_Jet_DeepCSVCvsBDiscN);
   fChain->SetBranchAddress("Jet_DeepCSVCvsBDiscP", Jet_DeepCSVCvsBDiscP, &b_Jet_DeepCSVCvsBDiscP);
   fChain->SetBranchAddress("Jet_nFirstSE", Jet_nFirstSE, &b_Jet_nFirstSE);
   fChain->SetBranchAddress("Jet_nLastSE", Jet_nLastSE, &b_Jet_nLastSE);
   fChain->SetBranchAddress("Jet_nFirstSM", Jet_nFirstSM, &b_Jet_nFirstSM);
   fChain->SetBranchAddress("Jet_nLastSM", Jet_nLastSM, &b_Jet_nLastSM);
   fChain->SetBranchAddress("Jet_nFirstTrack", Jet_nFirstTrack, &b_Jet_nFirstTrack);
   fChain->SetBranchAddress("Jet_nLastTrack", Jet_nLastTrack, &b_Jet_nLastTrack);
   fChain->SetBranchAddress("CTag_Jet_CvsB", CTag_Jet_CvsB, &b_CTag_Jet_CvsB);
   fChain->SetBranchAddress("CTag_Jet_CvsBN", CTag_Jet_CvsBN, &b_CTag_Jet_CvsBN);
   fChain->SetBranchAddress("CTag_Jet_CvsBP", CTag_Jet_CvsBP, &b_CTag_Jet_CvsBP);
   fChain->SetBranchAddress("CTag_Jet_CvsL", CTag_Jet_CvsL, &b_CTag_Jet_CvsL);
   fChain->SetBranchAddress("CTag_Jet_CvsLN", CTag_Jet_CvsLN, &b_CTag_Jet_CvsLN);
   fChain->SetBranchAddress("CTag_Jet_CvsLP", CTag_Jet_CvsLP, &b_CTag_Jet_CvsLP);
   fChain->SetBranchAddress("Jet_ptPruned", Jet_ptPruned, &b_Jet_ptPruned);
   fChain->SetBranchAddress("Jet_etaPruned", Jet_etaPruned, &b_Jet_etaPruned);
   fChain->SetBranchAddress("Jet_phiPruned", Jet_phiPruned, &b_Jet_phiPruned);
   fChain->SetBranchAddress("Jet_massPruned", Jet_massPruned, &b_Jet_massPruned);
   Notify();
}

Bool_t JetTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
//Int_t JetTree::Cut(Long64_t entry)
//{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   //return 1;
//}
#endif // #ifdef JetTree_cxx
