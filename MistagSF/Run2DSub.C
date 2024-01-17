/*-------------------------------------------------------------------------------
*
* Parameter "Type" can be
*   -"DATA12A", "DATA12B" ...
*   -"JetTree_Pt-15to30", "JetTree_Pt-30to50"...
*   -"JetTree_Pt-15to30_TP", "JetTree_Pt-30to50_TP"...
*
* Parameter "tTagger" can be "All", "CSV", "CSVJP"...
*   if it is set to all, it will loop on all defined taggers.
*
*-------------------------------------------------------------------------------*/


#include "TString.h"
#include "TROOT.h"
#include <iostream>
#include "TChain.h"
#include "TFileCollection.h"
#include "TSystem.h"
#include "JetsAna2D.h"

//sampleType = "data","QCD","QCD_recodebug"
void Run2DSub(TString sampleType="data", TString samplePt="All", TString tTrigger="All", TString tTagger="All", bool WeightTracks = false, TString TrigType = "2022", TString period = "2022", bool runCondor=false)  //period is used to choose PU reweight file
{

 cout << "\n Input setting: " ;
 cout << "\n Sample type:   " << sampleType ;
 cout << "\n Sample_pt:     " << samplePt ;
 cout << "\n Trigger:       " << tTrigger ;
 cout << "\n Tagger:        " << tTagger ;
 cout << "\n Weight track:  " << WeightTracks ;
 cout << "\n TrigType:      " << TrigType ;
 cout << "\n Period:        " << period ;
 cout << "\n Run on condor: " << runCondor ;
 

 bool useInputFileList = false ;
 if (samplePt.Contains(".txt")) useInputFileList = true ; 

 gROOT->Reset() ; 

 // Other variables that should be FIX:

  vector<TString> TagNames;
  TagNames.push_back("DeepFlavour");
  TagNames.push_back("PNet");
  TagNames.push_back("ParT");

// -- General directory where we can find the NTUPLES
  TString ntDir="FileLists_2022/"; //in effect when not using file list as input
  if (period == "2022EE") ntDir = "FileLists_2022EE/";
// -- General directory for the outputs.
  TString oDir = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTagAnalyzer/Output_2022/";
  if (period == "2022EE") oDir = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTagAnalyzer/Output_2022EE/";
  if (runCondor) oDir = "./" ;

  TString iDir = "MC/QCD/" ;
  if (sampleType == "data") iDir = "Data/JetHT_JetMET/" ;   

  //FIXME
  TString weightPU_file = "Weights/nPV_weights_All.root" ;
  weightPU_file = "one" ; //no PU weight
  
  //FIXME
  TString weightPthat_file = "Weights/QCD_2022_fixPuppi.weightPthat" ;
  if (period == "2022EE")  weightPthat_file = "Weights/QCD_2022EE_fixPuppi.weightPthat" ;
  //weightPthat_file = "one"; //no PtHat applied

  
  if (sampleType == "data") {
     weightPU_file = "one" ;
     weightPthat_file = "one" ;
  }

  std::string trigger_ps_folder = "Weights/TriggerPS/" ;

  TString JSONFile = "one"; //no certification JSON used
  if(sampleType == "data") JSONFile = "Weights/Cert_Collisions2022_355100_362760_Golden.json"; 
  //FIXME
  JSONFile = "one"; //turn this off for now
 
  //NOTE: since run on condor PthatRanges do not affect things. If run locally need to fix PthatRange
  vector<TString> PthatRanges;
  if (samplePt == "All" || samplePt == "15to30")  PthatRanges.push_back("/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "30to50")  PthatRanges.push_back("/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "50to80")  PthatRanges.push_back("/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "80to120")  PthatRanges.push_back("/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "120to170")  PthatRanges.push_back("/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "170to300")  PthatRanges.push_back("/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "300to470")  PthatRanges.push_back("/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "470to600") PthatRanges.push_back("/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "600to800")  PthatRanges.push_back("/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "800to1000") PthatRanges.push_back("/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "1000to1400") PthatRanges.push_back("/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "1400to1800") PthatRanges.push_back("/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "1800to2400") PthatRanges.push_back("/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "2400to3200")  PthatRanges.push_back("/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/");
  if (samplePt == "All" || samplePt == "3200toInf")  PthatRanges.push_back("/QCD_Pt_3200toInf_TuneCP5_13TeV_pytlthia8/");

 /////////////for data PthatRanges is used as looping over data set period///////////////
//FIXME  
  if (sampleType == "data") {
    PthatRanges.clear() ;
    if (period == "2022") {
      if (samplePt == "All" || samplePt == "JetMETRun2022C-27Jun2023-v1") PthatRanges.push_back("JetMETRun2022C-27Jun2023-v1.txt") ;
      if (samplePt == "All" || samplePt == "JetHTRun2022C-27Jun2023-v2") PthatRanges.push_back("JetHTRun2022C-27Jun2023-v2.txt") ;
      if (samplePt == "All" || samplePt == "JetMETRun2022D-27Jun2023-v2") PthatRanges.push_back("JetMETRun2022D-27Jun2023-v2.txt") ;
    }
    if (period == "2022EE") {
      if (samplePt == "All" || samplePt == "JetMETRun2022G-PromptReco-v1") PthatRanges.push_back("JetMETRun2022G-PromptReco-v1.txt") ;
      if (samplePt == "All" || samplePt == "JetMETRun2022F-PromptReco-v1") PthatRanges.push_back("JetMETRun2022F-PromptReco-v1.txt") ;
    }
  }   


  if (sampleType == "data") {
    weightPU_file = "one" ;
    weightPthat_file = "one" ;
  }

  bool truePU = false;

  float maxCutJetPtMax = 1000000.;
  
  int TrigVal_12[11] =      {  0, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500 }; //500 will not use trigger prescale 
  TString STrigVal_12[11] = { "0", "40", "60", "80", "140", "200", "260", "320", "400", "450", "500"};
  float TrigCut_12[11] =    {  0., 50., 70., 100., 160., 220., 300., 360., 450., 500., 550.}; //not using trigger still use fatJet pT > 400 as minimum cut
  int NTriggers = 11 ; //FIXME
  
  int *TrigVal = TrigVal_12;
  TString *STrigVal = STrigVal_12;
  float *TrigCut = TrigCut_12;

  //FIXME, no trigger in MC?
  int iTrigMin = 0; 
  int iTrigMax = NTriggers;

  if (gROOT->GetClass("JetTree")==0) return;
  if (gROOT->GetClass("PS")==0) return;

 
  TChain c("btagana/ttree");

////////////////////////////////////////////////////////////////////////////////
// Main program start here
////////////////////////////////////////////////////////////////////////////////
//
  cout << "\n TagNames size: " << TagNames.size() ;
  cout << "\n trig: " << iTrigMin << "  " << iTrigMax ;
  cout << "\n Pthat: " << PthatRanges.size() ;
  for( int iCut = 0; iCut < TagNames.size(); ++iCut ) //loop over taggers
  { 
    for( int iTrig = iTrigMin ; iTrig < iTrigMax; ++iTrig )
    {    
      int nPtHat = PthatRanges.size() ;
      if (useInputFileList) nPtHat = 1 ;
      for( int iPthat = 0; iPthat < nPtHat; ++iPthat )
      {
        // Add RootTuple to the chain.
        c.Reset();  // But first reset it otherwise, you accumulate loop after loop.

       //////////////configure input sample//////////////////
        TString fileList = "" ;
        if (!useInputFileList) {
          TString tmp = PthatRanges.at(iPthat) ;
          tmp = tmp.ReplaceAll("/", "") ;
          fileList = ntDir + "/" + iDir + tmp + ".txt" ;
        }
        else { fileList = samplePt ; }
        cout << "\n File list is: " << fileList << endl ;
        TFileCollection fc("fc","list of input root files", fileList) ;
        c.AddFileInfoList((TCollection*)fc.GetList()) ;         

        /////////////configure output////////////////////////
        TString oFileBase = oDir;
        if (!runCondor) {
          oFileBase += iDir;  // use the same path as from the input directory
          oFileBase += "/2D/TrigType_" + TrigType;  // Add a subdir to differentiate various TrigType
          if (!useInputFileList) oFileBase += PthatRanges.at(iPthat);      // Add subdir for pthat ranges
          else oFileBase += "_" + samplePt ;
          gSystem->mkdir(oFileBase,kTRUE);	// Create output Directory if it does not exist
          cout << "Created dir" << endl;
        }
        if (sampleType != "data") {
          oFileBase += "/JetTree_mc_";
        }
        else {
          oFileBase += "/JetTree_data_";
        }
        if (runCondor) {
          TString tmp = samplePt ;
          tmp = tmp.ReplaceAll(".txt", "") ;
          oFileBase += tmp + "_" ;
        }

        oFileBase += period + "_" ;

        cout << "\n Number of entries of btagana: " << c.GetEntries() ;

        JetTree * t = new JetTree(&c);
        std::cout << "\n >>>>>>>>>: " << t << " " << t->fChain << " " << t->fChain->GetEntries();

        TString oFile = oFileBase;
        if( WeightTracks ) oFile += "TW_";	// For TrackWeight
        else oFile += "NW_";			// For NoWeight
        oFile += STrigVal[iTrig] + "_";	// Add the HLTCut in the fileName
        oFile = oFile + TagNames.at(iCut) + ".root";	// Add tagger type
        cout << "\n TagNames.at(iCut) = " << TagNames.at(iCut) << endl;
        cout << "tTagger = " << tTagger << endl;

        cout << "STrigVal[iTrig] = " << STrigVal[iTrig] << endl;
        cout << "tTrigger = " << tTrigger << endl;
        cout << "tTagger = " << tTagger << endl;
        cout << "TagNames.at(iCut) = " << TagNames.at(iCut) << endl;
        if( (TagNames.at(iCut) == tTagger || tTagger == "All") &&
            (STrigVal[iTrig] == tTrigger || tTrigger == "All") )
        {
          cout << "Trying to call loop with following parameters " << endl;
          cout << TrigVal[iTrig] << "\t" << TrigCut[iTrig] << endl;
          cout << oFile << endl;
          cout << "Json: " << JSONFile << "\tPileup: " << weightPU_file << "\tPtHat: " << weightPthat_file << "\tTruePU: " << truePU << endl;
          //
          //TEMP change to pt cut of 30
          cout << "\n 2. Chain t is " << t << " " << t->fChain->GetEntries() ;
          //FIXME trigger selection
          t->Loop(sampleType, 0,1000,TagNames.at(iCut), period, 20.,1000.,
                  TrigVal[iTrig],     //0: not apply trigger selection!!!!!!
                  0.,                 //not apply cuts on leading jet otherwise use TrigCut[iTrig]
                  maxCutJetPtMax, 
                  oFile, weightPU_file, weightPthat_file,
                  JSONFile, truePU, WeightTracks, TrigType) ;

          cout << "\n >>>>>>>>>> Done" << std::endl; 
        }

        delete t ;
      
      }
    }
  }
}
