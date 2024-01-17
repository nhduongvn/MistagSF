//need to change input and output directory around line 252
{
//  gSystem->Load("libFWCoreFWLite.so");

//  gSystem->Load("libCondFormatsBTauObjects.so") ;
//  gSystem->Load("libCondToolsBTau.so") ; 
//  AutoLibraryLoader::enable();
  gSystem->Load("Bases/BTagCalibrationStandalone_cc.so") ;
  gROOT->ProcessLine(".L mistag.C++g") ;
  TH1F retSF, retMistag ;
  LoopPlot(retSF, retMistag, "DeepFlavour", "2022EE", 4.) ;
  LoopPlot(retSF, retMistag, "PNet", "2022EE", 4.) ;
  LoopPlot(retSF, retMistag, "ParT", "2022EE", 4.) ;
}

//this is a test
