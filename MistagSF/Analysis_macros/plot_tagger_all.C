{
  gROOT->ProcessLine(".L plot_tagger.C++") ;
  gROOT->ProcessLine(".L plot_tagger_C.so") ;
  TString period = "2022";
  plot("DeepFlavour", "0", 1,  25, 50, period, 0.1, 0.4) ;
  plot("DeepFlavour", "40", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "60", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "80", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "140", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "200", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "260", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "320", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "400", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "450", 1,  25, 50, period, 0.1, 0.3) ;
  plot("DeepFlavour", "500", 1,  25, 50, period, 0.1, 0.3) ;

  plot("PNet", "0", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "40", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "60", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "80", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "140", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "200", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "260", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "320", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "400", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "450", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "500", 1,  25, 50, period, 0.1, 0.3) ;
  
  plot("ParT", "0", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "40", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "60", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "80", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "140", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "200", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "260", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "320", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "400", 1,  25, 50, period, 0.1, 0.3) ;
  plot("ParT", "450", 1,  25, 50, period, 0.1, 0.3) ;
  plot("PNet", "500", 1,  25, 50, period, 0.1, 0.3) ;
}
