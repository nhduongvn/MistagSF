#ifndef LTANAPS_H
#define LTANAPS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#define NMAX 100000

//namespace LTANA
//{
   
   class PS
     {

      public:
	
	PS();
	virtual ~PS();
	
	void init(std::string path, bool doPS);
	
	double getPSWeight(std::string trigName,int run,int lb);

      private:
		
	void fillPS();

	int _psPFJet40[NMAX], _runPFJet40[NMAX], _lbPFJet40[NMAX];
	int _psPFJet60[NMAX], _runPFJet60[NMAX], _lbPFJet60[NMAX];
	int _psPFJet80[NMAX], _runPFJet80[NMAX], _lbPFJet80[NMAX];
	int _psPFJet140[NMAX], _runPFJet140[NMAX], _lbPFJet140[NMAX];
	int _psPFJet200[NMAX], _runPFJet200[NMAX], _lbPFJet200[NMAX];
	int _psPFJet260[NMAX], _runPFJet260[NMAX], _lbPFJet260[NMAX];
	int _psPFJet320[NMAX], _runPFJet320[NMAX], _lbPFJet320[NMAX];
	int _psPFJet400[NMAX], _runPFJet400[NMAX], _lbPFJet400[NMAX];
	int _psPFJet500[NMAX], _runPFJet500[NMAX], _lbPFJet500[NMAX];
	
	std::ifstream _fPFJet40;
	std::ifstream _fPFJet60;
	std::ifstream _fPFJet80;
	std::ifstream _fPFJet140;
	std::ifstream _fPFJet200;
	std::ifstream _fPFJet260;
	std::ifstream _fPFJet320;
	std::ifstream _fPFJet400;
	std::ifstream _fPFJet500;
	
	bool _doPS;
	int _nElem[100];
     };
//}

#endif
