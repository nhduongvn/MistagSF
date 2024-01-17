#include "PS.h"

#include "TSystem.h"

//LTANA::PS::PS()
PS::PS()
{
}

//LTANA::PS::~PS()
PS::~PS()
{
   if( _doPS )
     {
	_fPFJet40.close();
	_fPFJet60.close();
	_fPFJet80.close();
	_fPFJet140.close();
	_fPFJet200.close();
	_fPFJet260.close();
	_fPFJet320.close();
	_fPFJet400.close();
	_fPFJet500.close();
     }   	     
}

//void LTANA::PS::init(std::string path, bool doPS)
void PS::init(std::string path, bool doPS)
{
   std::cout << "Trigger prescale folder in PS::init: " << path << std::endl ;
   std::string fnamePFJet40 = path+"/HLT_PFJet40.csv";
   std::string fnamePFJet60 = path+"/HLT_PFJet60.csv";
   std::string fnamePFJet80 = path+"/HLT_PFJet80.csv";
   std::string fnamePFJet140 = path+"/HLT_PFJet140.csv";
   std::string fnamePFJet200 = path+"/HLT_PFJet200.csv";
   std::string fnamePFJet260 = path+"/HLT_PFJet260.csv";
   std::string fnamePFJet320 = path+"/HLT_PFJet320.csv";
   std::string fnamePFJet400 = path+"/HLT_PFJet400.csv";
   std::string fnamePFJet500 = path+"HLT_PFJet500.csv";
   
   bool existPFJet40 = !(gSystem->AccessPathName(fnamePFJet40.c_str()));
   bool existPFJet60 = !(gSystem->AccessPathName(fnamePFJet60.c_str()));
   bool existPFJet80 = !(gSystem->AccessPathName(fnamePFJet80.c_str()));
   bool existPFJet140 = !(gSystem->AccessPathName(fnamePFJet140.c_str()));
   bool existPFJet200 = !(gSystem->AccessPathName(fnamePFJet200.c_str()));
   bool existPFJet260 = !(gSystem->AccessPathName(fnamePFJet260.c_str()));
   bool existPFJet320 = !(gSystem->AccessPathName(fnamePFJet320.c_str()));
   bool existPFJet400 = !(gSystem->AccessPathName(fnamePFJet400.c_str()));
   bool existPFJet500 = !(gSystem->AccessPathName(fnamePFJet500.c_str()));

   if( !existPFJet40  || !existPFJet60  || !existPFJet80 || !existPFJet140
    || !existPFJet200 || !existPFJet260 || !existPFJet320 || !existPFJet400 || !existPFJet500)
     {
	std::cout << "Input file with PS factors can not be opened, no reweighting will be applied" << std::endl;
	doPS = 0;
     }
   else 
     {
	std::cout << "Read PS" << std::endl;
	_fPFJet40.open(fnamePFJet40.c_str());
	_fPFJet60.open(fnamePFJet60.c_str());
	_fPFJet80.open(fnamePFJet80.c_str());
	_fPFJet140.open(fnamePFJet140.c_str());
	_fPFJet200.open(fnamePFJet200.c_str());
	_fPFJet260.open(fnamePFJet260.c_str());
	_fPFJet320.open(fnamePFJet320.c_str());
	_fPFJet400.open(fnamePFJet400.c_str());
	_fPFJet500.open(fnamePFJet500.c_str());
     }   
   
   _doPS = doPS;

   for(int i=0;i<100;i++) _nElem[i] = 0;
   
   for(int i=0;i<NMAX;i++)
     {	
	_psPFJet40[i] = 0.;
	_runPFJet40[i] = 0.;
	_lbPFJet40[i] = 0.;
	
	_psPFJet60[i] = 0.;
	_runPFJet60[i] = 0.;
	_lbPFJet60[i] = 0.;

	_psPFJet80[i] = 0.;
	_runPFJet80[i] = 0.;
	_lbPFJet80[i] = 0.;

	_psPFJet140[i] = 0.;
	_runPFJet140[i] = 0.;
	_lbPFJet140[i] = 0.;

	_psPFJet200[i] = 0.;
	_runPFJet200[i] = 0.;
	_lbPFJet200[i] = 0.;

	_psPFJet260[i] = 0.;
	_runPFJet260[i] = 0.;
	_lbPFJet260[i] = 0.;

	_psPFJet320[i] = 0.;
	_runPFJet320[i] = 0.;
	_lbPFJet320[i] = 0.;

	_psPFJet400[i] = 0.;
	_runPFJet400[i] = 0.;
	_lbPFJet400[i] = 0.;
	
  _psPFJet500[i] = 0.;
	_runPFJet500[i] = 0.;
	_lbPFJet500[i] = 0.;
     }   
   
   if( _doPS ) fillPS();
   std::cout << "\n End fillPS " ;
}

//void LTANA::PS::fillPS()
void PS::fillPS()
{
   std::string trigName[9] = {"PFJet40",
      "PFJet60","PFJet80","PFJet140",
      "PFJet200","PFJet260","PFJet320","PFJet400","PFJet500"};
   
   int run,cmsls,totprescval;
   int prev_run, prev_cmsls;
   
   for(int i=0;i<9;i++)
     {	
	std::string name = trigName[i];
	
	std::ifstream *fps;
	if( strcmp(name.c_str(),"PFJet40") == 0 ) fps = &_fPFJet40;
	else if( strcmp(name.c_str(),"PFJet60") == 0 ) fps = &_fPFJet60;
	else if( strcmp(name.c_str(),"PFJet80") == 0 ) fps = &_fPFJet80;
	else if( strcmp(name.c_str(),"PFJet140") == 0 ) fps = &_fPFJet140;
	else if( strcmp(name.c_str(),"PFJet200") == 0 ) fps = &_fPFJet200;
	else if( strcmp(name.c_str(),"PFJet260") == 0 ) fps = &_fPFJet260;
	else if( strcmp(name.c_str(),"PFJet320") == 0 ) fps = &_fPFJet320;
	else if( strcmp(name.c_str(),"PFJet400") == 0 ) fps = &_fPFJet400;
	else if( strcmp(name.c_str(),"PFJet500") == 0 ) fps = &_fPFJet500;
	else
	  {
	     std::cout << "PS file not known: " << name << std::endl;
	     exit(1);
	  }	
	
	int c = 0;

        prev_run = 0;
        prev_cmsls = 0;
	while( !fps->eof() )
	  {
	     std::string line;
	     
	     *fps >> line;

	     if( c > 0 )
	       {		  
		  std::stringstream ss(line);
		  std::string item;
		  
		  int idx = 0;
		  while( std::getline(ss,item,',') )
		    {
		       if( idx == 0 ) run = atoi(item.c_str());
		       if( idx == 1 ) cmsls = atoi(item.c_str());
		       if( idx == 3 ) totprescval = atoi(item.c_str());
                       if( run == prev_run )
                       {
                         if( cmsls < prev_cmsls )
                         {
                           cout << "Error: run = " << run << " and cmsls = " << cmsls << endl;
                           cout << "while previous values were: run = " << prev_run << " and cmsls = " << prev_cmsls << endl;
                         }
                       }
                       if( run < prev_run )
                       {
                           //cout << "\tWarning: run = " << run << " and cmsls = " << cmsls << endl;
                           //cout << "\twhile previous values were: run = " << prev_run << " and cmsls = " << prev_cmsls << endl;
                       }
          
		       idx++;
		       
		       if( strcmp(name.c_str(),"PFJet40") == 0 )
			 {		       
			    _psPFJet40[c] = totprescval;
			    _runPFJet40[c] = run;
			    _lbPFJet40[c] = cmsls;
			 }		  
		       else if( strcmp(name.c_str(),"PFJet60") == 0 )
			 {		       
			    _psPFJet60[c] = totprescval;
			    _runPFJet60[c] = run;
			    _lbPFJet60[c] = cmsls;
			 }		  
		       else if( strcmp(name.c_str(),"PFJet80") == 0 )
			 {		       
			    _psPFJet80[c] = totprescval;
			    _runPFJet80[c] = run;
			    _lbPFJet80[c] = cmsls;
			 }		  
		       else if( strcmp(name.c_str(),"PFJet140") == 0 )
			 {		       
			    _psPFJet140[c] = totprescval;
			    _runPFJet140[c] = run;
			    _lbPFJet140[c] = cmsls;
			 }		  
		       else if( strcmp(name.c_str(),"PFJet200") == 0 )
			 {		       
			    _psPFJet200[c] = totprescval;
			    _runPFJet200[c] = run;
			    _lbPFJet200[c] = cmsls;
                         }
		       else if( strcmp(name.c_str(),"PFJet260") == 0 )
			 {		       
			    _psPFJet260[c] = totprescval;
			    _runPFJet260[c] = run;
			    _lbPFJet260[c] = cmsls;
			 }
		       else if( strcmp(name.c_str(),"PFJet320") == 0 )
			 {		       
			    _psPFJet320[c] = totprescval;
			    _runPFJet320[c] = run;
			    _lbPFJet320[c] = cmsls;
			 }
		       else if( strcmp(name.c_str(),"PFJet400") == 0 )
			 {		       
			    _psPFJet400[c] = totprescval;
			    _runPFJet400[c] = run;
			    _lbPFJet400[c] = cmsls;
			 }
		       else if( strcmp(name.c_str(),"PFJet500") == 0 )
			 {		       
			    _psPFJet500[c] = totprescval;
			    _runPFJet500[c] = run;
			    _lbPFJet500[c] = cmsls;
			 }
		    }
                    prev_run = run;
                    prev_cmsls = cmsls;
                 //cout << "File " << name.c_str() <<  " run = " << run << " ; cmsls = " << cmsls << " ; totps = " << totprescval << endl;
	       }
	     
	     c++;
	  }

	_nElem[i] = c;
     }   
}

double PS::getPSWeight(std::string trigName,int run,int lb)
{
   double w = 1.;
   
   if( _doPS )
     {
	if( strcmp(trigName.c_str(),"PFJet40") == 0 )
	  {
	     for(int i=0;i<_nElem[0];i++)
	       {
		  if( run == _runPFJet40[i] &&
		      lb >= _lbPFJet40[i] )
		    {
		       int j = i;
		       while( j < _nElem[0]-1 )
			 {			    
			    if( !(run == _runPFJet40[j+1] &&
				  lb >= _lbPFJet40[j+1] &&
				   _lbPFJet40[j+1] > _lbPFJet40[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet40[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet60") == 0 )
	  {
	     for(int i=0;i<_nElem[1];i++)
	       {
		  if( run == _runPFJet60[i] &&
		      lb >= _lbPFJet60[i] )
		    {
		       int j = i;
		       while( j < _nElem[1]-1 )
			 {			    
			    if( !(run == _runPFJet60[j+1] &&
				  lb >= _lbPFJet60[j+1]   &&
				   _lbPFJet60[j+1] > _lbPFJet60[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet60[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet80") == 0 )
	  {
	     for(int i=0;i<_nElem[2];i++)
	       {
		  if( run == _runPFJet80[i] &&
		      lb >= _lbPFJet80[i] )
		    {
		       int j = i;
		       while( j < _nElem[2]-1 )
			 {			    
			    if( !(run == _runPFJet80[j+1] &&
				  lb >= _lbPFJet80[j+1]   &&
				   _lbPFJet80[j+1] > _lbPFJet80[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet80[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet140") == 0 )
	  {
	     for(int i=0;i<_nElem[3];i++)
	       {
		  if( run == _runPFJet140[i] &&
		      lb >= _lbPFJet140[i] )
		    {
		       int j = i;
		       while( j < _nElem[3]-1 )
			 {			    
			    if( !(run == _runPFJet140[j+1] &&
				  lb >= _lbPFJet140[j+1]   &&
				   _lbPFJet140[j+1] > _lbPFJet140[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet140[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet200") == 0 )
	  {
	     for(int i=0;i<_nElem[4];i++)
	       {
		  if( run == _runPFJet200[i] &&
		      lb >= _lbPFJet200[i] )
		    {
		       int j = i;
		       while( j < _nElem[4]-1 )
			 {			    
			    if( !(run == _runPFJet200[j+1] &&
				  lb >= _lbPFJet200[j+1]   &&
				   _lbPFJet200[j+1] > _lbPFJet200[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet200[j];
		       
		       break;
		    }		  
	       }	     
	  }
	else if( strcmp(trigName.c_str(),"PFJet260") == 0 )
	  {
	     for(int i=0;i<_nElem[5];i++)
	       {
		  if( run == _runPFJet260[i] &&
		      lb >= _lbPFJet260[i] )
		    {
		       int j = i;
		       while( j < _nElem[5]-1 )
			 {			    
			    if( !(run == _runPFJet260[j+1] &&
				  lb >= _lbPFJet260[j+1]   &&
				   _lbPFJet260[j+1] > _lbPFJet260[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet260[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet320") == 0 )
	  {
	     for(int i=0;i<_nElem[6];i++)
	       {
		  if( run == _runPFJet320[i] &&
		      lb >= _lbPFJet320[i] )
		    {
		       int j = i;
		       while( j < _nElem[6]-1 )
			 {			    
			    if( !(run == _runPFJet320[j+1] &&
				  lb >= _lbPFJet320[j+1]  &&
				   _lbPFJet320[j+1] > _lbPFJet320[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet320[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet400") == 0 )
	  {
	     for(int i=0;i<_nElem[7];i++)
	       {
		  if( run == _runPFJet400[i] &&
		      lb >= _lbPFJet400[i] )
		    {
		       int j = i;
		       while( j < _nElem[7]-1 )
			 {			    
			    if( !(run == _runPFJet400[j+1] &&
				  lb >= _lbPFJet400[j+1]  &&
				   _lbPFJet400[j+1] > _lbPFJet400[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet400[j];
		       
		       break;
		    }		  
	       }	     
	  }	
	else if( strcmp(trigName.c_str(),"PFJet500") == 0 )
	  {
	     for(int i=0;i<_nElem[8];i++)
	       {
		  if( run == _runPFJet500[i] &&
		      lb >= _lbPFJet500[i] )
		    {
		       int j = i;
		       while( j < _nElem[8]-1 )
			 {			    
			    if( !(run == _runPFJet500[j+1] &&
				  lb >= _lbPFJet500[j+1]  &&
				   _lbPFJet500[j+1] > _lbPFJet500[i]) )
			      break;
			    j = i++;
			 }
		       
		       w = _psPFJet500[j];
		       
		       break;
		    }		  
	       }	     
	  }	
     }
   
   return w;
}
