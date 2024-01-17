#ifndef BTagEntryN_H
#define BTagEntryN_H
/**
*
* BTagEntryN
*
* Represents one pt- or discriminator-dependent calibration function.
*
* measurement_type: e.g. comb, ttbar, di-mu, boosted, ...
* sys_type: e.g. central, plus, minus, plus_JEC, plus_JER, ...
*
* Everything is converted into a function, as it is easiest to store it in a
* txt or json file.
*
************************************************************/
#include <string>
#include <TF1.h>
#include <TH1.h>
class BTagEntryN
{
public:
enum OperatingPoint {
OP_LOOSE=0,
OP_MEDIUM=1,
OP_TIGHT=2,
OP_XTIGHT=3,
OP_XXTIGHT=4,
OP_RESHAPING=5,
};
enum JetFlavor {
FLAV_B=5,
FLAV_C=4,
FLAV_UDSG=0,
};
struct Parameters {
OperatingPoint operatingPoint;
std::string measurementType;
std::string sysType;
JetFlavor jetFlavor;
float etaMin;
float etaMax;
float ptMin;
float ptMax;
float discrMin;
float discrMax;
// default constructor
Parameters(
OperatingPoint op=OP_TIGHT,
std::string measurement_type="comb",
std::string sys_type="central",
JetFlavor jf=FLAV_B,
float eta_min=-99999.,
float eta_max=99999.,
float pt_min=0.,
float pt_max=99999.,
float discr_min=0.,
float discr_max=99999.
);
};
//BTagEntryN() {}
BTagEntryN() ;
BTagEntryN(const std::string &csvLine);
BTagEntryN(const std::string &func, Parameters p);
BTagEntryN(const TF1* func, Parameters p);
BTagEntryN(const TH1* histo, Parameters p);
~BTagEntryN() {}
static std::string makeCSVHeader();
std::string makeCSVLine() const;
static std::string trimStr(std::string str);
// public, no getters needed
std::string formula;
Parameters params;
};
#endif // BTagEntryN_H
#ifndef BTagCalibrationN_H
#define BTagCalibrationN_H
/**
* BTagCalibrationN
*
* The 'hierarchy' of stored information is this:
* - by tagger (BTagCalibrationN)
* - by operating point or reshape bin
* - by jet parton flavor
* - by type of measurement
* - by systematic
* - by eta bin
* - as 1D-function dependent of pt or discriminant
*
************************************************************/
#include <map>
#include <vector>
#include <string>
#include <istream>
#include <ostream>
class BTagCalibrationN
{
public:
//BTagCalibrationN() {}
BTagCalibrationN() ; 
BTagCalibrationN(const std::string &tagger);
BTagCalibrationN(const std::string &tagger, const std::string &filename);
~BTagCalibrationN() {}
std::string tagger() const {return tagger_;}
void addEntry(const BTagEntryN &entry);
const std::vector<BTagEntryN>& getEntries(const BTagEntryN::Parameters &par) const;
void readCSV(std::istream &s);
void readCSV(const std::string &s);
void makeCSV(std::ostream &s) const;
std::string makeCSV() const;
protected:
static std::string token(const BTagEntryN::Parameters &par);
std::string tagger_;
std::map<std::string, std::vector<BTagEntryN> > data_;
};
#endif // BTagCalibrationN_H
#ifndef BTagCalibrationNReader_H
#define BTagCalibrationNReader_H
/**
* BTagCalibrationNReader
*
* Helper class to pull out a specific set of BTagEntryN's out of a
* BTagCalibrationN. TF1 functions are set up at initialization time.
*
************************************************************/
#include <map>
#include <string>
#include <vector>
#include <TF1.h>
class BTagCalibrationNReader
{
public:
BTagCalibrationNReader() {}
BTagCalibrationNReader(const BTagCalibrationN* c,
BTagEntryN::OperatingPoint op,
std::string measurementType="comb",
std::string sysType="central");
~BTagCalibrationNReader() {}
double eval(BTagEntryN::JetFlavor jf,
float eta,
float pt,
float discr=0.) const;
struct TmpEntry {
  float etaMin;
  float etaMax;
  float ptMin;
  float ptMax;
  float discrMin;
  float discrMax;
  TF1 func;
};
protected:
void setupTmpData(const BTagCalibrationN* c);
BTagEntryN::Parameters params;
std::map<BTagEntryN::JetFlavor, std::vector<TmpEntry> > tmpData_;
std::vector<bool> useAbsEta;
};
#endif // BTagCalibrationNReader_H

