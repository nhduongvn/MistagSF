#include "BTagCalibrationStandalone.h"
#include <iostream>
#include <exception>
#include <algorithm>
#include <sstream>
BTagEntryN::Parameters::Parameters(
OperatingPoint op,
//std::string op,
std::string measurement_type,
std::string sys_type,
JetFlavor jf,
float eta_min,
float eta_max,
float pt_min,
float pt_max,
float discr_min,
float discr_max
):
operatingPoint(op),
measurementType(measurement_type),
sysType(sys_type),
jetFlavor(jf),
etaMin(eta_min),
etaMax(eta_max),
ptMin(pt_min),
ptMax(pt_max),
discrMin(discr_min),
discrMax(discr_max)
{}
BTagEntryN::BTagEntryN(const std::string &csvLine)
{
// make tokens
std::stringstream buff(csvLine);
std::vector<std::string> vec;
std::string token;
while (std::getline(buff, token, ","[0])) {
token = BTagEntryN::trimStr(token);
if (token.empty()) {
continue;
}
vec.push_back(token);
}
if (vec.size() != 11) {
std::cerr << "ERROR in BTagCalibrationN: "
<< "Invalid csv line; num tokens != 11: "
<< csvLine;
throw std::exception();
}
// clean string values
char chars[] = " \"\n";
for (unsigned int i = 0; i < strlen(chars); ++i) {
vec[1].erase(remove(vec[1].begin(),vec[1].end(),chars[i]),vec[1].end());
vec[2].erase(remove(vec[2].begin(),vec[2].end(),chars[i]),vec[2].end());
vec[10].erase(remove(vec[10].begin(),vec[10].end(),chars[i]),vec[10].end());
}
// make formula
formula = vec[10];
TF1 f1("", formula.c_str()); // compile formula to check validity
if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibrationN: "
<< "Invalid csv line; formula does not compile: "
<< csvLine;
throw std::exception();
}
// make parameters
if (stoi(vec[0]) > 3) {
std::cerr << "ERROR in BTagCalibrationN: "
<< "Invalid csv line; OperatingPoint > 3: "
<< csvLine;
throw std::exception();
}
if (stoi(vec[3]) > 2) {
std::cerr << "ERROR in BTagCalibrationN: "
<< "Invalid csv line; JetFlavor > 2: "
<< csvLine;
throw std::exception();
}
params = BTagEntryN::Parameters(
BTagEntryN::OperatingPoint(stoi(vec[0])),
vec[1],
vec[2],
BTagEntryN::JetFlavor(stoi(vec[3])),
stof(vec[4]),
stof(vec[5]),
stof(vec[6]),
stof(vec[7]),
stof(vec[8]),
stof(vec[9])
);
}
BTagEntryN::BTagEntryN(const std::string &func, BTagEntryN::Parameters p):
formula(func),
params(p)
{}
BTagEntryN::BTagEntryN(const TF1* func, BTagEntryN::Parameters p):
formula(std::string(func->GetExpFormula("p").Data())),
params(p)
{}
// Creates chained step functions like this:
// "<prevous_bin> : x<bin_high_bound ? bin_value : <next_bin>"
// e.g. "x<0 ? 1 : x<1 ? 2 : x<2 ? 3 : 4"
BTagEntryN::BTagEntryN(const TH1* hist, BTagEntryN::Parameters p):
params(p)
{
int nbins = hist->GetNbinsX();
auto axis = hist->GetXaxis();
// overwrite bounds with histo values
if (params.operatingPoint == BTagEntryN::OP_RESHAPING) {
params.discrMin = axis->GetBinLowEdge(1);
params.discrMax = axis->GetBinUpEdge(nbins);
} else {
params.ptMin = axis->GetBinLowEdge(1);
params.ptMax = axis->GetBinUpEdge(nbins);
}
std::stringstream buff;
buff << "x<" << axis->GetBinLowEdge(1) << " ? 0. : "; // default value
for (int i=1; i<nbins+1; ++i) {
char tmp_buff[100];
sprintf(tmp_buff,
"x<%g ? %g : ", // %g is the smaller one of %e or %f
axis->GetBinUpEdge(i),
hist->GetBinContent(i));
buff << tmp_buff;
}
buff << 0.; // default value
formula = buff.str();
}
std::string BTagEntryN::makeCSVHeader()
{
return "wp,"
"type,"
"syst,"
"flav,"
"etaMin,"
"etaMax,"
"ptMin,"
"ptMax,"
"discrMin,"
"discrMax,"
"formula\n";
}
std::string BTagEntryN::makeCSVLine() const
{
std::stringstream buff;
std::string wp;
if (params.operatingPoint == BTagEntryN::OP_LOOSE) wp = "L";  
if (params.operatingPoint == BTagEntryN::OP_MEDIUM) wp = "M";  
if (params.operatingPoint == BTagEntryN::OP_TIGHT) wp = "T";  
if (params.operatingPoint == BTagEntryN::OP_XTIGHT) wp = "XT";  
if (params.operatingPoint == BTagEntryN::OP_XXTIGHT) wp = "XXT";  
buff << wp
<< "," << params.measurementType
<< "," << params.sysType
<< "," << params.jetFlavor
<< "," << params.etaMin
<< "," << params.etaMax
<< "," << params.ptMin
<< "," << params.ptMax
<< "," << params.discrMin
<< "," << params.discrMax
<< ",\"" << formula
<< "\"\n";
return buff.str();
}
std::string BTagEntryN::trimStr(std::string str) {
size_t s = str.find_first_not_of(" \n\r\t");
size_t e = str.find_last_not_of (" \n\r\t");
if((std::string::npos == s) || (std::string::npos == e))
return "";
else
return str.substr(s, e-s+1);
}
#include <fstream>
#include <sstream>
BTagCalibrationN::BTagCalibrationN(const std::string &taggr):
tagger_(taggr)
{}
BTagCalibrationN::BTagCalibrationN(const std::string &taggr,
const std::string &filename):
tagger_(taggr)
{
std::ifstream ifs(filename);
readCSV(ifs);
ifs.close();
}
void BTagCalibrationN::addEntry(const BTagEntryN &entry)
{
data_[token(entry.params)].push_back(entry);
}
const std::vector<BTagEntryN>& BTagCalibrationN::getEntries(
const BTagEntryN::Parameters &par) const
{
auto tok = token(par);
if (!data_.count(tok)) {
std::cerr << "ERROR in BTagCalibrationN: "
<< "(OperatingPoint, measurementType, sysType) not available: "
<< tok;
throw std::exception();
}
return data_.at(tok);
}
void BTagCalibrationN::readCSV(const std::string &s)
{
std::stringstream buff(s);
readCSV(buff);
}
void BTagCalibrationN::readCSV(std::istream &s)
{
std::string line;
// firstline might be the header
getline(s,line);
if (line.find("OperatingPoint") == std::string::npos) {
addEntry(BTagEntryN(line));
}
while (getline(s,line)) {
line = BTagEntryN::trimStr(line);
if (line.empty()) { // skip empty lines
continue;
}
addEntry(BTagEntryN(line));
}
}
void BTagCalibrationN::makeCSV(std::ostream &s) const
{
s << BTagEntryN::makeCSVHeader();
for (auto i = data_.cbegin(); i != data_.cend(); ++i) {
auto vec = i->second;
for (auto j = vec.cbegin(); j != vec.cend(); ++j) {
s << j->makeCSVLine();
}
}
}
std::string BTagCalibrationN::makeCSV() const
{
std::stringstream buff;
makeCSV(buff);
return buff.str();
}
std::string BTagCalibrationN::token(const BTagEntryN::Parameters &par)
{
std::stringstream buff;
buff << par.operatingPoint << ", "
<< par.measurementType << ", "
<< par.sysType;
return buff.str();
}
BTagCalibrationNReader::BTagCalibrationNReader(const BTagCalibrationN* c,
BTagEntryN::OperatingPoint op,
std::string measurementType,
std::string sysType):
params(BTagEntryN::Parameters(op, measurementType, sysType)),
useAbsEta(true)
{
setupTmpData(c);
}
double BTagCalibrationNReader::eval(BTagEntryN::JetFlavor jf,
float eta,
float pt,
float discr) const
{
bool use_discr = (params.operatingPoint == BTagEntryN::OP_RESHAPING);
if (useAbsEta[jf] && eta < 0) {
eta = -eta;
}
// search linearly through eta, pt and discr ranges and eval
// future: find some clever data structure based on intervals
const auto &entries = tmpData_.at(jf);
for (unsigned i=0; i<entries.size(); ++i) {
const BTagCalibrationNReader::TmpEntry &e = entries.at(i);
if (
e.etaMin <= eta && eta < e.etaMax // find eta
&& e.ptMin <= pt && pt < e.ptMax // check pt
){
if (use_discr) { // discr. reshaping?
if (e.discrMin <= discr && discr < e.discrMax) { // check discr
return e.func.Eval(discr);
}
} else {
return e.func.Eval(pt);
}
}
}
return 0.; // default value
}
void BTagCalibrationNReader::setupTmpData(const BTagCalibrationN* c)
{
useAbsEta = std::vector<bool>(4, true);
const auto &entries = c->getEntries(params);
for (unsigned i=0; i<entries.size(); ++i) {
const BTagEntryN &be = entries[i];
BTagCalibrationNReader::TmpEntry te;
te.etaMin = be.params.etaMin;
te.etaMax = be.params.etaMax;
te.ptMin = be.params.ptMin;
te.ptMax = be.params.ptMax;
te.discrMin = be.params.discrMin;
te.discrMax = be.params.discrMax;
if (params.operatingPoint == BTagEntryN::OP_RESHAPING) {
te.func = TF1("", be.formula.c_str(),
be.params.discrMin, be.params.discrMax);
} else {
te.func = TF1("", be.formula.c_str(),
be.params.ptMin, be.params.ptMax);
}
tmpData_[be.params.jetFlavor].push_back(te);
if (te.etaMin < 0) {
useAbsEta[be.params.jetFlavor] = false;
}
}
}

