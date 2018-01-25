#include "BTagSFUtil.C"

void prueba(){
 BTagSFUtil *fBTagSFnom = new BTagSFUtil("comb", "./", "CSVv2", "Medium",  0);
  vector<string> sysTypes;
  sysTypes.push_back("up");
  sysTypes.push_back("down");


  const BTagCalibration calib("CSVv2", "CSVv2.csv");
  BTagCalibrationReader *FastSimReader_b = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", sysTypes);
  FastSimReader_b -> load(calib, BTagEntry::FLAV_B, "fastsim");
  BTagCalibrationReader *FastSimReader_c = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", sysTypes);
  FastSimReader_c -> load(calib, BTagEntry::FLAV_C, "fastsim");
  BTagCalibrationReader *FastSimReader_c = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", sysTypes);
  FastSimReader_l -> load(calib, BTagEntry::FLAV_UDSG, "fastsim");



}

Float_t 
  return FSSF;
}
