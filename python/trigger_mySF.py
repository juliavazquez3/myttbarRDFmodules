import os

from analysis_tools.utils import import_root
ROOT = import_root()

import correctionlib
correctionlib.register_pyroot_binding()

json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{}/btagging.json.gz"

class btag_mySFRDFProducer():
    def __init__(self, *args, **kwargs):
        self.year = kwargs.pop("year")
        self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if not self.isUL:
                raise ValueError("Only implemented for UL datasets")
            if self.year == 2018:
                self.tag = "2018_UL"
                filename = json_path.format(self.tag)
            else:
                raise ValueError("2016 and 2017 not yet implemented")

            if "/libCorrectionsWrapper.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gInterpreter.Load("libCorrectionsWrapper.so")
            ROOT.gInterpreter.Declare(os.path.expandvars(
                '#include "$CMSSW_BASE/src/Corrections/Wrapper/interface/custom_sf.h"'))

            # Declaring two objects, one for tightID and the other for tightRelIso.
            # Of course this should be changed if more WPs are needed, maybe looping over 
            # an input parameter.
            ROOT.gInterpreter.ProcessLine(
                'auto corr_combbtag = '
                    'MyCorrections("%s", "deepJet_comb");'
                % filename)
            ROOT.gInterpreter.ProcessLine(
                'auto corr_inclbtag = '
                    'MyCorrections("%s", "deepJet_incl");'
                % filename)

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint = const ROOT::RVec<int>&;
                ROOT::RVec<double> get_btag_comb_MED_sf(std::string syst, Vfloat pt, Vint flav, Vfloat eta) {
                    ROOT::RVec<double> sf;
                    bool cond1 = false;
                    bool cond2 = false;
                    std::string letter = "M";
                    for (size_t i = 0; i < pt.size(); i++) {
                        cond1 = pt[i] > 20.;
                        cond2 = (flav[i] == 4 ||flav[i] == 5);
                        if (fabs(eta[i]) < 2.5 && cond1 && cond2) {
                             sf.push_back(corr_combbtag.eval({syst, letter, flav[i], fabs(eta[i]), pt[i]}));
                        }else{
                             sf.push_back(1.);
                        } 
                    }
                    return sf;
                }
                ROOT::RVec<double> get_btag_comb_LOO_sf(std::string syst, Vfloat pt, Vint flav, Vfloat eta) {
                    ROOT::RVec<double> sf;
                    bool cond1 = false;
                    bool cond2 = false;
                    std::string letter = "L";
                    for (size_t i = 0; i < pt.size(); i++) {
                        cond1 = pt[i] > 20.;
                        cond2 = (flav[i] == 4 ||flav[i] == 5);
                        if (fabs(eta[i]) < 2.5 && cond1 && cond2) {
                             sf.push_back(corr_combbtag.eval({syst, letter, flav[i], fabs(eta[i]), pt[i]}));
                        }else{
                             sf.push_back(1.);
                        }
                    }
                    return sf;
                }
                ROOT::RVec<double> get_btag_incl_MED_sf(std::string syst, Vfloat pt, Vint flav, Vfloat eta) {
                    ROOT::RVec<double> sf;
                    bool cond1 = false;
                    bool cond2 = false;
                    std::string letter = "M";
                    for (size_t i = 0; i < pt.size(); i++) {
                        cond1 = pt[i] > 20.;
                        cond2 = flav[i] == 0;
                        if (fabs(eta[i]) < 2.5 && cond1 && cond2) {
                             sf.push_back(corr_inclbtag.eval({syst, letter, flav[i], fabs(eta[i]), pt[i]}));
                        }else{
                             sf.push_back(1.);
                        }
                    }
                    return sf;
                }
                ROOT::RVec<double> get_btag_incl_LOO_sf(std::string syst, Vfloat pt, Vint flav, Vfloat eta) {
                    ROOT::RVec<double> sf;
                    bool cond1 = false;
                    bool cond2 = false;
                    std::string letter = "L";
                    for (size_t i = 0; i < pt.size(); i++) {
                        cond1 = pt[i] > 20.;
                        cond2 = flav[i] == 0;
                        if (fabs(eta[i]) < 2.5 && cond1 && cond2) {
                             sf.push_back(corr_inclbtag.eval({syst, letter, flav[i], fabs(eta[i]), pt[i]}));
                        }else{
                             sf.push_back(1.);
                        }
                    }
                    return sf;
                }
            """)

    def run(self, df):
        if not self.isMC:
            return df, []

        df = df.Define("btag_MED_sf" ,'get_btag_comb_MED_sf("central", Jet_pt, Jet_hadronFlavour, Jet_eta)')
        df = df.Define("btag_LOO_sf",'get_btag_comb_LOO_sf("central", Jet_pt, Jet_hadronFlavour, Jet_eta)')
        df = df.Define("btag_MED_incl_sf" ,'get_btag_incl_MED_sf("central", Jet_pt, Jet_hadronFlavour, Jet_eta)')
        df = df.Define("btag_LOO_incl_sf",'get_btag_incl_LOO_sf("central", Jet_pt, Jet_hadronFlavour, Jet_eta)')
#        branches = []
#        for syst_name, syst in [("", "central"), ("_up", "up"), ("_down", "down")]:
#                df = df.Define("btag_MED_sf%s" % (syst_name),'get_btag_comb_MED_sf("%s", Jet_pt, Jet_hadronFlavour, Jet_eta)' % (syst))
#                branches.append("btag_MED_sf%s" % (syst_name))
#                df = df.Define("btag_LOO_sf%s" % (syst_name),'get_btag_comb_LOO_sf("%s", Jet_pt, Jet_hadronFlavour, Jet_eta)' % (syst))
#                branches.append("btag_LOO_sf%s" % (syst_name))

        return df, ["btag_MED_sf" ,"btag_LOO_sf", "btag_MED_incl_sf" ,"btag_LOO_incl_sf"]


def btag_mySFRDF(**kwargs):
    return lambda: btag_mySFRDFProducer(**kwargs)

