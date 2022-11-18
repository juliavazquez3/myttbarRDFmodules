###################################################
################   DEFINITIONS   ##################
###################################################

print('SV CHANNEL')

import ROOT, os, sys
from ROOT import *

### LEPTONS

## Funciones para seleccionar muones y electrones
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (pt[i]>cutpt && fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i]){
                        vb.push_back(i);
                }
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                if (pt[i]>cutpt && fabs(eta[i])<2.5 && mva80[i]){
                        vb.push_back(i);
                }
            }
            return vb;
      };
""")

## Numero de muones dentro de un jet
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muoninjet(UInt_t nmu, Vint mu_jetid, Vint mu_good) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (mu_good.size()>0){
                        if (i!=mu_good[0] && mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
                } else {
                        if (mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
                }

            }
            return vb;
      };
""")

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi, Vint puID, Vint jetID) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = false;
            for (unsigned int i=0; i<njet; ++i) {
                if(mu_good.size()>0){
                        pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
                if(el_good.size()>0){
                        pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };

      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
      auto InvariantM3(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1, const float pt2, const float eta2, const float phi2, const float mass2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);
            auto e2 = std::sqrt(x2*x2+y2*y2+z2*z2+mass2*mass2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

######### bottom jets ordering

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      template <typename T>
      vector<size_t> sort_indexes(const vector<T> &v) {

          // initialize original index locations
          vector<size_t> idx(v.size());
          iota(idx.begin(), idx.end(), 0);

          // sort indexes based on comparing values in v
          // using std::stable_sort instead of std::sort
          // to avoid unnecessary index re-orderings
          // when v contains elements of equal values
          stable_sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

          return idx;
      };
      auto bottomjets(UInt_t njetgood, Vint jetgood, UInt_t njet, Vfloat jet_btag) {
          vector<float> vb;
          vector<int> fin;
          if (njetgood>0) {
                for (unsigned int i=0; i<njetgood; ++i){
                    vb.push_back(jet_btag[jetgood[i]]);
                }
          }
          for (auto i: sort_indexes(vb)){
                fin.push_back(i);
          }
          return fin;
      };
""")

## Muon dentro del jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetMuonIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = good[0];
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = good[j];
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = i;
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = i;
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
""")

######### pedimos algunas condiciones al muon seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonCond( Vint mu_jet,Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_iso, Vbool mu_softid) {
            bool cond = false;
            if (mu_jet.size()>0){
                if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4 && mu_softid[mu_jet[0]]) {
                //if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4) {
                        cond = true;
                }
            }
            return cond;
      };
""")

## Secondary vertex

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetSVIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nSV; ++i){
                                cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[good[j]],sv_phi[i],phi[good[j]]) < 0.4;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                if(cond && pt[good[j]] > ptJ && condb){
                                        ind = good[j];
                                        ptJ = pt[good[j]];
                                }
                        }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
      auto JetSVIndSV(UInt_t njet, Vint good, Vint jet_sv, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  if (jet_sv.size() > 0) cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jet_sv[0]],sv_phi[i],phi[jet_sv[0]]) < 0.4;
                  if(cond && sv_pt[i] > ptM){
                               ind = i;
                               ptM = sv_pt[i];
                  }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
""")

######### pedimos algunas condiciones al SV seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto svCond( Vint sv_jet,Vfloat sv_pt, Vfloat sv_eta, Vint sv_charge) {
            bool cond = false;
            if (sv_jet.size()>0){
                if (sv_charge[sv_jet[0]] != 0) {
                        cond = true;
                }
            }
            return cond;
      };
""")

####### Cantidad de svs total

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto nSVtotal(UInt_t njet, Vint jetbot, Vint jetgood, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint sv_charge) {
            vector<int> vb;
            bool condb = false;
            bool cond = false;
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  ind = -1;
                  for (unsigned int j=0; j<jetgood.size(); ++j){
                             if (jetbot.size() > 1) condb = (jetgood[j] != jetgood[jetbot[0]] && jetgood[j] != jetgood[jetbot[1]]);
                             cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jetgood[j]],sv_phi[i],phi[jetgood[j]]) < 0.4;
                             if(cond && sv_charge[i] != 0 && condb){
                                     ind = i;
                             }

                  }
                  if(ind>0){
                             vb.push_back(ind);
                  }

            }
            return vb;
      };
""")

####### Cantidad de svs apropiados

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetSVMult(UInt_t njet, Vint jetsv, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint sv_charge, int n) {
            vector<int> vb;
            bool cond = false;
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  if (jetsv.size() > 0) cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jetsv[n]],sv_phi[i],phi[jetsv[n]]) < 0.4;
                  if(cond && sv_charge[i] != 0){
                          vb.push_back(i);
                  }
            }
            return vb;
      };
""")

######### Segundo jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto secondjet( UInt_t njetgood,Vint good, Vint jetmuon, Vint jetbotind, Vfloat jet_pt) {
            int ind = -1;
            float pT = -10.;
            for (unsigned int i=0; i<njetgood; ++i){
              if (jet_pt[good[i]] > pT && good[i]!=jetmuon[0] && good[i]!=good[jetbotind[0]] && good[i]!=good[jetbotind[1]]) {
                    ind = good[i];
                    pT = jet_pt[good[i]];
              }
            }
            return ind;
      };
""")

####   SSOS   ####

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto SSOS(Vint mu_good, Vint el_good, Vint sv_jet,Vint el_charge, Vint mu_charge, Vint sv_charge) {
            int ind = 0;
            if(sv_jet.size()>0){
                if(mu_good.size()>0){
                        ind = mu_charge[mu_good[0]]*sv_charge[sv_jet[0]];
                }
                if(el_good.size()>0){
                        ind = el_charge[el_good[0]]*sv_charge[sv_jet[0]];
                }
            }
            return ind;
      };
""")

## pT component calculations
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto ptot = plep1+plep2+phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py(); auto ptotX = ptot.Px(); auto ptotY = ptot.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = std::sqrt(ptotX*ptotX + ptotY*ptotY);
            return ptotmod/(plepmod+phadmod);
      };
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector plep1;
            float deltaRlepjet1;
            float deltaRlepjet2;
            float deltaPhilepjet1;
            float deltaPhilepjet2;
            float deltaEtalepjet1;
            float deltaEtalepjet2;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[0]],mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(el_phi[el_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(el_phi[el_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(fabs(plep2.Phi()-phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(fabs(plep.Phi()-phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(fabs(plep1.Phi()-phad.Phi()));
            vb.push_back(fabs(plep.Eta()-phad.Eta()));
            vb.push_back(fabs(plep1.Eta()-phad.Eta()));
            vb.push_back(fabs(plep.Pt()-phad.Pt()));
            vb.push_back(fabs(plep1.Pt()-phad.Pt()));
            vb.push_back(phad2.Pt()*std::sin(phad1.Phi()-phad2.Phi()));
            vb.push_back(phad.Pt()/(phad1.Pt()+phad2.Pt()));
            vb.push_back(plep.Pt());
            vb.push_back(deltaRlepjet1);
            vb.push_back(deltaRlepjet2);
            vb.push_back(deltaPhilepjet1);
            vb.push_back(deltaPhilepjet2);
            vb.push_back(deltaEtalepjet1);
            vb.push_back(deltaEtalepjet2);
            return vb;
      };
""")

## Trigger function for 2017 electron data

gInterpreter.Declare("""
      #include <iomanip>
      #include <math.h>
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto calDeltaR(float eta1, float phi1, float eta2, float phi2) {
          float dphi = phi1-phi2;
          float deta = eta1-eta2;
          if (dphi <= - M_PI) {
             dphi = dphi + M_PI;
          }
          if (dphi >= M_PI) {
             dphi = dphi - M_PI;
          }
          float prod = deta*deta + dphi*dphi;
          return prod;
      };
      auto triggeremulator(UInt_t nel, UInt_t ntrig, Vfloat el_eta, Vfloat el_deltaeta, Vfloat el_phi, Vfloat trig_eta, Vfloat trig_phi, Vint trig_bits, Vint trig_id) {
          bool cond = false;
          bool cond1 = false;
          bool cond2 = false;
          for (unsigned int j=0; j<nel; ++j) {
            for (unsigned int i=0; i<ntrig; ++i) {
                cond1 = calDeltaR(el_eta[j]+el_deltaeta[j],el_phi[j], trig_eta[i], trig_phi[i]) < 0.01;
                cond2 = trig_bits[i] & (0x1 << 10);
                if( trig_id[i]==11 && cond1 && cond2) {
                    cond = true;
                }
            }
          }
          return cond;
      };
""")

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=30
muon_pt["2016B"]=30
muon_pt["2017"]=30
muon_pt["2018"]=30

el_pt = {}

el_pt["2016"]=35
el_pt["2016B"]=35
el_pt["2017"]=35
el_pt["2018"]=35

## Triggers per year

muon_trig = {}

muon_trig["2016"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2016B"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2017"]="HLT_IsoMu27"
muon_trig["2018"]="HLT_IsoMu24"

el_trig = {}

el_trig["2016"]="HLT_Ele27_WPTight_Gsf"
el_trig["2016B"]="HLT_Ele27_WPTight_Gsf"
el_trig["2017"]="HLT_Ele32_WPTight_Gsf"
el_trig["2018"]="HLT_Ele32_WPTight_Gsf"

### MET Filters per year

met_filter = {}

met_filter["2016"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2016B"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2017"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_BadPFMuonFilter")
met_filter["2018"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")

#######################    end of C++ code    #########################

class selectionSV():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")
        self.isMU = kwargs.pop("isMU")

    def run(self, df):

        if not self.isMC:
           df = df.Define('Jet_pt_nom','Jet_pt')
           df = df.Define('Jet_mass_nom','Jet_mass')
           df = df.Define('MET_smeared_phi','MET_phi')
           df = df.Define('MET_smeared_pt','MET_pt')

        df = df.Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[str(self.year)])+')')
        df = df.Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[str(self.year)])+')')
        df = df.Define('nMuonGood','MuonGoodInd.size()')
        df = df.Define('nElectronGood','ElectronGoodInd.size()')
        df = df.Define('JetGoodInd','JetInd(nJet, Jet_pt_nom, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
        df = df.Define('nJetGood','JetGoodInd.size()')
        df = df.Define('JetBotInd','bottomjets(nJetGood, JetGoodInd, nJet, Jet_btagDeepFlavB)')
        df = df.Define('JetSVInd','JetSVIndJet(nJet, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, JetBotInd)')
        df = df.Define('SVJetInd','JetSVIndSV(nJet, JetGoodInd, JetSVInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta)')
        df = df.Define('SVJetGood','svCond(SVJetInd, SV_pt, SV_eta, SV_charge)')
        df = df.Define('SVLepSign','SSOS(MuonGoodInd, ElectronGoodInd, SVJetInd, Electron_charge, Muon_charge, SV_charge)')
        df = df.Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx, JetBotInd)')
        df = df.Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx, JetBotInd)')
        df = df.Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all,Muon_tightId)')
        df = df.Define('weightSSOS','-1*SVLepSign/std::abs(SVLepSign)')

        ########### Filtering and further definitions
        df = df.Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood >= 4')
        #df = df.Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood >= 4').Filter('SVJetInd.size() >= 1').Filter('SVJetGood')
        #df = df.Filter('!(MuonJetInd.size() >= 1 && MuonJetGood)')
        ## Differentiated definitions
        df = df.Filter(met_filter[str(self.year)])
        df = df.Define('lepton_type','nMuonGood>0 ? 0 : 1')
        if self.isMC:
                df = df.Filter('nMuonGood>0 ? '+ str(muon_trig[str(self.year)]) +' : '+ str(el_trig[str(self.year)]))
        else:
                if self.isMU:
                        df = df.Filter(muon_trig[str(self.year)])
                else:
                        if str(self.year) == "2017":
                                df = df.Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
                        else:
                                df = df.Filter(el_trig[str(self.year)])
        df = df.Define('lepton_pt','nMuonGood>0 ? Muon_pt[MuonGoodInd[0]] : Electron_pt[ElectronGoodInd[0]]')
        df = df.Define('lepton_eta','nMuonGood>0 ? Muon_eta[MuonGoodInd[0]] : Electron_eta[ElectronGoodInd[0]]')
        df = df.Define('transverse_mass', 'nMuonGood>0 ? std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_smeared_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_smeared_phi))): std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_smeared_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_smeared_phi)))')
        df = df.Define('lepton_iso', 'nMuonGood>0 ? Muon_pfRelIso04_all[MuonGoodInd[0]] : Electron_pfRelIso03_all[ElectronGoodInd[0]]')
        ## New cuts
        #df = df.Filter('transverse_mass > 50')
        #df = df.Filter('sv_jet_ntracks > 2')
        # BWP: cutting on the highest btag jets
        #df = df.Filter('jet_bot1_btag >'+str(cuts_btag[str(self.year)][1]))
        #df = df.Filter('jet_bot2_btag >'+str(cuts_btag[str(self.year)][0]))

        variables = ['MuonGoodInd','ElectronGoodInd','nMuonGood','nElectronGood','JetGoodInd','nJetGood','JetBotInd','JetSVInd','SVJetInd','SVJetGood','SVLepSign','weightSSOS',
           'MuonJetInd','JetMuonInd','MuonJetGood','lepton_type','lepton_pt','lepton_eta','transverse_mass','lepton_iso'        
        ]

        branches = variables

        return df,branches

def selectionSVRDF(**kwargs):
    """
    Module to compute jet pt and mass after applying smearing factors

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: jetVarRDF
            path: Base.Modules.smearing
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                isUL: self.dataset.has_tag('ul')
                ¿¿ proc: self.dataset.process ??
    """
    return lambda: selectionSV(**kwargs)


