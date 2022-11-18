###################################################
################   DEFINITIONS   ##################
###################################################

print('SL CHANNEL')

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

## Eleccion de los jets que no sean bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInds(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, Vint jetbotind) {
            vector<int> vb;
            bool condb = true;
            int ind1 = -1;
            int ind2 = -1;
            float ptJ{-10.};
            for (unsigned int j=0; j<good.size(); ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        if(pt[good[j]]>ptJ && condb){
                                ind1 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind1 > -1) {
                vb.push_back(ind1);
            }
            ptJ = -10.;
            for (unsigned int j=0; j<good.size(); ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        if(pt[good[j]]>ptJ && condb && good[j] != ind1){
                                ind2 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind2 > -1) {
                vb.push_back(ind2);
            }
            return vb;
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
                  if(ind>-1){
                             vb.push_back(ind);
                  }

            }
            return vb;
      };
""")

######### Masa invariante con los bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMassBot(Vint jetgood, Vint jetbot, Vint jetmuon, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, UInt_t jet_notmuon, Vfloat Jet_mass) {
            int ind = -1;
            vector<float> vb;
            float dR1 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }
            return vb;
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

class selectionWQQ():
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
        df = df.Define('JetQInd','JetInds(nJet, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, JetBotInd)')
        df = df.Define('weightSSOS','1')

        ########### Filtering and further definitions
        df = df.Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>=4')
        df = df.Define('JetnotMuonInd','nJetGood>=4 ? JetQInd[1] : -1')
        ### hists definitions
        df = df.Define('jet_1_pt','Jet_pt_nom[JetQInd[0]]')
        df = df.Define('jet_1_nmu','Jet_nMuons[JetQInd[0]]')
        df = df.Define('jet_2_pt','nJetGood>1 ? Jet_pt_nom[JetnotMuonInd] : 0')
        df = df.Define('jet_1_eta','Jet_eta[JetQInd[0]]')
        df = df.Define('jet_2_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
        df = df.Define('jet_2_mass','nJetGood>1 ? Jet_mass_nom[JetnotMuonInd] : 0')
        df = df.Define('jet_2_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
        df = df.Define('jet_2_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
        df = df.Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],Jet_mass_nom[JetQInd[0]],Jet_pt_nom[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass_nom[JetnotMuonInd]) : 0')
        df = df.Define('deltaR_jet1_jet2','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetQInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetQInd[0]])  : 10')
        df = df.Define('deltaphi_jet1_jet2','fabs(Jet_phi[JetQInd[0]]-Jet_phi[JetnotMuonInd])')
        df = df.Define('deltaeta_jet1_jet2','fabs(Jet_eta[JetQInd[0]]-Jet_eta[JetnotMuonInd])')
        df = df.Define('deltapt_jet1_jet2','fabs(Jet_pt_nom[JetQInd[0]]-Jet_pt_nom[JetnotMuonInd])')
        df = df.Define('tracks_jet1','Jet_nConstituents[JetQInd[0]]')
        df = df.Define('tracks_jet2','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
        df = df.Define('EMN_jet1','Jet_neEmEF[JetQInd[0]]')
        df = df.Define('EMC_jet1','Jet_chEmEF[JetQInd[0]]')
        df = df.Define('EMtotal_jet1','Jet_chEmEF[JetQInd[0]]+Jet_neEmEF[JetQInd[0]]')
        df = df.Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
        df = df.Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
        df = df.Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
        df = df.Define('deltaR_lep_2jets','aux_various[0]')
        df = df.Define('deltaphi_MET_2jets','aux_various[1]')
        df = df.Define('eta_2jets','aux_various[2]')
        df = df.Define('pt_2jets','aux_various[3]')
        df = df.Define('deltaphi_lephad','aux_various[4]')
        df = df.Define('deltaR_lephad','aux_various[5]')
        df = df.Define('deltaphi_lep_2jets','aux_various[6]')
        df = df.Define('deltaeta_lephad','aux_various[7]')
        df = df.Define('deltaeta_lep_2jets','aux_various[8]')
        df = df.Define('deltapt_lephad','aux_various[9]')
        df = df.Define('deltapt_lep_2jets','aux_various[10]')
        df = df.Define('pT_proy','aux_various[11]')
        df = df.Define('pT_sum_2J','aux_various[12]')
        df = df.Define('pT_Wlep','aux_various[13]')
        df = df.Define('deltaR_lep_jet1','aux_various[14]')
        df = df.Define('deltaR_lep_jet2','aux_various[15]')
        df = df.Define('deltaPhi_lep_jet1','aux_various[16]')
        df = df.Define('deltaPhi_lep_jet2','aux_various[17]')
        df = df.Define('deltaEta_lep_jet1','aux_various[18]')
        df = df.Define('deltaEta_lep_jet2','aux_various[19]')
        df = df.Define('jet_bot1_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]]')
        df = df.Define('jet_bot2_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]]')
        df = df.Define('jet_1_btag','Jet_btagDeepFlavB[JetQInd[0]]')
        df = df.Define('jet_2_btag','Jet_btagDeepFlavB[JetnotMuonInd]')
        df = df.Define('jet_2_deeptagG','Jet_btagDeepFlavQG[JetnotMuonInd]')
        df = df.Define('jet_2_deeptagC','Jet_btagDeepFlavCvL[JetnotMuonInd]')
        df = df.Define('nLooseLepton','nMuon+nElectron-1')
        df = df.Define('nSV_nobot_aux','nSVtotal(nJet, JetBotInd, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, SV_charge)')
        df = df.Define('nSV_nobot','nSV_nobot_aux.size()')
        df = df.Define('InvM3_aux','InvMassBot(JetGoodInd, JetBotInd, JetQInd, Jet_pt_nom, Jet_eta, Jet_phi, JetnotMuonInd, Jet_mass_nom)')
        df = df.Define('InvM_bot_closer','InvM3_aux[0]')
        df = df.Define('InvM_bot_farther','InvM3_aux[1]')


        ############ Gen level definitions
        if self.isMC:
                df = df.Define('jet_1_flavourH','Jet_hadronFlavour[JetQInd[0]]')
                df = df.Define('jet_2_flavourH','Jet_hadronFlavour[JetnotMuonInd]')
                df = df.Define('jet_1_flavourP','Jet_partonFlavour[JetQInd[0]]')
                df = df.Define('jet_2_flavourP','Jet_partonFlavour[JetnotMuonInd]') 
        else:
                df = df.Define('jet_muon_flavourH','-1')
                df = df.Define('jet_notmuon_flavourH','-1')
                df = df.Define('jet_muon_flavourP','-1')
                df = df.Define('jet_notmuon_flavourP','-1')

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
        df = df.Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetQInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetQInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetQInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetQInd[0]])')
        df = df.Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],Jet_mass_nom[JetQInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_nom[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],Jet_mass_nom[JetQInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df = df.Define('transverse_mass', 'nMuonGood>0 ? std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_smeared_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_smeared_phi))): std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_smeared_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_smeared_phi)))')
        df = df.Define('InvM_muon_jet','nMuonGood>0 ? InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : 50.')
        df = df.Define('lepton_iso', 'nMuonGood>0 ? Muon_pfRelIso04_all[MuonGoodInd[0]] : Electron_pfRelIso03_all[ElectronGoodInd[0]]')
        ## New cuts
        df = df.Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
        df = df.Filter('transverse_mass > 50')
        df = df.Filter('muon_jet_relpt<0.5')
        df = df.Filter('EMtotal_jetM<0.6')
        ## BWP: cutting on btagging working points
        df = df.Filter('jet_bot1_btag >'+str(cuts_btag[str(self.year)][1]))
        df = df.Filter('jet_bot2_btag >'+str(cuts_btag[str(self.year)][0]))


#        variables = ['MuonGoodInd','ElectronGoodInd','nMuonGood','nElectronGood','JetGoodInd','nJetGood','JetBotInd','MuonJetInd','JetMuonInd','MuonJetGood','MuonLepSign','weightSSOS']

        variables = ['MuonGoodInd','ElectronGoodInd','nMuonGood','nElectronGood','JetGoodInd','nJetGood','JetBotInd','JetQInd','weightSSOS',
           'JetnotMuonInd','jet_1_pt','jet_1_nmu','jet_1_mass','jet_2_pt','jet_1_eta','jet_2_eta','jet_2_mass','jet_2_qgl','jet_2_nmu',
           'InvM_2jets','deltaR_jet1_jet2','deltaphi_jet1_jet2','deltaeta_jet1_jet2','deltapt_jet1_jet2','tracks_jet1','tracks_jet2','EMN_jet1','EMC_jet1','EMtotal_jet1','pT_sum','pT_product',
           'aux_various','deltaR_lep_2jets','deltaphi_MET_2jets','eta_2jets','pt_2jets','deltaphi_lephad','deltaR_lephad','deltaphi_lep_2jets','deltaeta_lephad','deltaeta_lep_2jets','deltapt_lephad',
           'deltapt_lep_2jets','pT_proy','pT_sum_2J','pT_Wlep','deltaR_lep_jet1','deltaR_lep_jet2','deltaPhi_lep_jet1','deltaPhi_lep_jet2','deltaEta_lep_jet1','deltaEta_lep_jet2',
           'jet_bot1_btag','jet_bot2_btag','jet_1_btag','jet_2_btag','jet_2_deeptagG','jet_2_deeptagC','nLooseLepton','jet_1_flavourH','jet_2_flavourH','jet_1_flavourP',
           'jet_2_flavourP','lepton_type','lepton_pt','lepton_eta','transverse_mass','lepton_iso','nSV_nobot','InvM_bot_closer','InvM_bot_farther'
        ]

        branches = variables

        return df,branches

def selectionWQQRDF(**kwargs):
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
    return lambda: selectionWQQ(**kwargs)


