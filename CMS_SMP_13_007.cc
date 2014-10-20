#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include <iostream>
#include <ostream>

namespace Rivet {


  class CMS_SMP_13_007 : public Analysis {
  public:

    CMS_SMP_13_007()
      : Analysis("CMS_SMP_13_007")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }

  public:
 // Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      ZFinder zfinder_zee(fs,-5.,5., 0.0*GeV, ELECTRON, 71.0*GeV, 111.0*GeV, 0.1,true,false);
      ZFinder zfinder_zmm(fs,-5.,5., 0.0*GeV, MUON, 71.0*GeV, 111.0*GeV, 0.1,true,false);

      addProjection(zfinder_zee, "ZFinder_zee");
      addProjection(zfinder_zmm, "ZFinder_zmm");


    // Define veto FS 
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(zfinder_zee);
      vfs.addVetoOnThisFinalState(zfinder_zmm);

      FastJets fastjets(vfs, FastJets::ANTIKT, 0.5);
      addProjection(fastjets, "Jets");


      vector<double> jetPt_Zinc1jet ;
      vector<double> jetPt_Zinc2jet ;
      vector<double> jetPt_Zinc3jet ;
      vector<double> jetPt_Zinc4jet ;
      vector<double> jetPt_Zinc5jet ;

      vector<double> jetHT_Zinc1jet ;
      vector<double> jetHT_Zinc2jet ;
      vector<double> jetHT_Zinc3jet ;
      vector<double> jetHT_Zinc4jet ;
      vector<double> jetHT_Zinc5jet ;

      jetPt_Zinc1jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000;
      jetPt_Zinc2jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 800;
      jetPt_Zinc3jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 700;
      jetPt_Zinc4jet += 20, 24, 30, 39, 49, 62, 78, 96, 150;
      jetPt_Zinc5jet += 20, 24, 30, 39, 49, 62, 100;

      jetHT_Zinc1jet += 30, 40, 60, 90, 120, 160, 210, 260, 330, 410, 500, 620, 750, 900, 1200;
      jetHT_Zinc2jet += 60, 90, 120, 160, 210, 260, 330, 410, 500, 620, 750, 1000;
      jetHT_Zinc3jet += 90, 120, 160, 210, 260, 330, 410, 500, 620, 800;
      jetHT_Zinc4jet += 120, 140, 167, 203, 253, 320, 410, 530, 690, 910;
      jetHT_Zinc5jet += 150, 180, 222, 282, 365, 485, 650, 880;

     _hist_ee_JetPt1j =bookHistogram1D("jet_ee_Pt1jetcase", jetPt_Zinc1jet);
     _hist_ee_JetPt2j =bookHistogram1D("jet_ee_Pt2jetcase", jetPt_Zinc2jet);
     _hist_ee_JetPt3j =bookHistogram1D("jet_ee_Pt3jetcase", jetPt_Zinc3jet);
     _hist_ee_JetPt4j =bookHistogram1D("jet_ee_Pt4jetcase", jetPt_Zinc4jet);
     _hist_ee_JetPt5j =bookHistogram1D("jet_ee_Pt5jetcase", jetPt_Zinc5jet);

     _hist_ee_Jeteta1j =bookHistogram1D("jet_ee_eta1jetcase", 32, 0, 2.4);
     _hist_ee_Jeteta2j =bookHistogram1D("jet_ee_eta2jetcase", 32, 0, 2.4);
     _hist_ee_Jeteta3j =bookHistogram1D("jet_ee_eta3jetcase", 24, 0, 2.4);
     _hist_ee_Jeteta4j =bookHistogram1D("jet_ee_eta4jetcase", 12, 0, 2.4);
     _hist_ee_Jeteta5j =bookHistogram1D("jet_ee_eta5jetcase", 6,  0, 2.4);
     _hist_ee_Jeteta6j =bookHistogram1D("jet_ee_eta6jetcase", 6,  0, 2.4);

     _hist_mm_JetPt1j =bookHistogram1D("jet_mm_Pt1jetcase", jetPt_Zinc1jet);
     _hist_mm_JetPt2j =bookHistogram1D("jet_mm_Pt2jetcase", jetPt_Zinc2jet);
     _hist_mm_JetPt3j =bookHistogram1D("jet_mm_Pt3jetcase", jetPt_Zinc3jet);
     _hist_mm_JetPt4j =bookHistogram1D("jet_mm_Pt4jetcase", jetPt_Zinc4jet);
     _hist_mm_JetPt5j =bookHistogram1D("jet_mm_Pt5jetcase", jetPt_Zinc5jet);

     _hist_mm_Jeteta1j =bookHistogram1D("jet_mm_eta1jetcase", 32, 0, 2.4);
     _hist_mm_Jeteta2j =bookHistogram1D("jet_mm_eta2jetcase", 32, 0, 2.4);
     _hist_mm_Jeteta3j =bookHistogram1D("jet_mm_eta3jetcase", 24, 0, 2.4);
     _hist_mm_Jeteta4j =bookHistogram1D("jet_mm_eta4jetcase", 12, 0, 2.4);
     _hist_mm_Jeteta5j =bookHistogram1D("jet_mm_eta5jetcase", 6, 0, 2.4);
     _hist_mm_Jeteta6j =bookHistogram1D("jet_mm_eta6jetcase", 6, 0, 2.4);

     _hist_ee_Ht_1j =bookHistogram1D("JetsHT_ee_inc1jet", jetHT_Zinc1jet);
     _hist_ee_Ht_2j =bookHistogram1D("JetsHT_ee_inc2jet", jetHT_Zinc2jet);
     _hist_ee_Ht_3j =bookHistogram1D("JetsHT_ee_inc3jet", jetHT_Zinc3jet);
     _hist_ee_Ht_4j =bookHistogram1D("JetsHT_ee_inc4jet", jetHT_Zinc4jet);
     _hist_ee_Ht_5j =bookHistogram1D("JetsHT_ee_inc5jet", jetHT_Zinc5jet);

     _hist_mm_Ht_1j =bookHistogram1D("JetsHT_mm_inc1jet", jetHT_Zinc1jet);
     _hist_mm_Ht_2j =bookHistogram1D("JetsHT_mm_inc2jet", jetHT_Zinc2jet);
     _hist_mm_Ht_3j =bookHistogram1D("JetsHT_mm_inc3jet", jetHT_Zinc3jet);
     _hist_mm_Ht_4j =bookHistogram1D("JetsHT_mm_inc4jet", jetHT_Zinc4jet);
     _hist_mm_Ht_5j =bookHistogram1D("JetsHT_mm_inc5jet", jetHT_Zinc5jet);

     _hist_inc_ZeeJetMult = bookHistogram1D("njetZee_incl", 8, -0.5, 7.5);
     _hist_excl_ZeeJetMult= bookHistogram1D("njetZee_excl", 8, -0.5, 7.5);

     _hist_inc_ZmmJetMult = bookHistogram1D("njetZmm_incl", 8, -0.5, 7.5);
     _hist_excl_ZmmJetMult= bookHistogram1D("njetZmm_excl", 8, -0.5, 7.5);


  }

  void Fill(AIDA::IHistogram1D*& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){
       _histJetMult->fill(0, weight);
       for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
         if (i==6) break;
         _histJetMult->fill(i+1, weight);  // inclusive multiplicity
      }
   }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      const ZFinder& zfinder_zee = applyProjection<ZFinder>(event, "ZFinder_zee");
      const ZFinder& zfinder_zmm = applyProjection<ZFinder>(event, "ZFinder_zmm");

      if (zfinder_zee.bosons().size()==1) {

      FourMomentum zeemom(zfinder_zee.particles()[0].momentum());

      Particle lepton0 = zfinder_zee.particles().at(0);
      Particle lepton1 = zfinder_zee.particles().at(1);

      if (lepton0.pdgId() != -lepton1.pdgId())
           vetoEvent;

      double pt0 = lepton0.momentum().pT()/GeV;
      double pt1 = lepton1.momentum().pT()/GeV;
      double eta0 = lepton0.momentum().eta();
      double eta1 = lepton1.momentum().eta();
      double phi0 = lepton0.momentum().azimuthalAngle();
      double phi1 = lepton1.momentum().azimuthalAngle();

      bool inAcceptance = pt0>20 && pt1>20 && (fabs(eta0) < 1.4442 || fabs(eta0) > 1.566) && (fabs(eta1) < 1.4442 || fabs(eta1) > 1.566);

      if (!inAcceptance)return;

      //Obtain the jets.

      vector<FourMomentum> finaljet_list;
      double HT = 0.0;

      foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
        const double jeta = j.momentum().eta();
        const double jpt = j.momentum().pT();

        if (fabs(jeta) < 2.4)
        if(jpt>30){
                finaljet_list.push_back(j.momentum());
                HT += j.momentum().pT();
           }
      }
//Multiplicity plots.  


    _hist_excl_ZeeJetMult->fill(finaljet_list.size(), weight);

    Fill(_hist_inc_ZeeJetMult, weight, finaljet_list);

    if(finaljet_list.size()>=1) {
                  _hist_ee_JetPt1j->fill(finaljet_list[0].pT(), weight);
                  _hist_ee_Jeteta1j->fill(finaljet_list[0].eta(), weight);
                  _hist_ee_Ht_1j->fill(HT, weight);

    }

    if(finaljet_list.size()>=2) {
                  _hist_ee_JetPt2j->fill(finaljet_list[1].pT(), weight);
                  _hist_ee_Jeteta2j->fill(finaljet_list[1].eta(), weight);
                  _hist_ee_Ht_2j->fill(HT, weight);

    }

    if(finaljet_list.size()>=3) {
                  _hist_ee_JetPt3j->fill(finaljet_list[2].pT(), weight);
                  _hist_ee_Jeteta3j->fill(finaljet_list[2].eta(), weight);
                  _hist_ee_Ht_3j->fill(HT, weight);

    }

    if(finaljet_list.size()>=4) {
                  _hist_ee_JetPt4j->fill(finaljet_list[3].pT(), weight);
                  _hist_ee_Jeteta4j->fill(finaljet_list[3].eta(), weight);
                  _hist_ee_Ht_4j->fill(HT, weight);

    }

    if(finaljet_list.size()>=5) {
                  _hist_ee_JetPt5j->fill(finaljet_list[4].pT(), weight);
                  _hist_ee_Jeteta5j->fill(finaljet_list[4].eta(), weight);
                  _hist_ee_Ht_5j->fill(HT, weight);

    }

    if(finaljet_list.size()>=6) {
                  _hist_ee_Jeteta6j->fill(finaljet_list[5].eta(), weight);

    }


} ////close the electron loop 

 if (zfinder_zmm.bosons().size()==1) {

      FourMomentum zmmmom(zfinder_zmm.particles()[0].momentum());

      Particle lepton0 = zfinder_zmm.particles().at(0);
      Particle lepton1 = zfinder_zmm.particles().at(1);

       if (lepton0.pdgId() != -lepton1.pdgId())
           vetoEvent;

      double pt0 = lepton0.momentum().pT()/GeV;
      double pt1 = lepton1.momentum().pT()/GeV;
      double eta0 = lepton0.momentum().eta();
      double eta1 = lepton1.momentum().eta();
      double phi0 = lepton0.momentum().azimuthalAngle();
      double phi1 = lepton1.momentum().azimuthalAngle();

     bool inAcceptance = fabs(eta0)<2.4 && fabs(eta1)<2.4 && pt0>20 && pt1>20;

     if (!inAcceptance)return;


      //Obtain the jets.

      vector<FourMomentum> finaljet_list;
      double HT = 0.0;

      foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
  
        const double jeta = j.momentum().eta();
        const double jphi = j.momentum().phi();                                                                                                                         260,1         57%
        const double jpt = j.momentum().pT();

        if (fabs(jeta) < 2.4)
        if(jpt>30){
                finaljet_list.push_back(j.momentum());
                HT += j.momentum().pT();
           }
     }


      //Multiplicity plots.     
      _hist_excl_ZmmJetMult->fill(finaljet_list.size(), weight);
      Fill(_hist_inc_ZmmJetMult, weight, finaljet_list);

      if(finaljet_list.size()>=1) {
            _hist_mm_JetPt1j->fill(finaljet_list[0].pT(), weight);
            _hist_mm_Jeteta1j->fill(finaljet_list[0].eta(), weight);
            _hist_mm_Ht_1j->fill(HT, weight);
      }

      if(finaljet_list.size()>=2) {
            _hist_mm_JetPt2j->fill(finaljet_list[1].pT(), weight);
            _hist_mm_Jeteta2j->fill(finaljet_list[1].eta(), weight);
            _hist_mm_Ht_2j->fill(HT, weight);
      }

      if(finaljet_list.size()>=3) {
           _hist_mm_JetPt3j->fill(finaljet_list[2].pT(), weight);
           _hist_mm_Jeteta3j->fill(finaljet_list[2].eta(), weight);
           _hist_mm_Ht_3j->fill(HT, weight);
      }

      if(finaljet_list.size()>=4) {
           _hist_mm_JetPt4j->fill(finaljet_list[3].pT(), weight);
           _hist_mm_Jeteta4j->fill(finaljet_list[3].eta(), weight);
           _hist_mm_Ht_4j->fill(HT, weight);
     }


      if(finaljet_list.size()>=5) {
           _hist_mm_JetPt5j->fill(finaljet_list[4].pT(), weight);
           _hist_mm_Jeteta5j->fill(finaljet_list[4].eta(), weight);
           _hist_mm_Ht_5j->fill(HT, weight);
     }

      if(finaljet_list.size()>=6) {
           _hist_mm_Jeteta6j->fill(finaljet_list[5].eta(), weight);
     }


  }// close the muon loop

}


    /// Normalise histograms etc., after the run
    void finalize() {

        double crossSec=3503.71;

        scale(_hist_inc_ZeeJetMult, crossSec/sumOfWeights());
        scale(_hist_inc_ZmmJetMult, crossSec/sumOfWeights());
        scale(_hist_excl_ZeeJetMult, crossSec/sumOfWeights());
        scale(_hist_excl_ZmmJetMult, crossSec/sumOfWeights());

        scale(_hist_ee_JetPt1j, crossSec/sumOfWeights());
        scale(_hist_ee_JetPt2j, crossSec/sumOfWeights());
        scale(_hist_ee_JetPt3j, crossSec/sumOfWeights());
        scale(_hist_ee_JetPt4j, crossSec/sumOfWeights());
        scale(_hist_ee_JetPt5j, crossSec/sumOfWeights());
      
        scale(_hist_mm_JetPt1j, crossSec/sumOfWeights());
        scale(_hist_mm_JetPt2j, crossSec/sumOfWeights());
        scale(_hist_mm_JetPt3j, crossSec/sumOfWeights());
        scale(_hist_mm_JetPt4j, crossSec/sumOfWeights());
        scale(_hist_mm_JetPt5j, crossSec/sumOfWeights());

        scale(_hist_ee_Jeteta1j, crossSec/sumOfWeights());
        scale(_hist_ee_Jeteta2j, crossSec/sumOfWeights());
        scale(_hist_ee_Jeteta3j, crossSec/sumOfWeights());
        scale(_hist_ee_Jeteta4j, crossSec/sumOfWeights());
        scale(_hist_ee_Jeteta5j, crossSec/sumOfWeights());
        scale(_hist_ee_Jeteta6j, crossSec/sumOfWeights());
      
        scale(_hist_mm_Jeteta1j, crossSec/sumOfWeights());
        scale(_hist_mm_Jeteta2j, crossSec/sumOfWeights());
        scale(_hist_mm_Jeteta3j, crossSec/sumOfWeights());
        scale(_hist_mm_Jeteta4j, crossSec/sumOfWeights());
        scale(_hist_mm_Jeteta5j, crossSec/sumOfWeights());
        scale(_hist_mm_Jeteta6j, crossSec/sumOfWeights());

        scale(_hist_ee_Ht_1j, crossSec/sumOfWeights());
        scale(_hist_ee_Ht_2j, crossSec/sumOfWeights());
        scale(_hist_ee_Ht_3j, crossSec/sumOfWeights());
        scale(_hist_ee_Ht_4j, crossSec/sumOfWeights());
        scale(_hist_ee_Ht_5j, crossSec/sumOfWeights());
      
        scale(_hist_mm_Ht_1j, crossSec/sumOfWeights());
        scale(_hist_mm_Ht_2j, crossSec/sumOfWeights());
        scale(_hist_mm_Ht_3j, crossSec/sumOfWeights());
        scale(_hist_mm_Ht_4j, crossSec/sumOfWeights());
        scale(_hist_mm_Ht_5j, crossSec/sumOfWeights());

 private:

    // Data members like post-cuts event weight counters go here

  private:

        AIDA::IHistogram1D* _hist_inc_ZeeJetMult;
        AIDA::IHistogram1D* _hist_inc_ZmmJetMult;
        AIDA::IHistogram1D* _hist_excl_ZeeJetMult;
        AIDA::IHistogram1D* _hist_excl_ZmmJetMult;

        AIDA::IHistogram1D* _hist_ee_JetPt1j;
        AIDA::IHistogram1D* _hist_ee_JetPt2j;
        AIDA::IHistogram1D* _hist_ee_JetPt3j;
        AIDA::IHistogram1D* _hist_ee_JetPt4j;
        AIDA::IHistogram1D* _hist_ee_JetPt5j;

        AIDA::IHistogram1D* _hist_mm_JetPt1j;
        AIDA::IHistogram1D* _hist_mm_JetPt2j;
        AIDA::IHistogram1D* _hist_mm_JetPt3j;
        AIDA::IHistogram1D* _hist_mm_JetPt4j;
        AIDA::IHistogram1D* _hist_mm_JetPt5j;

        AIDA::IHistogram1D* _hist_ee_Jeteta1j;
        AIDA::IHistogram1D* _hist_ee_Jeteta2j;
        AIDA::IHistogram1D* _hist_ee_Jeteta3j;
        AIDA::IHistogram1D* _hist_ee_Jeteta4j;
        AIDA::IHistogram1D* _hist_ee_Jeteta5j;
        AIDA::IHistogram1D* _hist_ee_Jeteta6j;

        AIDA::IHistogram1D* _hist_mm_Jeteta1j;
        AIDA::IHistogram1D* _hist_mm_Jeteta2j;
        AIDA::IHistogram1D* _hist_mm_Jeteta3j;
        AIDA::IHistogram1D* _hist_mm_Jeteta4j;
        AIDA::IHistogram1D* _hist_mm_Jeteta5j;
        AIDA::IHistogram1D* _hist_mm_Jeteta6j;

        AIDA::IHistogram1D* _hist_ee_Ht_1j;
        AIDA::IHistogram1D* _hist_ee_Ht_2j;
        AIDA::IHistogram1D* _hist_ee_Ht_3j;
        AIDA::IHistogram1D* _hist_ee_Ht_4j;
        AIDA::IHistogram1D* _hist_ee_Ht_5j;

        AIDA::IHistogram1D* _hist_mm_Ht_1j;
        AIDA::IHistogram1D* _hist_mm_Ht_2j;
        AIDA::IHistogram1D* _hist_mm_Ht_3j;
        AIDA::IHistogram1D* _hist_mm_Ht_4j;
        AIDA::IHistogram1D* _hist_mm_Ht_5j;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_SMP_13_007> plugin_CMS_SMP_13_007;


}
                                                                                                                             421,1         Bot
}


