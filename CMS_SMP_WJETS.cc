// -*- C++ -*-

//**********Bhawandeep**********


//********WJets@8TeV************


#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, WFinder...
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"

#include <iostream>

namespace Rivet {
    
    // define bool used in sorting
    bool orderByIncRap(const FourMomentum& a, const FourMomentum& b) {
        return (a.rapidity() < b.rapidity());
    }
    
    class CMS_SMP_WJETS: public Analysis {
    public:
        
        // Constructors
        CMS_SMP_WJETS() : Analysis("CMS_SMP_WJETS") {
            setBeams(PROTON, PROTON);
            setNeedsCrossSection(true);
        }
        
    public:
        
        // Book histograms and initialise projections before the run
        void init() {
            
            FinalState fs;
            //WFinder wfinder_mu(fs, -2.4, 2.4, 0*GeV, MUON, 50*GeV, 1000000*GeV, 0*GeV, 0.1, 1, true, false);
            //WFinder wfinder_mu(fs, -2.4, 2.4, 0*GeV, MUON, 50*GeV, 1000000*GeV, 0*GeV, 0.1, true, false, 80.4, true);
            WFinder wfinder_mu(fs, -2.4, 2.4, 0*GeV, MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, true, false);
            addProjection(wfinder_mu, "WFinder_mu");
            
            
            // Define veto FS
            VetoedFinalState vfs;
            vfs.addVetoOnThisFinalState(wfinder_mu);
            vfs.addVetoPairId(MUON);
            vfs.vetoNeutrinos();
            
            FastJets fastjets(vfs, FastJets::ANTIKT, 0.5);
            addProjection(fastjets, "Jets");
            
            // define binning and book Histos
            vector<double> addjetPt_Winc1jet ;
            vector<double> addjetPt_Winc2jet ;
            vector<double> addjetPt_Winc3jet ;
            vector<double> addjetPt_Winc4jet ;
            vector<double> addjetPt_Winc5jet ;
            
            vector<double> jetPt_Winc1jet ;
            vector<double> jetPt_Winc2jet ;
            vector<double> jetPt_Winc3jet ;
            vector<double> jetPt_Winc4jet ;
            vector<double> jetPt_Winc5jet ;
            
            vector<double> addjetHT_Winc1jet ;
            vector<double> addjetHT_Winc2jet ;
            vector<double> addjetHT_Winc3jet ;
            vector<double> addjetHT_Winc4jet ;
            vector<double> addjetHT_Winc5jet ;
            
            vector<double> jetHT_Winc1jet ;
            vector<double> jetHT_Winc2jet ;
            vector<double> jetHT_Winc3jet ;
            vector<double> jetHT_Winc4jet ;
            vector<double> jetHT_Winc5jet ;
            
            vector<double> dijetM_Winc2jet ;
            vector<double> dijetM_Winc3jet ;
            vector<double> dijetM_Winc4jet ;
            
            addjetPt_Winc1jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000;
            addjetPt_Winc2jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 800;
            addjetPt_Winc3jet += 20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300;
            addjetPt_Winc4jet += 20, 24, 30, 39, 49, 62, 78, 96, 150;
            addjetPt_Winc5jet += 20, 24, 30, 39, 49, 62, 100;
            
            jetPt_Winc1jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 900, 1400;
            jetPt_Winc2jet += 20, 24, 30, 39, 49, 60, 72, 85, 100, 117, 136, 157, 187, 220, 258, 300, 350, 400, 450, 500, 590, 700, 1000;
            jetPt_Winc3jet += 20, 24, 30, 39, 49, 62, 78, 105, 142, 185, 235, 300, 500;
            jetPt_Winc4jet += 20, 24, 30, 39, 49, 62, 78, 96, 120, 160, 250;
            jetPt_Winc5jet += 20, 24, 30, 39, 49, 62, 90, 140;
            
            addjetHT_Winc1jet += 30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1500;
            addjetHT_Winc2jet += 60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1200;
            addjetHT_Winc3jet += 90, 105, 125, 151, 185, 230, 290, 366, 466, 586, 767, 990;
            addjetHT_Winc4jet += 120, 140, 167, 203, 253, 320, 410, 530, 690, 910;
            addjetHT_Winc5jet += 150, 180, 222, 282, 365, 485, 650, 880;
            
            jetHT_Winc1jet += 30, 39, 49, 62, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1300, 2000;
            jetHT_Winc2jet += 60, 78, 96, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1000, 1300, 2000;
            jetHT_Winc3jet += 90, 118, 150, 190, 240, 300, 370, 450, 540, 650, 800, 1100, 1800;
            jetHT_Winc4jet += 120, 140, 165, 200, 250, 310, 380, 480, 600, 800, 1100, 1800;
            jetHT_Winc5jet += 150, 180, 222, 282, 365, 485, 650, 880, 1300;
            
            dijetM_Winc2jet += 0, 25, 52, 81, 112, 145, 180, 217, 256, 297, 340, 385, 432, 481, 532, 585, 640, 700;
            dijetM_Winc3jet += 0, 30, 62, 96, 132, 170, 210, 252, 296, 342, 390, 440, 492, 546, 602, 660;
            //dijetM_Winc4jet += 0, 34, 70, 108, 148, 190, 234, 280, 328, 278, 430, 484, 540, 598, 660;
            dijetM_Winc4jet += 0, 34, 70, 108, 148, 190, 234, 280, 328, 378, 430, 484, 540, 598, 660; 

            //-------------
            _hist_addJetPt1j =bookHistogram1D("addjet_Pt1jetcase", addjetPt_Winc1jet);
            _hist_addJetPt2j =bookHistogram1D("addjet_Pt2jetcase", addjetPt_Winc2jet);
            _hist_addJetPt3j =bookHistogram1D("addjet_Pt3jetcase", addjetPt_Winc3jet);
            _hist_addJetPt4j =bookHistogram1D("addjet_Pt4jetcase", addjetPt_Winc4jet);
            _hist_addJetPt5j =bookHistogram1D("addjet_Pt5jetcase", addjetPt_Winc5jet);
            
            _hist_JetPt1j =bookHistogram1D("jet_Pt1jetcase", jetPt_Winc1jet);
            _hist_JetPt2j =bookHistogram1D("jet_Pt2jetcase", jetPt_Winc2jet);
            _hist_JetPt3j =bookHistogram1D("jet_Pt3jetcase", jetPt_Winc3jet);
            _hist_JetPt4j =bookHistogram1D("jet_Pt4jetcase", jetPt_Winc4jet);
            _hist_JetPt5j =bookHistogram1D("jet_Pt5jetcase", jetPt_Winc5jet);
            
            //-------------
            _hist_addHt_1j =bookHistogram1D("addJetsHT_inc1jet", addjetHT_Winc1jet);
            _hist_addHt_2j =bookHistogram1D("addJetsHT_inc2jet", addjetHT_Winc2jet);
            _hist_addHt_3j =bookHistogram1D("addJetsHT_inc3jet", addjetHT_Winc3jet);
            _hist_addHt_4j =bookHistogram1D("addJetsHT_inc4jet", addjetHT_Winc4jet);
            _hist_addHt_5j =bookHistogram1D("addJetsHT_inc5jet", addjetHT_Winc5jet);
            
            _hist_Ht_1j =bookHistogram1D("JetsHT_inc1jet", jetHT_Winc1jet);
            _hist_Ht_2j =bookHistogram1D("JetsHT_inc2jet", jetHT_Winc2jet);
            _hist_Ht_3j =bookHistogram1D("JetsHT_inc3jet", jetHT_Winc3jet);
            _hist_Ht_4j =bookHistogram1D("JetsHT_inc4jet", jetHT_Winc4jet);
            _hist_Ht_5j =bookHistogram1D("JetsHT_inc5jet", jetHT_Winc5jet);
            
            //-------------
            _hist_inc_WJetMult = bookHistogram1D("njetWJet_incl", 8, -0.5, 7.5);
            _hist_excl_WJetMult= bookHistogram1D("njetWJet_excl", 8, -0.5, 7.5);
            _hist_Mult_exc = bookHistogram1D("njet_exc_fbin", 8, -0.5, 7.5);
            
            //-------------
            _hist_Jeteta1j =bookHistogram1D("jet_eta1jetcase", 32, 0, 2.4);
            _hist_Jeteta2j =bookHistogram1D("jet_eta2jetcase", 32, 0, 2.4);
            _hist_Jeteta3j =bookHistogram1D("jet_eta3jetcase", 12, 0, 2.4);
            _hist_Jeteta4j =bookHistogram1D("jet_eta4jetcase", 12, 0, 2.4);
            _hist_Jeteta5j =bookHistogram1D("jet_eta5jetcase", 6,  0, 2.4);
            
            //-------------
            _hist_dyj1j2_2j =bookHistogram1D("dyj1j2_inc2jet", 20, 0, 4.8);
            _hist_dyj1j2_3j =bookHistogram1D("dyj1j2_inc3jet", 20, 0, 4.8);
            _hist_dyj1j2_4j =bookHistogram1D("dyj1j2_inc4jet", 16, 0, 4.8);
            
            _hist_dyjFjB_2j =bookHistogram1D("dyjFjB_inc2jet", 20, 0, 4.8);
            _hist_dyjFjB_3j =bookHistogram1D("dyjFjB_inc3jet", 20, 0, 4.8);
            _hist_dyjFjB_4j =bookHistogram1D("dyjFjB_inc4jet", 16, 0, 4.8);
            
            _hist_dyj1j3_3j =bookHistogram1D("dyj1j3_inc3jet", 20, 0, 4.8);
            _hist_dyj2j3_3j =bookHistogram1D("dyj2j3_inc3jet", 20, 0, 4.8);
            
            _hist_dphij1j2_2j =bookHistogram1D("dphij1j2_inc2jet", 20, 0, 3.14159265359);
            _hist_dphijFjB_2j =bookHistogram1D("dphijFjB_inc2jet", 20, 0, 3.14159265359);
            
            _hist_dRj1j2_2j =bookHistogram1D("dRj1j2_inc2jet", 30, 0, 6.);
            
            _hist_dijetM_2j =bookHistogram1D("dijetM_inc2jet", dijetM_Winc2jet);
            _hist_dijetM_3j =bookHistogram1D("dijetM_inc3jet", dijetM_Winc3jet);
            _hist_dijetM_4j =bookHistogram1D("dijetM_inc4jet", dijetM_Winc4jet);
            
            _hist_diJetPt_2j =bookHistogram1D("diJetPt_inc2jet", addjetPt_Winc2jet);
            _hist_diJetPt_3j =bookHistogram1D("diJetPt_inc3jet", addjetPt_Winc3jet);
            _hist_diJetPt_4j =bookHistogram1D("diJetPt_inc4jet", addjetPt_Winc4jet);
            
            //-------------
            _hist_dphij1mu_1j =bookHistogram1D("dphij1mu_inc1jet", 20, 0, 3.14159265359);
            _hist_dphij2mu_2j =bookHistogram1D("dphij2mu_inc2jet", 20, 0, 3.14159265359);
            _hist_dphij3mu_3j =bookHistogram1D("dphij3mu_inc3jet", 16, 0, 3.14159265359);
            _hist_dphij4mu_4j =bookHistogram1D("dphij4mu_inc4jet", 16, 0, 3.14159265359);
            _hist_dphij5mu_5j =bookHistogram1D("dphij5mu_inc5jet", 12, 0, 3.14159265359);
            
            
            //------- MeanNJ ------
            _hist_MeanNJht_1j =bookHistogram1D("MeanNJht_inc1jet", addjetHT_Winc1jet);
            _hist_MeanNJht_2j =bookHistogram1D("MeanNJht_inc2jet", addjetHT_Winc2jet);
            _hist_MeanNJdyj1j2_2j = bookHistogram1D("MeanNJdyj1j2_inc2jet", 20, 0, 4.8);
            _hist_MeanNJdyjFjB_2j = bookHistogram1D("MeanNJdyjFjB_inc2jet", 20, 0, 4.8);
            
            //-------------------------------------
            _hist_ht_Exc1jbin1 =bookHistogram1D("ht_Exc1jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc2jbin1 =bookHistogram1D("ht_Exc2jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc3jbin1 =bookHistogram1D("ht_Exc3jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc4jbin1 =bookHistogram1D("ht_Exc4jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc5jbin1 =bookHistogram1D("ht_Exc5jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc6jbin1 =bookHistogram1D("ht_Exc6jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc7jbin1 =bookHistogram1D("ht_Exc7jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc8jbin1 =bookHistogram1D("ht_Exc8jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc9jbin1 =bookHistogram1D("ht_Exc9jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc10jbin1 =bookHistogram1D("ht_Exc10jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc11jbin1 =bookHistogram1D("ht_Exc11jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc12jbin1 =bookHistogram1D("ht_Exc12jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc13jbin1 =bookHistogram1D("ht_Exc13jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc14jbin1 =bookHistogram1D("ht_Exc14jetbin1", addjetHT_Winc1jet);
            _hist_ht_Exc15jbin1 =bookHistogram1D("ht_Exc15jetbin1", addjetHT_Winc1jet);
            
            _hist_ht_Exc2jbin2 =bookHistogram1D("ht_Exc2jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc3jbin2 =bookHistogram1D("ht_Exc3jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc4jbin2 =bookHistogram1D("ht_Exc4jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc5jbin2 =bookHistogram1D("ht_Exc5jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc6jbin2 =bookHistogram1D("ht_Exc6jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc7jbin2 =bookHistogram1D("ht_Exc7jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc8jbin2 =bookHistogram1D("ht_Exc8jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc9jbin2 =bookHistogram1D("ht_Exc9jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc10jbin2 =bookHistogram1D("ht_Exc10jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc11jbin2 =bookHistogram1D("ht_Exc11jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc12jbin2 =bookHistogram1D("ht_Exc12jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc13jbin2 =bookHistogram1D("ht_Exc13jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc14jbin2 =bookHistogram1D("ht_Exc14jetbin2", addjetHT_Winc2jet);
            _hist_ht_Exc15jbin2 =bookHistogram1D("ht_Exc15jetbin2", addjetHT_Winc2jet);
            
            _hist_dyj1j2_Exc2j = bookHistogram1D("dyj1j2_Exc2jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc3j = bookHistogram1D("dyj1j2_Exc3jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc4j = bookHistogram1D("dyj1j2_Exc4jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc5j = bookHistogram1D("dyj1j2_Exc5jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc6j = bookHistogram1D("dyj1j2_Exc6jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc7j = bookHistogram1D("dyj1j2_Exc7jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc8j = bookHistogram1D("dyj1j2_Exc8jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc9j = bookHistogram1D("dyj1j2_Exc9jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc10j = bookHistogram1D("dyj1j2_Exc10jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc11j = bookHistogram1D("dyj1j2_Exc11jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc12j = bookHistogram1D("dyj1j2_Exc12jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc13j = bookHistogram1D("dyj1j2_Exc13jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc14j = bookHistogram1D("dyj1j2_Exc14jet", 20, 0, 4.8);
            _hist_dyj1j2_Exc15j = bookHistogram1D("dyj1j2_Exc15jet", 20, 0, 4.8);
            
            _hist_dyjFjB_Exc2j = bookHistogram1D("dyjFjB_Exc2jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc3j = bookHistogram1D("dyjFjB_Exc3jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc4j = bookHistogram1D("dyjFjB_Exc4jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc5j = bookHistogram1D("dyjFjB_Exc5jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc6j = bookHistogram1D("dyjFjB_Exc6jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc7j = bookHistogram1D("dyjFjB_Exc7jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc8j = bookHistogram1D("dyjFjB_Exc8jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc9j = bookHistogram1D("dyjFjB_Exc9jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc10j = bookHistogram1D("dyjFjB_Exc10jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc11j = bookHistogram1D("dyjFjB_Exc11jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc12j = bookHistogram1D("dyjFjB_Exc12jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc13j = bookHistogram1D("dyjFjB_Exc13jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc14j = bookHistogram1D("dyjFjB_Exc14jet", 20, 0, 4.8);
            _hist_dyjFjB_Exc15j = bookHistogram1D("dyjFjB_Exc15jet", 20, 0, 4.8);
        }
        
        // define function used for filiing inc Njets histo
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
            const WFinder& wfinder_mu = applyProjection<WFinder>(event, "WFinder_mu");
            
            if (wfinder_mu.bosons().size() != 1) {
                vetoEvent;
            }
            
            if (wfinder_mu.bosons().size() == 1) {
                
                const FourMomentum& lepton0 = wfinder_mu.constituentLeptons()[0].momentum();
                const FourMomentum& neutrino = wfinder_mu.constituentNeutrinos()[0].momentum();
                double WmT = sqrt( 2 * lepton0.pT() * neutrino.pT() * (1 - cos(deltaPhi(lepton0, neutrino))) );
                
                if (WmT < 50.0*GeV) vetoEvent;
                
                double pt0 = lepton0.pT();
                double eta0 = lepton0.eta();
                
                if ( (fabs(eta0) > 2.1) || (pt0 < 25.0*GeV) ) vetoEvent;
                
                //Obtain the jets::::::::::::::
                vector<FourMomentum> finaljet_list;
                double HT = 0.0;
                
                // loop over jets in a event, pushback in finaljet_list collection
                foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
                    
                    const double jeta = j.momentum().eta();
                    const double jpt = j.momentum().pT();
                    
                    if ( (fabs(jeta) < 2.4) && (deltaR(lepton0, j.momentum()) > 0.5) ) {
                        
                        if(jpt > 30.0*GeV) {
                            finaljet_list.push_back(j.momentum());
                            HT += j.momentum().pT();
                        }
                    }
                } // end looping over jets
                
                // new jet_list sorted by increasing rapidity
                vector<FourMomentum> jListRap = finaljet_list;
                std::sort(jListRap.begin(), jListRap.end(), orderByIncRap);
                
                // Multiplicity exc plot.
                _hist_excl_WJetMult->fill(finaljet_list.size(), weight);
                
                if(finaljet_list.size()<=7) {
                    _hist_Mult_exc->fill(finaljet_list.size(), weight);
                }
                else if (finaljet_list.size()>7){
                    _hist_Mult_exc->fill(7., weight);
                }
                
                
                
                // Multiplicity inc plot.
                Fill(_hist_inc_WJetMult, weight, finaljet_list);
                
                if(finaljet_list.size()>=1) {
                    _hist_addJetPt1j->fill(finaljet_list[0].pT(), weight);
                    _hist_JetPt1j->fill(finaljet_list[0].pT(), weight);
                    _hist_Jeteta1j->fill(fabs(finaljet_list[0].eta()), weight);
                    _hist_addHt_1j->fill(HT, weight);
                    _hist_Ht_1j->fill(HT, weight);
                    _hist_dphij1mu_1j->fill( deltaPhi(finaljet_list[0].phi(), lepton0.phi()), weight );
                    _hist_MeanNJht_1j->fill( HT, weight*finaljet_list.size());
                }
                
                if(finaljet_list.size()>=2) {
                    _hist_addJetPt2j->fill(finaljet_list[1].pT(), weight);
                    _hist_JetPt2j->fill(finaljet_list[1].pT(), weight);
                    _hist_Jeteta2j->fill(fabs(finaljet_list[1].eta()), weight);
                    _hist_addHt_2j->fill(HT, weight);
                    _hist_Ht_2j->fill(HT, weight);
                    
                    _hist_dyj1j2_2j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_2j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                    _hist_dphij1j2_2j ->fill( deltaPhi(finaljet_list[0].phi(), finaljet_list[1].phi()), weight);
                    _hist_dphijFjB_2j ->fill( deltaPhi(jListRap[0].phi(), jListRap[jListRap.size()-1].phi()) , weight);
                    
                    _hist_dijetM_2j   ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
                    _hist_diJetPt_2j  ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
                    _hist_dRj1j2_2j   ->fill( deltaR(finaljet_list[0].rapidity(), finaljet_list[0].phi(), finaljet_list[1].rapidity(), finaljet_list[1].phi()), weight);
                    
                    _hist_dphij2mu_2j ->fill( deltaPhi(finaljet_list[1].phi(), lepton0.phi()), weight );
                    
                    _hist_MeanNJht_2j->fill( HT, weight*finaljet_list.size());
                    _hist_MeanNJdyj1j2_2j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight*finaljet_list.size());
                    _hist_MeanNJdyjFjB_2j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight*finaljet_list.size());
                }
                
                if(finaljet_list.size()>=3) {
                    _hist_addJetPt3j->fill(finaljet_list[2].pT(), weight);
                    _hist_JetPt3j->fill(finaljet_list[2].pT(), weight);
                    _hist_Jeteta3j->fill(fabs(finaljet_list[2].eta()), weight);
                    _hist_addHt_3j->fill(HT, weight);
                    _hist_Ht_3j->fill(HT, weight);
                    
                    _hist_dyj1j2_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyj1j3_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[2].rapidity()), weight);
                    _hist_dyj2j3_3j     ->fill( fabs(finaljet_list[1].rapidity() - finaljet_list[2].rapidity()), weight);
                    _hist_dyjFjB_3j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                    
                    _hist_dijetM_3j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
                    _hist_diJetPt_3j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
                    
                    _hist_dphij3mu_3j->fill( deltaPhi(finaljet_list[2].phi(), lepton0.phi()), weight );
                }
                
                if(finaljet_list.size()>=4) {
                    _hist_addJetPt4j->fill(finaljet_list[3].pT(), weight);
                    _hist_JetPt4j->fill(finaljet_list[3].pT(), weight);
                    _hist_Jeteta4j->fill(fabs(finaljet_list[3].eta()), weight);
                    _hist_addHt_4j->fill(HT, weight);
                    _hist_Ht_4j->fill(HT, weight);
                    
                    _hist_dyj1j2_4j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_4j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                    
                    _hist_dijetM_4j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass(), weight);
                    _hist_diJetPt_4j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT(), weight);
                    _hist_dphij4mu_4j->fill( deltaPhi(finaljet_list[3].phi(), lepton0.phi()), weight );
                }
                
                if(finaljet_list.size()>=5) {
                    _hist_addJetPt5j->fill(finaljet_list[4].pT(), weight);
                    _hist_JetPt5j->fill(finaljet_list[4].pT(), weight);
                    _hist_Jeteta5j->fill(fabs(finaljet_list[4].eta()), weight);
                    _hist_addHt_5j->fill(HT, weight);
                    _hist_Ht_5j->fill(HT, weight);
                    _hist_dphij5mu_5j->fill( deltaPhi(finaljet_list[4].phi(), lepton0.phi()), weight );
                }
                
                
                //-------------------------------------
                if(finaljet_list.size()==1){
                    _hist_ht_Exc1jbin1->fill(HT, weight);
                }
                if(finaljet_list.size()==2){
                    _hist_ht_Exc2jbin1->fill(HT, weight);
                    _hist_ht_Exc2jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc2j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc2j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==3){
                    _hist_ht_Exc3jbin1->fill(HT, weight);
                    _hist_ht_Exc3jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc3j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc3j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==4){
                    _hist_ht_Exc4jbin1->fill(HT, weight);
                    _hist_ht_Exc4jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc4j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc4j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==5){
                    _hist_ht_Exc5jbin1->fill(HT, weight);
                    _hist_ht_Exc5jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc5j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc5j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==6){
                    _hist_ht_Exc6jbin1->fill(HT, weight);
                    _hist_ht_Exc6jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc6j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc6j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==7){
                    _hist_ht_Exc7jbin1->fill(HT, weight);
                    _hist_ht_Exc7jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc7j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc7j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==8){
                    _hist_ht_Exc8jbin1->fill(HT, weight);
                    _hist_ht_Exc8jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc8j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc8j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==9){
                    _hist_ht_Exc9jbin1->fill(HT, weight);
                    _hist_ht_Exc9jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc9j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc9j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==10){
                    _hist_ht_Exc10jbin1->fill(HT, weight);
                    _hist_ht_Exc10jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc10j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc10j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==11){
                    _hist_ht_Exc11jbin1->fill(HT, weight);
                    _hist_ht_Exc11jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc11j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc11j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==12){
                    _hist_ht_Exc12jbin1->fill(HT, weight);
                    _hist_ht_Exc12jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc12j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc12j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==13){
                    _hist_ht_Exc13jbin1->fill(HT, weight);
                    _hist_ht_Exc13jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc13j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc13j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==14){
                    _hist_ht_Exc14jbin1->fill(HT, weight);
                    _hist_ht_Exc14jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc14j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc14j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                if(finaljet_list.size()==15){
                    _hist_ht_Exc15jbin1->fill(HT, weight);
                    _hist_ht_Exc15jbin2->fill(HT, weight);
                    _hist_dyj1j2_Exc15j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), weight);
                    _hist_dyjFjB_Exc15j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), weight);
                }
                
            } //////
            
        } //void loop
        
        
        /// Normalise histograms etc., after the run
        void finalize() {
            
        //    double crossec=36702;
          
double crossec =36280;



//cout<<"crossection======="<<crossec<<endl;


  scale(_hist_inc_WJetMult, crossec/sumOfWeights());
            scale(_hist_excl_WJetMult, crossec/sumOfWeights());
            scale(_hist_Mult_exc, crossec/sumOfWeights());
            
            scale(_hist_addJetPt1j, crossec/sumOfWeights());
            scale(_hist_addJetPt2j, crossec/sumOfWeights());
            scale(_hist_addJetPt3j, crossec/sumOfWeights());
            scale(_hist_addJetPt4j, crossec/sumOfWeights());
            scale(_hist_addJetPt5j, crossec/sumOfWeights());
            
            scale(_hist_JetPt1j, crossec/sumOfWeights());
            scale(_hist_JetPt2j, crossec/sumOfWeights());
            scale(_hist_JetPt3j, crossec/sumOfWeights());
            scale(_hist_JetPt4j, crossec/sumOfWeights());
            scale(_hist_JetPt5j, crossec/sumOfWeights());
            
            scale(_hist_Jeteta1j, crossec/sumOfWeights());
            scale(_hist_Jeteta2j, crossec/sumOfWeights());
            scale(_hist_Jeteta3j, crossec/sumOfWeights());
            scale(_hist_Jeteta4j, crossec/sumOfWeights());
            scale(_hist_Jeteta5j, crossec/sumOfWeights());
            
            scale(_hist_addHt_1j, crossec/sumOfWeights());
            scale(_hist_addHt_2j, crossec/sumOfWeights());
            scale(_hist_addHt_3j, crossec/sumOfWeights());
            scale(_hist_addHt_4j, crossec/sumOfWeights());
            scale(_hist_addHt_5j, crossec/sumOfWeights());
            
            scale(_hist_Ht_1j, crossec/sumOfWeights());
            scale(_hist_Ht_2j, crossec/sumOfWeights());
            scale(_hist_Ht_3j, crossec/sumOfWeights());
            scale(_hist_Ht_4j, crossec/sumOfWeights());
            scale(_hist_Ht_5j, crossec/sumOfWeights());
            
            //-------------------------------------
            scale(_hist_dyj1j2_2j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_3j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_4j, crossec/sumOfWeights());
            
            scale(_hist_dyjFjB_2j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_3j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_4j, crossec/sumOfWeights());
            
            scale(_hist_dyj1j3_3j, crossec/sumOfWeights());
            scale(_hist_dyj2j3_3j, crossec/sumOfWeights());
            
            scale(_hist_dphij1j2_2j, crossec/sumOfWeights());
            scale(_hist_dphijFjB_2j, crossec/sumOfWeights());
            
            scale(_hist_dRj1j2_2j, crossec/sumOfWeights());
            
            scale(_hist_dijetM_2j, crossec/sumOfWeights());
            scale(_hist_dijetM_3j, crossec/sumOfWeights());
            scale(_hist_dijetM_4j, crossec/sumOfWeights());
            
            scale(_hist_diJetPt_2j, crossec/sumOfWeights());
            scale(_hist_diJetPt_3j, crossec/sumOfWeights());
            scale(_hist_diJetPt_4j, crossec/sumOfWeights());
            
            scale(_hist_dphij1mu_1j, crossec/sumOfWeights());
            scale(_hist_dphij2mu_2j, crossec/sumOfWeights());
            scale(_hist_dphij3mu_3j, crossec/sumOfWeights());
            scale(_hist_dphij4mu_4j, crossec/sumOfWeights());
            scale(_hist_dphij5mu_5j, crossec/sumOfWeights());
            
            //-------------------------------------
            scale(_hist_MeanNJht_1j, crossec/sumOfWeights());
            scale(_hist_MeanNJht_2j, crossec/sumOfWeights());
            scale(_hist_MeanNJdyj1j2_2j, crossec/sumOfWeights());
            scale(_hist_MeanNJdyjFjB_2j, crossec/sumOfWeights());
            
            //-------------------------------------
            scale(_hist_ht_Exc1jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc2jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc3jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc4jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc5jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc6jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc7jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc8jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc9jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc10jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc11jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc12jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc13jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc14jbin1, crossec/sumOfWeights());
            scale(_hist_ht_Exc15jbin1, crossec/sumOfWeights());
            
            scale(_hist_ht_Exc2jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc3jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc4jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc5jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc6jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc7jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc8jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc9jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc10jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc11jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc12jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc13jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc14jbin2, crossec/sumOfWeights());
            scale(_hist_ht_Exc15jbin2, crossec/sumOfWeights());
            
            scale(_hist_dyj1j2_Exc2j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc3j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc4j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc5j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc6j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc7j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc8j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc9j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc10j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc11j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc12j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc13j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc14j, crossec/sumOfWeights());
            scale(_hist_dyj1j2_Exc15j, crossec/sumOfWeights());
            
            scale(_hist_dyjFjB_Exc2j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc3j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc4j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc5j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc6j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc7j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc8j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc9j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc10j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc11j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc12j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc13j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc14j, crossec/sumOfWeights());
            scale(_hist_dyjFjB_Exc15j, crossec/sumOfWeights());
            
        }
        
    private:
        
        // Data members like post-cuts event weight counters go here
        
    private:
        
        AIDA::IHistogram1D* _hist_inc_WJetMult;
        AIDA::IHistogram1D* _hist_excl_WJetMult;
        AIDA::IHistogram1D* _hist_Mult_exc;
        
        AIDA::IHistogram1D* _hist_addJetPt1j;
        AIDA::IHistogram1D* _hist_addJetPt2j;
        AIDA::IHistogram1D* _hist_addJetPt3j;
        AIDA::IHistogram1D* _hist_addJetPt4j;
        AIDA::IHistogram1D* _hist_addJetPt5j;
        
        AIDA::IHistogram1D* _hist_JetPt1j;
        AIDA::IHistogram1D* _hist_JetPt2j;
        AIDA::IHistogram1D* _hist_JetPt3j;
        AIDA::IHistogram1D* _hist_JetPt4j;
        AIDA::IHistogram1D* _hist_JetPt5j;
        
        AIDA::IHistogram1D* _hist_Jeteta1j;
        AIDA::IHistogram1D* _hist_Jeteta2j;
        AIDA::IHistogram1D* _hist_Jeteta3j;
        AIDA::IHistogram1D* _hist_Jeteta4j;
        AIDA::IHistogram1D* _hist_Jeteta5j;
        AIDA::IHistogram1D* _hist_Jeteta6j;
        
        AIDA::IHistogram1D* _hist_addHt_1j;
        AIDA::IHistogram1D* _hist_addHt_2j;
        AIDA::IHistogram1D* _hist_addHt_3j;
        AIDA::IHistogram1D* _hist_addHt_4j;
        AIDA::IHistogram1D* _hist_addHt_5j;
        
        AIDA::IHistogram1D* _hist_Ht_1j;
        AIDA::IHistogram1D* _hist_Ht_2j;
        AIDA::IHistogram1D* _hist_Ht_3j;
        AIDA::IHistogram1D* _hist_Ht_4j;
        AIDA::IHistogram1D* _hist_Ht_5j;
        
        //-------------------------------------
        AIDA::IHistogram1D* _hist_dyj1j2_2j;
        AIDA::IHistogram1D* _hist_dyj1j2_3j;
        AIDA::IHistogram1D* _hist_dyj1j2_4j;
        
        AIDA::IHistogram1D* _hist_dyjFjB_2j;
        AIDA::IHistogram1D* _hist_dyjFjB_3j;
        AIDA::IHistogram1D* _hist_dyjFjB_4j;
        
        AIDA::IHistogram1D* _hist_dyj1j3_3j;
        AIDA::IHistogram1D* _hist_dyj2j3_3j;
        
        AIDA::IHistogram1D* _hist_dphij1j2_2j;
        AIDA::IHistogram1D* _hist_dphijFjB_2j;
        
        AIDA::IHistogram1D* _hist_dRj1j2_2j;
        
        AIDA::IHistogram1D* _hist_dijetM_2j;
        AIDA::IHistogram1D* _hist_dijetM_3j;
        AIDA::IHistogram1D* _hist_dijetM_4j;
        
        AIDA::IHistogram1D* _hist_diJetPt_2j;
        AIDA::IHistogram1D* _hist_diJetPt_3j;
        AIDA::IHistogram1D* _hist_diJetPt_4j;
        
        AIDA::IHistogram1D* _hist_dphij1mu_1j;
        AIDA::IHistogram1D* _hist_dphij2mu_2j;
        AIDA::IHistogram1D* _hist_dphij3mu_3j;
        AIDA::IHistogram1D* _hist_dphij4mu_4j;
        AIDA::IHistogram1D* _hist_dphij5mu_5j;
        
        //-------------------------------------
        AIDA::IHistogram1D* _hist_MeanNJht_1j;
        AIDA::IHistogram1D* _hist_MeanNJht_2j;
        AIDA::IHistogram1D* _hist_MeanNJdyj1j2_2j;
        AIDA::IHistogram1D* _hist_MeanNJdyjFjB_2j;
        
        //-------------------------------------
        AIDA::IHistogram1D* _hist_ht_Exc1jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc2jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc3jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc4jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc5jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc6jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc7jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc8jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc9jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc10jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc11jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc12jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc13jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc14jbin1;
        AIDA::IHistogram1D* _hist_ht_Exc15jbin1;
        
        AIDA::IHistogram1D* _hist_ht_Exc2jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc3jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc4jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc5jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc6jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc7jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc8jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc9jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc10jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc11jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc12jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc13jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc14jbin2;
        AIDA::IHistogram1D* _hist_ht_Exc15jbin2;
        
        AIDA::IHistogram1D* _hist_dyj1j2_Exc2j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc3j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc4j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc5j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc6j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc7j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc8j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc9j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc10j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc11j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc12j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc13j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc14j;
        AIDA::IHistogram1D* _hist_dyj1j2_Exc15j;
        
        AIDA::IHistogram1D* _hist_dyjFjB_Exc2j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc3j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc4j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc5j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc6j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc7j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc8j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc9j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc10j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc11j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc12j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc13j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc14j;
        AIDA::IHistogram1D* _hist_dyjFjB_Exc15j;
        
    };
    
    // This global object acts as a hook for the plugin system
    AnalysisBuilder<CMS_SMP_WJETS> plugin_CMS_SMP_WJETS;
    
}




