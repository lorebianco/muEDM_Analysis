#ifndef TRACKDATAMANAGER_HH
#define TRACKDATAMANAGER_HH

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVector3.h>

#include "CHeT/CHeTReader.hh"
#include "config.hh"

class TrackDataManager
{
  public:
    TrackDataManager();
    ~TrackDataManager();

    inline Int_t GetNEvents() const
    {
        return nEvents;
    };
    inline Int_t GetProcessedEvents() const
    {
        return processedEvents;
    };

    void ProcessAndFilterEvent(Long64_t eventID = -1);
    inline TChain *GetChain() const
    {
        return simChain;
    }
    Bool_t SetTrueDecayData(Int_t eventID);
    void InitRecoTree(bool isMichel);
    void SaveRecoTree();

    // Data (Linked to input TTree branches)
    Int_t truedecay_EventID;
    Double_t trueDecay_posX, trueDecay_posY, trueDecay_posZ, trueDecay_momX, trueDecay_momY,
        trueDecay_momZ, trueDecay_time;

    std::vector<Int_t> *all_hits;

    // Data (for analysis)
    Float_t trueDecayTime;
    TVector3 trueDecayPos, trueDecayMom;
    Double_t trueDecayThetaEDM;
    std::vector<std::vector<Double_t>> hitsCoordinates;
    std::vector<Int_t> hitsCylinderID;

    // Fitted Data
    TVector3 fittedDecayPos, fittedDecayMom;
    Double_t fittedDecayThetaEDM;

    // Output Data
    TFile *outputFile = nullptr;
    TTree *recTree = nullptr;

    // Cosmic rec parameters
    Double_t rec_x0, rec_z0, rec_sx, rec_sz;

    // Michel rec parameters
    Double_t rec_R, rec_cx, rec_cy, rec_dz_ds, rec_phi0, rec_t_min, rec_t_max;

    // Extrapolation parameters
    Double_t rec_extrap_x, rec_extrap_y, rec_extrap_z;
    Double_t rec_extrap_px, rec_extrap_py, rec_extrap_pz;

    // Additional parameters
    Double_t rec_chi2;
    Bool_t rec_converged;
    Int_t rec_n_candidates_2d, rec_n_candidates_z;

    // Flags
    Bool_t is_cosmic_valid;
    Bool_t is_michel_valid;

    std::vector<Int_t> rec_hits;
    std::vector<Int_t> rec_hough2d_idx;
    std::vector<Int_t> rec_houghz_idx;

    // Double_t trueMomentum;
    // Double_t polarAngle;
    // Double_t azimuthalAngle;
    // Double_t theThetaAngle;
    // Double_t spinAngle;
    // Double_t emissionAngle;
    // TVector3* fOrigin;
    // std::vector<std::vector<Double_t>>* hitsCoordinates;
    // std::vector<std::vector<Double_t>>* trackCoordinates;
    // std::vector<Int_t>* cylinderID;

    // Histos
    TEfficiency *accPhi, *accTheta, *effPhi, *effTheta;
    TH1I *histTurns, *histCylinders, *histFakeHits;
    TEfficiency *effTurns, *effCylinders;
    TProfile *histCylVMom, *histTurnsVMom;
    TH2D *graphMom, *graphTheta, *graphPhi;
    TH1D *histDiffX, *histDiffY, *histDiffZ;
    TH1D *histDiffMom, *histDiffTheta, *histDiffPhi;
    TH2D *hist2MomRes, *hist2ThetaRes, *hist2PhiRes;
    TProfile *profMomRes, *profThetaRes, *profPhiRes;
    TH1F *histTime;

  private:
    TChain *simChain;
    std::unique_ptr<CHeT::Data::Reader> reader;
    std::vector<ROOT::VecOps::RVec<int>> allEventsHits;

    Int_t nEvents, processedEvents;
    Long64_t currentEventIndex = 0;
};

#endif // TRACKDATAMANAGER_HH
