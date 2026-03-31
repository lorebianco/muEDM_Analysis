#ifndef TRACKDATAMANAGER_HH
#define TRACKDATAMANAGER_HH

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include <TChain.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TMath.h>

#include "config.hh"


class TrackDataManager
{
  public:
    TrackDataManager();
    ~TrackDataManager();
    
    inline TChain* GetChain() const { return tracksChain; };
    inline Int_t GetNEvents() const { return nEvents; };
    inline Int_t GetProcessedEvents() const { return processedEvents; };

    void ProcessAndFilterEvent();
    Bool_t SetTrueDecayData(Int_t eventID);
    
    // Data (Linked to input TTree branches)
    Int_t   truedecay_EventID;
    Float_t trueDecay_posX,
            trueDecay_posY,
            trueDecay_posZ,
            trueDecay_momX,
            trueDecay_momY,
            trueDecay_momZ,
            trueDecay_time;

    std::vector<Int_t>       *hits_EventID,
                             *hits_TrackID,
                             *hits_ParticleID,
                             *hits_DetectorID;
    std::vector<std::string> *hits_DetName;
    std::vector<Float_t>     *hits_posX,
                             *hits_posY,
                             *hits_posZ,
                             *hits_time;

    // Data (for analysis)
    Float_t                            trueDecayTime;
    TVector3                           trueDecayPos,
                                       trueDecayMom;
    Double_t                           trueDecayThetaEDM;
    std::vector<std::vector<Double_t>> hitsCoordinates;
    std::vector<Int_t>                 hitsCylinderID;
    
    // Fitted Data
    TVector3 fittedDecayPos,
             fittedDecayMom;
    Double_t fittedDecayThetaEDM;

    //Double_t trueMomentum;
    //Double_t polarAngle;
    //Double_t azimuthalAngle;
    //Double_t theThetaAngle;
    //Double_t spinAngle;
    //Double_t emissionAngle;
    //TVector3* fOrigin;
    //std::vector<std::vector<Double_t>>* hitsCoordinates;
    //std::vector<std::vector<Double_t>>* trackCoordinates;
    //std::vector<Int_t>* cylinderID;

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
    TChain *tracksChain;
    TChain *decayChain;
    Int_t nEvents,
          processedEvents;

    // Map for synchro: EventID -> EntryIndex in decayChain
    std::map<Int_t, Long64_t> decayIndexMap;

    // Function to construct the map
    void BuildDecayIndex();

    // Useful methods for filter
    std::set<std::string> excludeStrings = {"EntranceTrigger", "VerticalScint","ScintDownEndCap","ScintUpEndCap","EndBar"};

    void FilterVectors(std::vector<std::string>& strVec, 
                       std::vector<std::vector<float>>& fltVecs, 
                       std::vector<std::vector<int>>& intVecs);

    int FindFirstTrack(const std::vector<int>& trackIDs, 
                       const std::vector<int>& particleIDs, 
                       int targetParticleID);

    void FilterParticleAndTrack(std::vector<std::string>& strVec, 
                                std::vector<std::vector<float>>& fltVecs, 
                                std::vector<std::vector<int>>& intVecs, 
                                int trackID, int particleID);
};
    
#endif  // TRACKDATAMANAGER_HH
