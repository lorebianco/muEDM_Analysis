#ifndef AUXILIARYALGORITHMS_HH
#define AUXILIARYALGORITHMS_HH

#include <iostream>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <numeric>

#include <TMath.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>

#include <Track.h>
#include <RKTrackRep.h>
#include <mySpacepointDetectorHit.h>

#include "config.hh"
#include "globalsettings.hh"


namespace AUXALG
{
    void DrawXYView_hits(TVector3* origin, std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas, std::vector<std::vector<Double_t>> virtualCoordinates = {});

    void DrawYZView_hits(TVector3* origin, std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas, std::vector<std::vector<Double_t>> virtualCoordinates = {});

    void DrawXYZView_hits(std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas, std::vector<std::vector<Double_t>> virtualCoordinates = {});
    
    void DrawXYView_arc(Double_t xC, Double_t yC, Double_t R, const std::vector<std::vector<Double_t>> &hitsCoordinates, TCanvas* canvas, Int_t nTurns, Float_t turnID);

    void DrawXYZView_helixFromHits(Double_t xC, Double_t yC, Double_t R,
        Double_t z0, Double_t phi0,
        Double_t tanLambda,
        std::vector<std::vector<Double_t>> &hitsCoordinates,
        TCanvas *canvas);
    
    void DrawZvsSFit(TGraph *graph, TF1 *fitFunc, TCanvas *canvas);

    std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> GetResults(genfit::Track *fitTrack, genfit::AbsTrackRep *rep, TVector3 truePos, Double_t trueMom, Double_t trueTheta, Double_t truePhi, Bool_t normalizedPulls = true);

    std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> GetResults(TVectorD &fittedState, TMatrixDSym &covFittedState, TVector3 truePos, Double_t trueMom, Double_t trueTheta, Double_t truePhi, Bool_t normalizedPulls = true);

    std::vector<Double_t> SmearMeasurement(Int_t cylID, std::vector<Double_t> hitCoords);

};


#endif  // AUXILIARYALGORITHMS_HH