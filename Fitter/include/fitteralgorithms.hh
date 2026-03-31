#ifndef FITTERALGORITHMS_HH
#define FITTERALGORITHMS_HH

#include <AbsBFieldDer.h>
#include <AbsMeasurement.h>
#include <ConstField.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <Exception.h>
#include <FieldManager.h>
#include <Fit/Fitter.h>
#include <KalmanFitter.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <MeasurementFactory.h>
#include <MeasurementProducer.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mySpacepointDetectorHit.h>
#include <mySpacepointMeasurement.h>
#include <numeric>
#include <random>
#include <unistd.h>
#include <vector>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <RKTrackRep.h>
#include <ROOT/RVec.hxx>
#include <TApplication.h>
#include <TArc.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TDecompChol.h>
#include <TDecompSVD.h>
#include <TEfficiency.h>
#include <TEllipse.h>
#include <TEveManager.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGeoManager.h>
#include <TGeoMaterialInterface.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSymEigen.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVector3.h>

#include "CHeT/CHeTGlobalSettings.hh"
#include "CHeT/CHeTVisualizer.hh"
#include "auxiliaryalgorithms.hh"
#include "config.hh"
#include "globalsettings.hh"
#include "patternalgorithms.hh"
#include "trackdatamanager.hh"

namespace FITALG
{
struct RecoTrack
{
    double x0;
    double z0;
    double sx;
    double sz;
    double chi2;
    bool converged;
};

struct FitOutput
{
    RecoTrack track;
    std::vector<CHeT::Vis::VisPoint3D> fittedPoints;
};

struct FitData
{
    std::vector<CHeT::Config::FiberProp> props;
};

double Track3DNeg2LogL(const double *par, const FitData &data, bool usePrior);
FitOutput Do3DFit(const std::vector<int> &hit_ids, bool usePrior = false);

void CosmicFitter();
void SpacepointFitter();
void HelixFitter();

std::array<Double_t, 6> HelixPrefitter(
    const std::vector<std::vector<Double_t>> &hitsCoordinates, const std::vector<Int_t> &cylinders);

void AddFakeHitFromHelix(genfit::TrackCand &trackCand, Int_t hitIndex, Double_t sortingParameter,
    Double_t s, Double_t xC, Double_t yC, Double_t R, Double_t z0, Double_t phi0,
    Double_t tanLambda, TClonesArray &chetHitArray, std::vector<TVector3> &virtualCoordinates,
    Double_t sigmaBig);
};

#endif // FITTERALGORITHMS_HH
