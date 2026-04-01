#ifndef PATTERNALGORITHMS_HH
#define PATTERNALGORITHMS_HH

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <unistd.h>
#include <utility>
#include <vector>

#include <TROOT.h>
#include <TRandom.h>

#include "config.hh"
#include "globalsettings.hh"

// Struttura per Hough Transform
struct LinearHoughResult
{
    Double_t rho;
    Double_t theta;
    Double_t weight;
};

// Struttura per Hough Circolare (Piano XY)
struct CircularHoughResult
{
    Double_t xc;
    Double_t yc;
    Double_t R;
    Double_t weight;
    std::vector<int> track_hits_idx;
};

struct Point3D
{
    double x, y, z;
};

struct ZHoughResult
{
    double z0;
    double dz_ds;
    double weight;
    std::vector<Point3D> track_hits;
    std::vector<int> track_hits_idx;
    double t_min;
    double t_max;
};

namespace PTTALG
{
std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SortVectorZ(
    std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> ShuffleVectorZ(
    std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

Int_t CountTurns(const std::vector<std::vector<Double_t>> hitsCoordinates);

std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SelectTurn(Float_t turnID,
    const std::vector<std::vector<Double_t>> &hitsCoordinates, const std::vector<Int_t> &cylinders);

std::vector<Int_t> SplitTurns(const std::vector<std::vector<Double_t>> &hitsCoordinates);

Int_t CountCylinders(const std::vector<Int_t> &cylinders);

std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SelectCylinders(
    std::vector<Int_t> targetCylinders, const std::vector<std::vector<Double_t>> &hitsCoordinates,
    const std::vector<Int_t> &cylinders);

std::pair<TVector3, TVector3> SmearSeed(TVector3 pos, TVector3 mom);

std::vector<CircularHoughResult> DoCircularHoughTransform(const std::vector<Int_t> &hit_ids,
    int nCandidates = 2, int combinatorial_threshold = 1000, int n_random_triplets = 10000,
    bool drawGraphs = true);

std::vector<ZHoughResult> DoZHoughTransform(const std::vector<Int_t> &hit_ids, double xc, double yc,
    double R_reco, int nCandidates = 1, bool drawGraphs = true, double tol_Z_fit = 15.0);
};

#endif // PATTERNALGORITHMS_HH
