#ifndef PATTERNALGORITHMS_HH
#define PATTERNALGORITHMS_HH

#include <iostream>
#include <unistd.h>
#include <vector>
#include <utility>
#include <numeric>
#include <cmath>
#include <random>
#include <set>

#include <TROOT.h>
#include <TRandom.h>

#include "config.hh"
#include "globalsettings.hh"

namespace PTTALG
{
    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SortVectorZ(std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> ShuffleVectorZ(std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

    Int_t CountTurns(const std::vector<std::vector<Double_t>> hitsCoordinates);

    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SelectTurn(Float_t turnID, const std::vector<std::vector<Double_t>>& hitsCoordinates, const std::vector<Int_t>& cylinders);

    std::vector<Int_t> SplitTurns(const std::vector<std::vector<Double_t>>& hitsCoordinates);

    Int_t CountCylinders(const std::vector<Int_t>& cylinders);

    std::pair<std::vector<std::vector<Double_t>> , std::vector<Int_t>> SelectCylinders(std::vector<Int_t> targetCylinders, const std::vector<std::vector<Double_t>>& hitsCoordinates, const std::vector<Int_t>& cylinders);

    std::pair<TVector3, TVector3> SmearSeed(TVector3 pos, TVector3 mom);
};


#endif  // PATTERNALGORITHMS_HH