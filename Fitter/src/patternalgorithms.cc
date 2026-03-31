#include "patternalgorithms.hh"

using namespace std;



pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SortVectorZ(vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Sort vectors according to z(third) value

    vector<Double_t> z_of_vectors;
    for(auto v : vectors)
        z_of_vectors.push_back(abs(v.at(2)));

    // initialize original index locations
    vector<size_t> idx(z_of_vectors.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using stable_sort instead of sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&z_of_vectors](size_t i1, size_t i2) {return z_of_vectors[i1] < z_of_vectors[i2];});

    vector<vector<Double_t>> copy_vector;
    for(auto i : idx)
        copy_vector.push_back(vectors.at(i));

    vector<Int_t> copy_cylinders;
    for(auto i : idx)
        copy_cylinders.push_back(cylinders.at(i));

    return {copy_vector, copy_cylinders};
}



pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::ShuffleVectorZ(vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Shuffle the vectors and cylinders randomly
    
    // Create a random engine and a distribution
    random_device rd;
    mt19937 g(rd());

    // Combine vectors and cylinders into a single vector of pairs
    vector<pair<vector<Double_t>, Int_t>> combined;
    for (size_t i = 0; i < vectors.size(); ++i) {
        combined.push_back({vectors[i], cylinders[i]});
    }

    // Shuffle the combined vector randomly
    shuffle(combined.begin(), combined.end(), g);

    // Extract the shuffled vectors and cylinders back
    vector<vector<Double_t>> shuffled_vectors;
    vector<Int_t> shuffled_cylinders;
    for (const auto& p : combined) {
        shuffled_vectors.push_back(p.first);
        shuffled_cylinders.push_back(p.second);
    }

    return {shuffled_vectors, shuffled_cylinders};
}



Int_t PTTALG::CountTurns(const vector<vector<Double_t>> hitsCoordinates) 
{
    if(hitsCoordinates.size() < 3) 
        return 0; // Servono almeno 3 punti per trovare un massimo o minimo

    Int_t nTurns = 0;
    Bool_t foundMax = false, foundMin = false;

    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se abbiamo sia un massimo che un minimo -> un giro completato
        if(foundMax && foundMin)
        {
            nTurns++;
            foundMax = false;
            foundMin = false;
        }
    }

    return nTurns;
}



pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SelectTurn(Float_t turnID, const vector<vector<Double_t>>& hitsCoordinates, const vector<Int_t>& cylinders)
{
    vector<vector<Double_t>> turnHits;
    vector<Int_t> turnCylinders;

    if (hitsCoordinates.size() < 3)
        return {turnHits, turnCylinders}; // Troppi pochi punti per definire un giro

    Int_t totalTurns = CountTurns(hitsCoordinates);
    if (turnID > totalTurns)
        return {hitsCoordinates, cylinders}; // Se voglio più giri di quelli presenti, prendo tutta la traccia

    Double_t nHalfTurns = 0.0;
    turnHits.push_back(hitsCoordinates.front()); // Includi il primo punto
    turnCylinders.push_back(cylinders.front());

    for (size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Identificazione di un massimo o minimo locale
        if ((yCurr > yPrev && yCurr > yNext) || (yCurr < yPrev && yCurr < yNext))
        {
            nHalfTurns += 1.0; // Ora conto direttamente i mezzi giri
        }

        // Se il numero di **giri completi** supera `turnID`, interrompo
        if (nHalfTurns / 2.0 >= turnID)
        {
            break;
        }

        turnHits.push_back(hitsCoordinates[i]);
        turnCylinders.push_back(cylinders[i]);
    }

    // Includi sempre l'ultimo punto del semigiro
    turnHits.push_back(hitsCoordinates[turnHits.size()]);
    turnCylinders.push_back(cylinders[turnCylinders.size()]);

    return {turnHits, turnCylinders};
}



vector<Int_t> PTTALG::SplitTurns(const vector<vector<Double_t>>& hitsCoordinates) 
{
    vector<Int_t> turnIndices;
    if(hitsCoordinates.size() < 3) 
        return turnIndices;
    
    turnIndices.push_back(0); // Il primo indice è sempre 0
    
    Bool_t foundMax = false, foundMin = false;
    
    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se troviamo un massimo e poi un minimo (o viceversa), aggiungiamo l'indice
        if(foundMax && foundMin)
        {
            turnIndices.push_back(i);
            foundMax = false;
            foundMin = false;
        }
    }
    
    return turnIndices;
}



Int_t PTTALG::CountCylinders(const vector<Int_t>& cylinders)
{
    set<Int_t> uniqueValues(cylinders.begin(), cylinders.end());
    return uniqueValues.size();
}



pair<vector<vector<Double_t>>, vector<Int_t>> PTTALG::SelectCylinders(vector<Int_t> targetCylinders, const vector<vector<Double_t>>& hitsCoordinates, const vector<Int_t>& cylinders)
{
    vector<vector<Double_t>> selectedHits;
    vector<Int_t> selectedCylinders;

    if(hitsCoordinates.size() != cylinders.size())
        return {selectedHits, selectedCylinders};

    for(size_t i = 0; i < hitsCoordinates.size(); i++)
    {
        Int_t currentID = cylinders[i];
        Bool_t isTarget = false;

        for(size_t j = 0; j < targetCylinders.size(); j++)
        {
            if(targetCylinders[j] == currentID)
            {
                isTarget = true;
                break;
            }
        }

        if(isTarget)
        {
            selectedHits.push_back(hitsCoordinates[i]);
            selectedCylinders.push_back(cylinders[i]);
        }
    }

    return {selectedHits, selectedCylinders};
}



pair<TVector3, TVector3> PTTALG::SmearSeed(TVector3 pos, TVector3 mom)
{
    TVector3 smearPos(pos);
    TVector3 smearMom(mom);

    smearPos.SetX(gRandom->Gaus(pos.X(), ANS::sigmaSeedPos));
    smearPos.SetY(gRandom->Gaus(pos.Y(), ANS::sigmaSeedPos));
    smearPos.SetZ(gRandom->Gaus(pos.Z(), ANS::sigmaSeedPos));

    smearMom.SetPhi(gRandom->Gaus(mom.Phi(), ANS::sigmaSeedMomDir));
    smearMom.SetTheta(gRandom->Gaus(mom.Theta(), ANS::sigmaSeedMomDir));
    smearMom.SetMag(gRandom->Gaus(mom.Mag(), ANS::sigmaSeedMomMag*mom.Mag()));

    return {smearPos, smearMom};
}
