#ifndef CHETREADER_HH
#define CHETREADER_HH

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "CHeT/CHeTGlobalSettings.hh"

namespace CHeT
{
/**
 * @namespace CHeT::Data
 * @brief Data processing namespace.
 * Contains the Reader class and data structures for analyzing CHeT data.
 */
namespace Data
{

/**
 * @class Reader
 * @brief Reads CHeT data from ROOT files and processes it into high-level
 * objects.
 *
 * The Reader class wraps RDataFrame to provide a convenient interface for:
 * - Loading raw TTree data.
 * - Applying time and energy cuts.
 * - Filtering specific events.
 * - Generating "Estimators" (calibrated and mapped hits).
 */
class Reader
{
  public:
    /**
     * @brief Constructor.
     * @param filename Path to the ROOT file (or glob pattern).
     * @param treeName Name of the tree (default: "Event").
     */
    Reader(const std::string &filename, const std::string &treeName = "Event");

    // Disable copy (RDataFrame has unique pointers)
    Reader(const Reader &) = delete;
    Reader &operator=(const Reader &) = delete;

    /**
     * @brief Set cuts for hit selection.
     * Must be called before GetEstimators().
     * @param toaMin Minimum Time of Arrival [ns]
     * @param toaMax Maximum Time of Arrival [ns]
     * @param totMin Minimum Time over Threshold [LSB]
     * @param totMax Maximum Time over Threshold [LSB]
     */
    void SetCuts(double toaMin, double toaMax, unsigned int totMin, unsigned int totMax);

    /**
     * @brief Restricts the analysis to a single event ID.
     * Replaces any previous range or filter on the event entry.
     * @param eventID The event entry index to process.
     */
    void SetSingleEvent(long eventID);

    /**
     * @brief Restricts the analysis to specific boards.
     * Only hits from these boards will be processed.
     * @param boards Vector of board IDs (0, 1, 2, 3).
     */
    void SetEnabledBoards(const std::vector<int> &boards);

    /**
     * @brief Restricts the analysis to specific cylinders.
     * Only hits belonging to these cylinders will be included in the output.
     * @param cylinders Vector of cylinder IDs.
     */
    void SetEnabledCylinders(const std::vector<int> &cylinders);

    /**
     * @brief Restricts the analysis to specific layers.
     * Only hits belonging to these layers will be included in the output.
     * @param layers Vector of layer IDs.
     */
    void SetEnabledLayers(const std::vector<int> &layers);

    /**
     * @brief Restricts the analysis to specific (cylinder, layer) combinations.
     * Only hits matching one of the provided pairs will be included.
     * @param geometries Vector of pairs {cylinderId, layerId}.
     */
    void SetEnabledGeometries(const std::vector<std::pair<int, int>> &geometries);

    /**
     * @brief Returns the Raw node (the original tree).
     */
    ROOT::RDF::RNode GetRaw();

    /**
     * @brief Returns the node with calculated high-level variables.
     *
     * Performs:
     * 1. Time correction (ToA alignment w.r.t. board 0).
     * 2. Hit filtering based on cuts (ToA/ToT).
     * 3. Geometric mapping (Channel ID -> Global Bundle ID -> Layer/Cylinder).
     * 4. Aggregation of all hits into unique vectors (All_Bundle, All_Lay...).
     */
    ROOT::RDF::RNode GetEstimators();

  private:
    std::string fFilename;
    std::string fTreeName;
    ROOT::RDataFrame fDF;
    ROOT::RDF::RNode fHeadNode; // The active head node (potentially filtered)

    // Cut parameters (default: open)
    double fToaMin = -1e9;
    double fToaMax = 1e9;
    unsigned int fTotMin = 0;
    unsigned int fTotMax = 99999;

    // Filters
    std::vector<int> fEnabledBoards;
    std::vector<int> fEnabledCylinders;
    std::vector<int> fEnabledLayers;
    std::vector<std::pair<int, int>> fEnabledGeometries;
};

} // namespace Data
} // namespace CHeT

#endif
