#ifndef CONFIG_HH
#define CONFIG_HH

#include <Rtypes.h>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "CHeT/CHeTGlobalSettings.hh"
#include "globalsettings.hh"

enum class FitterType
{
    Unknown,
    CosmicFitter,
    SpacepointFitter,
    HelixFitter
};

class Config
{
  public:
    static Config &get()
    {
        static Config instance;
        return instance;
    }

    // I/O files
    std::string inputDir = "../../data/input/toy/";
    std::string inputFileName = "cosmictoy0000.root";
    std::string inputTreeName = "auto";
    std::vector<std::string> inputDataFiles = { inputDir + inputFileName };

    std::string outputDir = "../../data/output/";
    std::string outputFile = outputDir + "fitted_" + inputFileName;

    // Event loop
    std::pair<Long_t, Long_t> rangeLoop = { 0, std::numeric_limits<Long_t>::max() };

    // Geometry and Fields
    std::string geomDir = "../geometry/";
    std::string geomFileName = "chet_sim_geometry.gdml";
    std::string geometryFile = geomDir + geomFileName;

    muEDM::Fields::FieldBType bFieldType = muEDM::Fields::FieldBType::Constant;
    std::string bFieldDir = "../fieldmaps/";
    std::string bFieldMapName = "Field24.5MeV/";
    std::string bFieldPath = bFieldDir + bFieldMapName;

    Double_t bScale = -1.0;
    Double_t bConstZ = -28.9;

    FitterType fitter = FitterType::Unknown;
    Bool_t usePrefitter = false;
    Int_t event = -1;
    Float_t turnID = -1.0;
    Bool_t turnMode = false;
    Bool_t useSmearing = false;
    Bool_t pttrecMode = false;
    Bool_t quietMode = false;
    Bool_t processSingle = false;

    double cutToAMin = 215.0;
    double cutToAMax = 250.0;
    double cutToTMin = 40.0;
    double cutToTMax = 220.0;

    // Geometry Overrides
    bool overrideGeom = false;
    double geom_tx = 0.0, geom_ty = 0.0, geom_tz = 0.0;
    double geom_rx = 0.0, geom_ry = 0.0, geom_rz = 0.0;
    double geom_offset = CHeT::Config::GetOffsetExp();
    std::vector<double> geom_deltas = CHeT::Config::GetDeltas();

    // MC Options
    bool useTrueMCGeom = false;
    bool useTrueMCHits = false;
    bool active_cyls[6] = { true, true, true, true, true, true };

    Config(Config const &) = delete;
    void operator=(Config const &) = delete;

  private:
    Config() { }
};

#endif // CONFIG_HH
