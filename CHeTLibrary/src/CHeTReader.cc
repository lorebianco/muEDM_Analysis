#include "CHeT/CHeTReader.hh"

using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

// Aliases for convenience
using RVecUI = RVec<unsigned int>;
using RVecUS = RVec<unsigned short>;
using RVecUC = RVec<unsigned char>;

namespace CHeT
{
namespace Data
{

Reader::Reader(const std::string &filename, const std::string &treeName)
    : fFilename(filename)
    , fTreeName(treeName)
    , fDF(treeName, filename)
    , fHeadNode(fDF)
{
    // Disable implicit MT if necessary, or manage it at higher level
}

void Reader::SetCuts(
    double toaMin, double toaMax, unsigned int totMin, unsigned int totMax)
{
    fToaMin = toaMin;
    fToaMax = toaMax;
    fTotMin = totMin;
    fTotMax = totMax;
}

void Reader::SetSingleEvent(long eventID)
{
    fHeadNode = fDF.Range(eventID, eventID + 1);
}

void Reader::SetEnabledBoards(const std::vector<int> &boards)
{
    fEnabledBoards = boards;
}

void Reader::SetEnabledCylinders(const std::vector<int> &cylinders)
{
    fEnabledCylinders = cylinders;
}

void Reader::SetEnabledLayers(const std::vector<int> &layers)
{
    fEnabledLayers = layers;
}

void Reader::SetEnabledGeometries(
    const std::vector<std::pair<int, int>> &geometries)
{
    fEnabledGeometries = geometries;
}

ROOT::RDF::RNode Reader::GetRaw()
{
    return fHeadNode;
}

ROOT::RDF::RNode Reader::GetEstimators()
{
    // Capture parameters for lambdas
    double tMin = fToaMin;
    double tMax = fToaMax;
    unsigned int totMin = fTotMin;
    unsigned int totMax = fTotMax;

    // Capture filters
    auto enabledBoards = fEnabledBoards;
    auto enabledCylinders = fEnabledCylinders;
    auto enabledLayers = fEnabledLayers;
    auto enabledGeometries = fEnabledGeometries;

    // 1. ToA Corrections (Alignment w.r.t. FD00)
    // FD00 is the reference: only LSB -> ns conversion
    auto df_corr = fHeadNode.Define("FD00_corrToA",
        [](const RVecUI &toa) { return RVecD(toa) / 2.0; }, { "FD00_ToA" });

    // Other boards (FD01..FD03) are aligned using fine timestamp
    for(const string &id : { "01", "02", "03" })
    {
        df_corr = df_corr.Define("FD" + id + "_corrToA",
            [](const RVecUI &toa, Double_t t, Double_t t0)
            {
                // Correction: ToA/2 + fine timestamp difference (ps -> ns?)
                // Note: original code had * 1e3, check units!
                return RVecD(toa) / 2.0 + (t - t0) * 1e3;
            },
            { "FD" + id + "_ToA", "FD" + id + "_fine_tstamp",
                "FD00_fine_tstamp" });
    }

    // 2. Loop Board -> Selection -> Bundle Mapping
    auto df_geo = df_corr;
    vector<string> toa_cols, tot_cols, cyl_cols, lay_cols, bund_cols;

    // Map string ID -> numeric ID for geometry
    map<string, int> board_id_map
        = { { "00", 0 }, { "01", 1 }, { "02", 2 }, { "03", 3 } };

    for(const auto &[id_str, id_num] : board_id_map)
    {
        string id = id_str; // "00", "01"...
        int b_id = id_num;

        string toa_n = "FD" + id + "_ToA_sel";
        string tot_n = "FD" + id + "_ToT_sel";
        string cyl_n = "FD" + id + "_cyl_sel";
        string lay_n = "FD" + id + "_lay_sel";
        string bund_n = "FD" + id + "_bund_sel";

        toa_cols.push_back(toa_n);
        tot_cols.push_back(tot_n);
        cyl_cols.push_back(cyl_n);
        lay_cols.push_back(lay_n);
        bund_cols.push_back(bund_n);

        // Check Board Filtering
        bool isBoardEnabled = true;
        if(!enabledBoards.empty())
        {
            if(std::find(enabledBoards.begin(), enabledBoards.end(), b_id)
                == enabledBoards.end())
            {
                isBoardEnabled = false;
            }
        }

        if(!isBoardEnabled)
        {
            // Define empty columns if board is disabled
            // We use specific types to match the expected types in Concatenate
            df_geo = df_geo.Define(toa_n, []() { return RVecD {}; })
                         .Define(tot_n, []() { return RVecUS {}; })
                         .Define(cyl_n, []() { return RVecI {}; })
                         .Define(lay_n, []() { return RVecI {}; })
                         .Define(bund_n, []() { return RVecI {}; });
            continue;
        }

        // If board is enabled, we first calculate everything into temporary
        // columns, then apply cylinder/layer filtering.
        string toa_tmp = toa_n + "_tmp";
        string tot_tmp = tot_n + "_tmp";
        string cyl_tmp = cyl_n + "_tmp";
        string lay_tmp = lay_n + "_tmp";
        string bund_tmp = bund_n + "_tmp";

        // Time cuts on ToA and ToT
        df_geo = df_geo.Define(toa_tmp,
            [=](const RVecD &toa, const RVecUS &tot)
            {
                return toa[(toa > tMin) && (toa < tMax) && (tot > totMin)
                    && (tot < totMax)];
            },
            { "FD" + id + "_corrToA", "FD" + id + "_ToT" });

        df_geo = df_geo.Define(tot_tmp,
            [=](const RVecD &toa, const RVecUS &tot)
            {
                return tot[(toa > tMin) && (toa < tMax) && (tot > totMin)
                    && (tot < totMax)];
            },
            { "FD" + id + "_corrToA", "FD" + id + "_ToT" });

        // Bundle ID Calculation (Geometry)
        df_geo = df_geo.Define(bund_tmp,
            [=](const RVecD &toa, const RVecUS &tot, const RVecUC &chans)
            {
                RVecI bunds;
                auto mask = (toa > tMin) && (toa < tMax) && (tot > totMin)
                    && (tot < totMax);
                auto sel_chans = chans[mask];
                for(auto c : sel_chans)
                {
                    int gid = CHeT::Config::GetGlobalBundleId(b_id, (int)c);
                    bunds.push_back(gid);
                }
                return bunds;
            },
            { "FD" + id + "_corrToA", "FD" + id + "_ToT",
                "FD" + id + "_channel" });

        // Derive Cylinder and Layer from Global IDs
        df_geo = df_geo.Define(cyl_tmp,
            [](const RVecI &bunds)
            {
                RVecI cyls;
                for(auto b : bunds)
                    cyls.push_back((b >= 0)
                            ? CHeT::Config::GetFiberProp(b).cylinderId
                            : -1);
                return cyls;
            },
            { bund_tmp });

        df_geo = df_geo.Define(lay_tmp,
            [](const RVecI &bunds)
            {
                RVecI lays;
                for(auto b : bunds)
                    lays.push_back(
                        (b >= 0) ? CHeT::Config::GetFiberProp(b).layerId : -1);
                return lays;
            },
            { bund_tmp });

        // Create Filter Mask based on Cylinder and Layer
        string mask_n = "FD" + id + "_geo_mask";
        df_geo = df_geo.Define(mask_n,
            [enabledCylinders, enabledLayers, enabledGeometries](
                const RVecI &cyls, const RVecI &lays)
            {
                // By default, keep everything (mask=1)
                RVec<int> mask(cyls.size(), 1);

                // If no filters are set, return all 1s
                if(enabledCylinders.empty() && enabledLayers.empty()
                    && enabledGeometries.empty())
                {
                    return mask;
                }

                for(size_t i = 0; i < cyls.size(); ++i)
                {
                    // Check Cylinder
                    if(!enabledCylinders.empty())
                    {
                        bool found = false;
                        for(auto c : enabledCylinders)
                            if(c == cyls[i])
                            {
                                found = true;
                                break;
                            }
                        if(!found)
                        {
                            mask[i] = 0;
                            continue;
                        } // Optimization: if fail cyl, fail hit
                    }

                    // Check Layer
                    if(!enabledLayers.empty())
                    {
                        bool found = false;
                        for(auto l : enabledLayers)
                            if(l == lays[i])
                            {
                                found = true;
                                break;
                            }
                        if(!found)
                            mask[i] = 0;
                    }

                    // Check Specific Geometries (Cyl, Lay)
                    if(!enabledGeometries.empty())
                    {
                        bool found = false;
                        for(const auto &geo : enabledGeometries)
                        {
                            if(geo.first == cyls[i] && geo.second == lays[i])
                            {
                                found = true;
                                break;
                            }
                        }
                        if(!found)
                            mask[i] = 0;
                    }
                }
                return mask;
            },
            { cyl_tmp, lay_tmp });

        // Apply Mask to final columns
        // We use explicit lambdas to ensure type correctness
        df_geo = df_geo
                     .Define(toa_n,
                         [](const RVecD &v, const RVecI &m) { return v[m]; },
                         { toa_tmp, mask_n })
                     .Define(tot_n,
                         [](const RVecUS &v, const RVecI &m) { return v[m]; },
                         { tot_tmp, mask_n })
                     .Define(cyl_n,
                         [](const RVecI &v, const RVecI &m) { return v[m]; },
                         { cyl_tmp, mask_n })
                     .Define(lay_n,
                         [](const RVecI &v, const RVecI &m) { return v[m]; },
                         { lay_tmp, mask_n })
                     .Define(bund_n,
                         [](const RVecI &v, const RVecI &m) { return v[m]; },
                         { bund_tmp, mask_n });
    }

    // 3. Global Concatenation (All Boards -> Single Event Vector)
    auto df_final
        = df_geo
              .Define(
                  "All_Cyl",
                  [](const RVecI &c0, const RVecI &c1, const RVecI &c2,
                      const RVecI &c3) {
                      return Concatenate(
                          c0, Concatenate(c1, Concatenate(c2, c3)));
                  },
                  cyl_cols)
              .Define(
                  "All_Lay",
                  [](const RVecI &l0, const RVecI &l1, const RVecI &l2,
                      const RVecI &l3) {
                      return Concatenate(
                          l0, Concatenate(l1, Concatenate(l2, l3)));
                  },
                  lay_cols)
              .Define(
                  "All_Bundle",
                  [](const RVecI &b0, const RVecI &b1, const RVecI &b2,
                      const RVecI &b3) {
                      return Concatenate(
                          b0, Concatenate(b1, Concatenate(b2, b3)));
                  },
                  bund_cols)

              // Counts and Event-Level variables
              .Define("nHits_Total", "(double)All_Cyl.size()")
              .Define("nHits_Cyl0", "Sum(All_Cyl == 0)")
              .Define("nHits_Cyl1", "Sum(All_Cyl == 1)")
              .Define("nHits_Cyl0_Inner", "Sum(All_Cyl == 0 && All_Lay == 0)")
              .Define("nHits_Cyl0_Outer", "Sum(All_Cyl == 0 && All_Lay == 1)")
              .Define("nHits_Cyl1_Inner", "Sum(All_Cyl == 1 && All_Lay == 0)")
              .Define("nHits_Cyl1_Outer", "Sum(All_Cyl == 1 && All_Lay == 1)")
              .Define(
                  "FirstToA",
                  [](const RVecD &t0, const RVecD &t1, const RVecD &t2,
                      const RVecD &t3)
                  {
                      auto all = Concatenate(
                          t0, Concatenate(t1, Concatenate(t2, t3)));
                      return all.empty()
                          ? std::numeric_limits<Double_t>::infinity()
                          : Min(all);
                  },
                  toa_cols)
              .Define(
                  "SumToT",
                  [](const RVecUS &t0, const RVecUS &t1, const RVecUS &t2,
                      const RVecUS &t3)
                  {
                      auto all = Concatenate(
                          t0, Concatenate(t1, Concatenate(t2, t3)));
                      return (Double_t)Sum(all);
                  },
                  tot_cols);

    return df_final;
}

} // namespace Data
} // namespace CHeT
