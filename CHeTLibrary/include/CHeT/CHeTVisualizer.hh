#ifndef CHETVISUALIZER_HH
#define CHETVISUALIZER_HH

// --- External Dependencies ---
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// --- ROOT Headers ---
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TStyle.h>

// --- Project Dependencies ---
#include "CHeT/CHeTGlobalSettings.hh"

/**
 * @namespace CHeT::Vis
 * @brief CHeT Visualization namespace.
 * Contains tools for 2D and 3D visualization of the detector, tracks, and hits
 * using ROOT.
 */
namespace CHeT
{
namespace Vis
{

// --- Visualization Data Structures ---

/**
 * @brief Structure representing a reconstructed 3D track for visualization.
 */
struct VisLineTrack
{
    double x0, y0, z0; ///< Origin point
    double ux, uy, uz; ///< Direction unit vector
    int color; ///< ROOT Color index
    int width; ///< Line width
    int style; ///< Line style

    /**
     * @brief Construct a new VisLineTrack.
     * Automatically normalizes the direction vector (dx, dy, dz).
     */
    VisLineTrack(double _x, double _y, double _z, double _dx, double _dy,
        double _dz, int _col = kYellow, int _w = 2, int _s = 1)
        : x0(_x)
        , y0(_y)
        , z0(_z)
        , color(_col)
        , width(_w)
        , style(_s)
    {
        double norm = std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
        if(norm > 0)
        {
            ux = _dx / norm;
            uy = _dy / norm;
            uz = _dz / norm;
        }
        else
        {
            ux = 0;
            uy = 0;
            uz = 0;
        }
    }
};

/**
 * @brief Structure representing a 2D point (e.g., for Hough space or transverse
 * projection).
 */
struct VisPoint2D
{
    double x, y; ///< Coordinates
    int color; ///< ROOT Color index
    int markerStyle; ///< ROOT Marker style (e.g., 20 = Full Circle)

    VisPoint2D(double _x, double _y, int _col = kBlue, int _mst = 1)
        : x(_x)
        , y(_y)
        , color(_col)
        , markerStyle(_mst)
    {
    }
};

/**
 * @brief Structure representing a 3D point (e.g., space points, clusters).
 */
struct VisPoint3D
{
    double x, y, z; ///< Coordinates
    int color; ///< ROOT Color index
    int markerStyle; ///< ROOT Marker style
    double size; ///< Marker size

    VisPoint3D(double _x, double _y, double _z, int _col = kRed, int _mst = 20,
        double _sz = 1.0)
        : x(_x)
        , y(_y)
        , z(_z)
        , color(_col)
        , markerStyle(_mst)
        , size(_sz)
    {
    }
};

// --- Visualization Functions ---

/**
 * @brief Draws 2D projections: Phi-Z Unrolled Map, ZX/ZY Side views, and XY
 * Transverse view.
 *
 * @param bundle_ids Vector of global bundle IDs that fired.
 * @param inters Vector of intersections (usually from
 * CHeT::Config::FindIntersections).
 * @param tracks Optional list of tracks to overlay.
 * @param extraPoints Optional list of extra 2D points to draw on the XY plane.
 */
void Draw2D(const std::vector<int> &bundle_ids,
    const std::vector<CHeT::Config::BundlesIntersection> &inters,
    const std::vector<VisLineTrack> &tracks = {},
    const std::vector<VisPoint2D> &extraPoints = {});

/**
 * @brief Draws the full 3D detector view.
 *
 * @param hit_ids Vector of global bundle IDs that fired (drawn as solid lines).
 * @param tracks Optional list of tracks to draw.
 * @param points Optional list of 3D points to draw.
 * @param drawSkeleton If true, draws all inactive fibers as faint transparent
 * lines (computationally heavy).
 */
void Draw3D(const std::vector<int> &hit_ids,
    const std::vector<VisLineTrack> &tracks = {},
    const std::vector<VisPoint3D> &points = {}, bool drawSkeleton = true);

} // namespace Vis
} // namespace CHeT

#endif // CHETVISUALIZER_HH
