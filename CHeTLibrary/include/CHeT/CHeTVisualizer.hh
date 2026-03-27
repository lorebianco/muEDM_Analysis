#ifndef CHETVISUALIZER_HH
#define CHETVISUALIZER_HH

// --- External Dependencies ---
#include <algorithm>
#include <cmath>
#include <functional>
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

namespace CHeT
{

/**
 * @namespace CHeT::Vis
 * @brief CHeT Visualization namespace.
 * Contains tools for 2D and 3D visualization of the detector, tracks, and hits
 * using ROOT.
 */
namespace Vis
{

// --- Visualization Data Structures ---

/**
 * @brief Structure representing a 2D point (e.g., for Hough space or transverse
 * projection).
 */
struct VisPoint2D
{
    double x, y; ///< Coordinates
    int color; ///< ROOT Color index
    int markerStyle; ///< ROOT Marker style (e.g., 20 = Full Circle)
    double size; ///< Marker size

    VisPoint2D() = default;

    VisPoint2D(double _x, double _y, int _col = kBlue, int _mst = 1, double _sz = 1.0)
        : x(_x)
        , y(_y)
        , color(_col)
        , markerStyle(_mst)
        , size(_sz)
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
    bool isLocalFrame; ///< If true, coordinates are in Detector Local Frame

    VisPoint3D() = default;

    VisPoint3D(double _x, double _y, double _z, int _col = kRed, int _mst = 20, double _sz = 1.0,
        bool _isLocal = false)
        : x(_x)
        , y(_y)
        , z(_z)
        , color(_col)
        , markerStyle(_mst)
        , size(_sz)
        , isLocalFrame(_isLocal)
    {
    }
};

/**
 * @brief Structure representing a reconstructed 3D generic track for visualization.
 */
struct VisGenericTrack
{
    std::vector<VisPoint3D> points; ///< Sequence of 3D points
    int color; ///< ROOT Color index
    int width; ///< Line width
    int style; ///< Line style
    bool isLocalFrame; ///< If true, coordinates are in Detector Local Frame

    VisGenericTrack() = default;

    VisGenericTrack(
        std::vector<VisPoint3D> _pts, int _c = 1, int _w = 2, int _s = 1, bool _loc = false)
        : points(std::move(_pts))
        , color(_c)
        , width(_w)
        , style(_s)
        , isLocalFrame(_loc)
    {
    }
};

/**
 * @brief Structure representing a reconstructed 3D line track for visualization.
 */
struct VisLineTrack
{
    double x0, y0, z0; ///< Origin point
    double ux, uy, uz; ///< Direction unit vector
    int color; ///< ROOT Color index
    int width; ///< Line width
    int style; ///< Line style
    bool isLocalFrame; ///< If true, coordinates are in Detector Local Frame

    /**
     * @brief Construct a new VisLineTrack.
     * Automatically normalizes the direction vector (dx, dy, dz).
     */
    VisLineTrack() = default;

    VisLineTrack(double _x, double _y, double _z, double _dx, double _dy, double _dz, int _col = 54,
        /* kYellow fallback */ int _w = 2, int _s = 1, bool _isLocal = false)
        : x0(_x)
        , y0(_y)
        , z0(_z)
        , color(_col)
        , width(_w)
        , style(_s)
        , isLocalFrame(_isLocal)
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
 * @brief Structure representing a reconstructed 3D helix track for visualization.
 */
struct VisHelixTrack
{
    double cx, cy; ///< Helix center in XY plane
    double radius; ///< Radius
    double z0; ///< Initial Z position (at t=0)
    double dz_dt; ///< Pitch: Z advancement per radian
    double t_min = 0; ///< Minimum angular parameter range
    double t_max = 6.28; ///< Maximum angular parameter range (default: 2*PI)
    int color; ///< ROOT Color index
    int width; ///< Line width
    int style; ///< Line style
    bool isLocalFrame; ///< If true, coordinates are in Detector Local Frame

    VisHelixTrack() = default;

    VisHelixTrack(double _cx, double _cy, double _r, double _z0, double _dz, double _tmin,
        double _tmax, int _c = 1, int _w = 2, int _s = 1, bool _loc = false)
        : cx(_cx)
        , cy(_cy)
        , radius(_r)
        , z0(_z0)
        , dz_dt(_dz)
        , t_min(_tmin)
        , t_max(_tmax)
        , color(_c)
        , width(_w)
        , style(_s)
        , isLocalFrame(_loc)
    {
    }
};

// --- Core Functions ---

void Draw2DCore(const std::vector<int> &bundle_ids, const std::vector<VisPoint2D> &extraPoints,
    const std::function<void(TMultiGraph *, TMultiGraph *, TMultiGraph *, TMultiGraph *)>
        &trackDrawer,
    bool wrap_phi = true);

void Draw3DCore(const std::vector<int> &hit_ids, const std::vector<VisPoint3D> &points,
    bool drawSkeleton, const std::function<void(double[6])> &trackDrawer);

// --- Rendering Functions ---
void RenderTrack2D(const VisLineTrack &tr, TMultiGraph *mg_xy);
void RenderTrack2D(const VisHelixTrack &tr, TMultiGraph *mg_xy);
void RenderTrack2D(const VisGenericTrack &tr, TMultiGraph *mg_xy);

void RenderTrackZX(const VisLineTrack &tr, TMultiGraph *mg_zx);
void RenderTrackZX(const VisHelixTrack &tr, TMultiGraph *mg_zx);
void RenderTrackZX(const VisGenericTrack &tr, TMultiGraph *mg_zx);

void RenderTrackZY(const VisLineTrack &tr, TMultiGraph *mg_zy);
void RenderTrackZY(const VisHelixTrack &tr, TMultiGraph *mg_zy);
void RenderTrackZY(const VisGenericTrack &tr, TMultiGraph *mg_zy);

void RenderTrackPhiZ(const VisLineTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi = true);
void RenderTrackPhiZ(const VisHelixTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi = true);
void RenderTrackPhiZ(const VisGenericTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi = true);

void RenderTrack3D(const VisLineTrack &tr, double bounds[6]);
void RenderTrack3D(const VisHelixTrack &tr, double bounds[6]);
void RenderTrack3D(const VisGenericTrack &tr, double bounds[6]);

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
template <typename TrackContainer = std::vector<VisGenericTrack>>
void Draw2D(const std::vector<int> &bundle_ids, const TrackContainer &tracks = {},
    const std::vector<VisPoint2D> &extraPoints = {}, bool wrap_phi = true)
{
    // Invoke core function with a lambda iterating over the generic container
    Draw2DCore(
        bundle_ids, extraPoints,
        [&](TMultiGraph *mg_xy, TMultiGraph *mg_zx, TMultiGraph *mg_zy, TMultiGraph *mg_phiz)
        {
            for(const auto &tr : tracks)
            {
                RenderTrack2D(tr, mg_xy); // Compiler chooses the correct overload!
                RenderTrackZX(tr, mg_zx);
                RenderTrackZY(tr, mg_zy);
                RenderTrackPhiZ(tr, mg_phiz, wrap_phi);
            }
        },
        wrap_phi);
}

/**
 * @brief Draws the full 3D detector view.
 *
 * @param hit_ids Vector of global bundle IDs that fired (drawn as solid lines).
 * @param tracks Optional list of tracks to draw.
 * @param points Optional list of 3D points to draw.
 * @param drawSkeleton If true, draws all inactive fibers as faint transparent
 * lines (computationally heavy).
 */
template <typename TrackContainer = std::vector<VisGenericTrack>>
void Draw3D(const std::vector<int> &hit_ids, const TrackContainer &tracks = {},
    const std::vector<VisPoint3D> &points = {}, bool drawSkeleton = true)
{
    Draw3DCore(hit_ids, points, drawSkeleton,
        [&](double bounds[6])
        {
            for(const auto &tr : tracks)
            {
                RenderTrack3D(tr, bounds);
            }
        });
}

} // namespace Vis
} // namespace CHeT

#endif // CHETVISUALIZER_HH
