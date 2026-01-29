#ifndef CHETGLOBALSETTINGS_HH
#define CHETGLOBALSETTINGS_HH

#include <algorithm> // for std::find if needed
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/**
 * @namespace GS
 * @brief Global Settings namespace containing configuration structures,
 * physical constants, and utility functions for the detector geometry.
 */
namespace GS
{

// --- Configuration Structures ---

/**
 * @brief Configuration for a single fiber layer.
 */
struct LayerConfig
{
    int nBundles; ///< Number of bundles in the layer
    double radius; ///< Radius in [mm]
    double phiOffset; ///< Angular offset in [rad]
    int direction; ///< Winding direction: +1 (CCW) or -1 (CW)
    int color; ///< ROOT color index for display
};

/**
 * @brief Configuration for a full cylinder (Inner + Outer layers).
 */
struct CylinderConfig
{
    int id; ///< Cylinder ID
    LayerConfig inner; ///< Configuration for the inner layer
    LayerConfig outer; ///< Configuration for the outer layer
};

/**
 * @brief Properties of a specific reconstructed fiber.
 */
struct FiberProp
{
    int cylinderId;
    int layerId; ///< 0 for inner, 1 for outer
    double r; ///< Radius [mm]
    double phi0; ///< Initial angle [rad]
    int dir; ///< Winding direction
    int color; ///< ROOT color index
};

/**
 * @brief Represents a geometric intersection between two bundles.
 */
struct BundlesIntersection
{
    double z; ///< Longitudinal position [mm]
    double x; ///< Transverse position X [mm]
    double y; ///< Transverse position Y [mm]
    int cylinderId; ///< ID of the cylinder where intersection occurs
    double widthZ; ///< Longitudinal uncertainty due to fiber thickness [mm]
};

// --- Global Physical Constants ---

// Kept as constexpr in header for compile-time optimization
constexpr double L_HALF = 150.0; ///< Half-length of the detector [mm]
constexpr double BUNDLE_WIDTH = 2.0; ///< Physical width of a bundle [mm]

// --- Experimental Offsets ---
constexpr double OFFSET_EXP = 40.0 * (M_PI / 180.0);
constexpr double DELTA1 = 1.267;
constexpr double DELTA2 = 1.420 + (18.0 * M_PI / 180.0);

// --- Function Declarations ---

/**
 * @brief Returns the global offset (bundle count) for a specific board.
 * @param board_id The hardware board ID.
 * @return The offset in bundle index, or -1 if invalid.
 */
int GetBoardGlobalOffset(int board_id);

/**
 * @brief Retrieves the configuration for all cylinders.
 * @return A constant reference to the vector of cylinder configurations.
 */
const std::vector<CylinderConfig> &GetCylinders();

/**
 * @brief Wraps an angle into the [0, 2*PI) range.
 * @param angle Input angle in radians.
 * @return Wrapped angle.
 */
double wrap0_2pi(double angle);

/**
 * @brief Retrieves fiber properties given a global bundle ID.
 * @param b_id Global bundle ID.
 * @return FiberProp structure containing geometry info.
 */
FiberProp GetFiberProp(int b_id);

/**
 * @brief Converts Hardware Board/Channel to Global Bundle ID.
 * @param board_id Hardware board ID.
 * @param channel_id Hardware channel ID.
 * @return Global Bundle ID, or -1 if mapping not found.
 */
int GetGlobalBundleId(int board_id, int channel_id);

/**
 * @brief Finds 3D intersections between a list of hit bundles.
 * @param hit_ids Vector of global bundle IDs that fired.
 * @return Vector of intersection points.
 */
std::vector<BundlesIntersection> FindIntersections(
    const std::vector<int> &hit_ids);

// --- Debug / Helper Functions ---

/**
 * @brief Helper: Converts Geometry (Cyl, Layer, Index) to Global ID.
 */
int GetGlobalIdFromGeometry(int cyl_id, int layer_id, int layer_idx);

/**
 * @brief Helper: Prints mapping details for a given Global ID to stdout.
 */
void PrintBundleMapping(int global_id);

/**
 * @brief Interactive CLI menu to explore the mapping.
 */
void MapExplorer();

} // namespace GS

#endif // CHETGLOBALSETTINGS_HH
