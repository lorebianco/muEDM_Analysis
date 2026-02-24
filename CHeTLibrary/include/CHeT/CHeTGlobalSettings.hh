#ifndef CHETGLOBALSETTINGS_HH
#define CHETGLOBALSETTINGS_HH

#include <algorithm> // for std::find if needed
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/**
 * @namespace CHeT
 * @brief Root namespace for the Cylindrical Helix Tracker (CHeT) project.
 */
namespace CHeT
{

/**
 * @namespace CHeT::Config
 * @brief Global Settings namespace containing configuration structures,
 * physical constants, and utility functions for the detector geometry.
 */
namespace Config
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
    double z_loc; ///< Longitudinal position [mm] (Local)
    double x_loc; ///< Transverse position X [mm] (Local)
    double y_loc; ///< Transverse position Y [mm] (Local)
    int cylinderId; ///< ID of the cylinder where intersection occurs
    double widthZ; ///< Longitudinal uncertainty due to fiber thickness [mm]
};

// --- Global Physical Constants ---

// Kept as constexpr in header for compile-time optimization
constexpr double L_HALF = 150.0; ///< Half-length of the detector [mm]
constexpr double BUNDLE_WIDTH = 2.0; ///< Physical width of a bundle [mm]

// --- Experimental Offsets ---

/**
 * @brief Retrieves the current experimental angular offset (OFFSET_EXP).
 * @return The offset value in radians.
 */
double GetOffsetExp();

/**
 * @brief Sets the experimental angular offset (OFFSET_EXP).
 * @note Setting this value invalidates the geometry cache, causing a lazy
 * recalculation.
 * @param val The new offset value in radians.
 */
void SetOffsetExp(double val);

/**
 * @brief Retrieves the current experimental parameter DELTA1.
 * Used in the configuration of Cylinder 1.
 * @return The value of DELTA1.
 */
double GetDelta1();

/**
 * @brief Sets the experimental parameter DELTA1.
 * @note Setting this value invalidates the geometry cache, causing a lazy
 * recalculation.
 * @param val The new value for DELTA1.
 */
void SetDelta1(double val);

/**
 * @brief Retrieves the current experimental parameter DELTA2.
 * Used in the configuration of Cylinder 2.
 * @return The value of DELTA2.
 */
double GetDelta2();

/**
 * @brief Sets the experimental parameter DELTA2.
 * @note Setting this value invalidates the geometry cache, causing a lazy
 * recalculation.
 * @param val The new value for DELTA2.
 */
void SetDelta2(double val);

// --- Rotation Settings ---

/**
 * @brief Sets the global rotation of the detector using Euler angles.
 * Rotations are applied in order: X, then Y, then Z.
 * @param rx Rotation around X axis [rad]
 * @param ry Rotation around Y axis [rad]
 * @param rz Rotation around Z axis [rad]
 */
void SetRotation(double rx, double ry, double rz);

/**
 * @brief Retrieves the current global rotation angles.
 * @param rx Output rotation around X axis [rad]
 * @param ry Output rotation around Y axis [rad]
 * @param rz Output rotation around Z axis [rad]
 */
void GetRotation(double &rx, double &ry, double &rz);

/**
 * @brief Sets the global translation of the detector origin.
 * @param tx Translation X [mm]
 * @param ty Translation Y [mm]
 * @param tz Translation Z [mm]
 */
void SetTranslation(double tx, double ty, double tz);

/**
 * @brief Retrieves the current global translation.
 * @param tx Output Translation X [mm]
 * @param ty Output Translation Y [mm]
 * @param tz Output Translation Z [mm]
 */
void GetTranslation(double &tx, double &ty, double &tz);

/**
 * @brief Rotates a 3D VECTOR in-place from the Detector Local Frame to the
 * Global Lab Frame. Uses the currently set rotation angles.
 * This applies ONLY rotation (suitable for direction vectors).
 * @param x X component
 * @param y Y component
 * @param z Z component
 */
void ApplyRotation(double &x, double &y, double &z);

/**
 * @brief Rotates a 3D VECTOR in-place from the Global Lab Frame to the
 * Detector Local Frame. (Inverse of ApplyRotation).
 * This applies ONLY rotation (suitable for direction vectors).
 * @param x X component
 * @param y Y component
 * @param z Z component
 */
void ApplyInverseRotation(double &x, double &y, double &z);

/**
 * @brief Transforms a 3D POINT in-place from the Detector Local Frame to the
 * Global Lab Frame. Applies Rotation then Translation.
 * P_global = R * P_local + T
 * @param x X coordinate [mm]
 * @param y Y coordinate [mm]
 * @param z Z coordinate [mm]
 */
void ApplyTransformation(double &x, double &y, double &z);

/**
 * @brief Transforms a 3D POINT in-place from the Global Lab Frame to the
 * Detector Local Frame. Applies Inverse Translation then Inverse Rotation.
 * P_local = R^T * (P_global - T)
 * @param x X coordinate [mm]
 * @param y Y coordinate [mm]
 * @param z Z coordinate [mm]
 */
void ApplyInverseTransformation(double &x, double &y, double &z);

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
std::vector<BundlesIntersection> FindIntersections(const std::vector<int> &hit_ids);

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

} // namespace Config

} // namespace CHeT

#endif // CHETGLOBALSETTINGS_HH
