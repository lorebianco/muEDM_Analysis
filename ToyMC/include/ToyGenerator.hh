#ifndef TOYGENERATOR_HH
#define TOYGENERATOR_HH

#include <vector>

namespace ToyMC
{

/**
 * @brief Structure representing a generated Cosmic Ray track (straight line).
 */
struct CosmicTrack
{
    double x0; ///< X origin
    double y0; ///< Y origin
    double z0; ///< Z origin
    double ux; ///< X direction cosine
    double uy; ///< Y direction cosine
    double uz; ///< Z direction cosine
    bool active; ///< Flag if track is active/accepted
};

/**
 * @brief Structure representing a generated Michel Electron track (helix).
 */
struct MichelTrack
{
    double E_kin; ///< Kinetic energy
    double x0, y0; ///< Emission point on the orbit plane
    double cx, cy, radius; ///< 2D Helix parameters
    double z0, dz_dt; ///< Longitudinal Helix parameters
    double t_min, t_max; ///< Angular limits of traversal
    double theta_rad, phi_dir; ///< Kinematic angles at vertex
    bool hits_boundary; ///< Did the track hit radial or Z boundaries?
};

// --- Generation Functions ---

/**
 * @brief Generates a random cosmic muon track traversing the detector.
 * @return A CosmicTrack object with truth parameters.
 */
CosmicTrack GenerateCosmic();

/**
 * @brief Finds the global bundle IDs intersected by a cosmic track.
 * @param track The true cosmic track.
 * @param efficiency The single-layer detection efficiency [0.0, 1.0].
 * @return Vector of global bundle IDs that recorded a hit.
 */
std::vector<int> FindCosmicHits(const CosmicTrack &track, double efficiency = 1.0);

/**
 * @brief Generates a random energy for a Michel electron according to the Michel spectrum.
 * @return Kinetic energy [MeV].
 */
double GenerateMichelEnergy();

/**
 * @brief Generates a random Michel electron track from the central target.
 * @param requireHitsBoundary If true, forces generation until the track leaves the tracking volume.
 * @return A MichelTrack object with truth parameters.
 */
MichelTrack GenerateMichelTrack(bool requireHitsBoundary = false);

/**
 * @brief Finds the global bundle IDs intersected by a Michel electron helix.
 * @param tr The true Michel track.
 * @param efficiency The single-layer detection efficiency [0.0, 1.0].
 * @return Vector of global bundle IDs that recorded a hit.
 */
std::vector<int> FindMichelHits(const MichelTrack &tr, double efficiency = 1.0);

} // namespace ToyMC

#endif // TOYGENERATOR_HH
