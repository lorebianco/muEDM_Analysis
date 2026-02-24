#include "CHeT/CHeTGlobalSettings.hh"

// Standard TMath or cmath defines check
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace CHeT
{
namespace Config
{

// --- Internal Data Structures (Hidden from Header) ---

// Backing variables for geometry settings
static double g_OFFSET_EXP = 40.0 * (M_PI / 180.0);
static double g_DELTA1 = 1.267;
static double g_DELTA2 = 1.420 + (18.0 * M_PI / 180.0);
static bool g_GeometryDirty = true;
static std::vector<CylinderConfig> g_Cylinders;

// Rotation Globals (Euler angles in radians)
static double g_RotX = 0.0;
static double g_RotY = 0.0;
static double g_RotZ = 0.0;
// Precomputed Rotation Matrix (Local -> Global)
static double g_R[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
static bool g_IsRotated = false;

// Translation Globals [mm]
static double g_TransX = 0.0;
static double g_TransY = 0.0;
static double g_TransZ = 0.0;

// Static map for board offsets
static const std::map<int, int> BOARD_OFFSETS = {
    { 0, 0 }, // Board 0: Cylinder 0 (Inner)
    { 1, 0 }, // Board 1: Cylinder 0 (Outer)
    { 2, 94 }, // Board 2: After Cylinder 0 (45+49)
    { 3, 94 } // Board 3: After Cylinder 0
    //{4, 213}, // Board 4: Hypothetical Cyl 2 (94+119)
};

// Huge map for channel mapping. Defined here to keep header clean.
static const std::map<int, std::map<int, int>> BOARD_MAPS = {
    { 0,
        { { 31, 0 }, { 27, 1 }, { 23, 2 }, { 19, 3 }, { 15, 4 }, { 11, 5 }, { 7, 6 }, { 3, 7 },
            { 29, 8 }, { 25, 9 }, { 21, 10 }, { 17, 11 }, { 13, 12 }, { 9, 13 }, { 5, 14 },
            { 1, 15 }, { 0, 16 }, { 4, 17 }, { 8, 18 }, { 12, 19 }, { 16, 20 }, { 20, 21 },
            { 24, 22 }, { 28, 23 }, { 2, 24 }, { 6, 25 }, { 10, 26 }, { 14, 27 }, { 18, 28 },
            { 22, 29 }, { 26, 30 }, { 30, 31 }, { 32, 32 }, { 36, 33 }, { 40, 34 }, { 44, 35 },
            { 48, 36 }, { 52, 37 }, { 56, 38 }, { 60, 39 }, { 34, 40 }, { 38, 41 }, { 42, 42 },
            { 46, 43 }, { 50, 44 }, { 54, 45 }, { 58, 46 }, { 62, 47 } } },
    { 1,
        { { 31, 48 }, { 27, 49 }, { 23, 50 }, { 19, 51 }, { 15, 52 }, { 11, 53 }, { 7, 54 },
            { 3, 55 }, { 29, 56 }, { 25, 57 }, { 21, 58 }, { 17, 59 }, { 13, 60 }, { 9, 61 },
            { 5, 62 }, { 1, 63 }, { 0, 64 }, { 4, 65 }, { 8, 66 }, { 12, 67 }, { 16, 68 },
            { 20, 69 }, { 24, 70 }, { 28, 71 }, { 2, 72 }, { 6, 73 }, { 10, 74 }, { 14, 75 },
            { 18, 76 }, { 22, 77 }, { 26, 78 }, { 30, 79 }, { 32, 80 }, { 36, 81 }, { 40, 82 },
            { 44, 83 }, { 48, 84 }, { 52, 85 }, { 56, 86 }, { 60, 87 }, { 34, 88 }, { 38, 89 },
            { 42, 90 }, { 46, 91 }, { 50, 92 }, { 54, 93 } } },
    { 2,
        { { 31, 0 }, { 27, 1 }, { 23, 2 }, { 19, 3 }, { 15, 4 }, { 11, 5 }, { 7, 6 }, { 3, 7 },
            { 29, 8 }, { 25, 9 }, { 21, 10 }, { 17, 11 }, { 13, 12 }, { 9, 13 }, { 5, 14 },
            { 1, 15 }, { 0, 16 }, { 4, 17 }, { 8, 18 }, { 12, 19 }, { 16, 20 }, { 20, 21 },
            { 24, 22 }, { 28, 23 }, { 2, 24 }, { 6, 25 }, { 10, 26 }, { 14, 27 }, { 18, 28 },
            { 22, 29 }, { 26, 30 }, { 30, 31 }, { 63, 32 }, { 59, 33 }, { 55, 34 }, { 51, 35 },
            { 47, 36 }, { 43, 37 }, { 39, 38 }, { 35, 39 }, { 61, 40 }, { 57, 41 }, { 53, 42 },
            { 49, 43 }, { 45, 44 }, { 41, 45 }, { 37, 46 }, { 33, 47 }, { 32, 48 }, { 36, 49 },
            { 40, 50 }, { 44, 51 }, { 48, 52 }, { 52, 53 }, { 56, 54 }, { 60, 55 }, { 34, 56 },
            { 38, 57 }, { 42, 58 }, { 46, 59 }, { 50, 60 }, { 54, 61 }, { 58, 62 }, { 62, 63 } } },
    { 3,
        { { 31, 64 }, { 27, 65 }, { 23, 66 }, { 19, 67 }, { 15, 68 }, { 11, 69 }, { 7, 70 },
            { 3, 71 }, { 29, 72 }, { 25, 73 }, { 21, 74 }, { 17, 75 }, { 13, 76 }, { 9, 77 },
            { 5, 78 }, { 1, 79 }, { 0, 80 }, { 4, 81 }, { 8, 82 }, { 12, 83 }, { 16, 84 },
            { 20, 85 }, { 24, 86 }, { 28, 87 }, { 2, 88 }, { 6, 89 }, { 10, 90 }, { 14, 91 },
            { 18, 92 }, { 22, 93 }, { 26, 94 }, { 30, 95 }, { 63, 96 }, { 59, 97 }, { 55, 98 },
            { 51, 99 }, { 47, 100 }, { 43, 101 }, { 39, 102 }, { 35, 103 }, { 61, 104 },
            { 57, 105 }, { 53, 106 }, { 49, 107 }, { 45, 108 }, { 41, 109 }, { 37, 110 },
            { 33, 111 }, { 32, 112 }, { 36, 113 }, { 40, 114 }, { 44, 115 }, { 48, 116 },
            { 52, 117 }, { 56, 118 } } }
};

// --- Implementation ---

double GetOffsetExp()
{
    return g_OFFSET_EXP;
}
void SetOffsetExp(double val)
{
    g_OFFSET_EXP = val;
    g_GeometryDirty = true;
}

double GetDelta1()
{
    return g_DELTA1;
}
void SetDelta1(double val)
{
    g_DELTA1 = val;
    g_GeometryDirty = true;
}

double GetDelta2()
{
    return g_DELTA2;
}
void SetDelta2(double val)
{
    g_DELTA2 = val;
    g_GeometryDirty = true;
}

void SetRotation(double rx, double ry, double rz)
{
    g_RotX = rx;
    g_RotY = ry;
    g_RotZ = rz;
    g_IsRotated = (std::abs(rx) > 1e-9 || std::abs(ry) > 1e-9 || std::abs(rz) > 1e-9);

    if(g_IsRotated)
    {
        double cx = std::cos(rx), sx = std::sin(rx);
        double cy = std::cos(ry), sy = std::sin(ry);
        double cz = std::cos(rz), sz = std::sin(rz);

        // Calculate Rotation Matrix Columns by transforming basis vectors
        auto rotate = [&](double x, double y, double z, double &ox, double &oy, double &oz)
        {
            // Rx
            double y1 = y * cx - z * sx;
            double z1 = y * sx + z * cx;
            double x1 = x;
            // Ry
            double z2 = z1 * cy - x1 * sy;
            double x2 = z1 * sy + x1 * cy;
            double y2 = y1;
            // Rz
            double x3 = x2 * cz - y2 * sz;
            double y3 = x2 * sz + y2 * cz;
            double z3 = z2;
            ox = x3;
            oy = y3;
            oz = z3;
        };

        // Col 0
        rotate(1, 0, 0, g_R[0][0], g_R[1][0], g_R[2][0]);
        // Col 1
        rotate(0, 1, 0, g_R[0][1], g_R[1][1], g_R[2][1]);
        // Col 2
        rotate(0, 0, 1, g_R[0][2], g_R[1][2], g_R[2][2]);
    }
    else
    {
        // Identity
        g_R[0][0] = 1;
        g_R[0][1] = 0;
        g_R[0][2] = 0;
        g_R[1][0] = 0;
        g_R[1][1] = 1;
        g_R[1][2] = 0;
        g_R[2][0] = 0;
        g_R[2][1] = 0;
        g_R[2][2] = 1;
    }
}

void GetRotation(double &rx, double &ry, double &rz)
{
    rx = g_RotX;
    ry = g_RotY;
    rz = g_RotZ;
}

void SetTranslation(double tx, double ty, double tz)
{
    g_TransX = tx;
    g_TransY = ty;
    g_TransZ = tz;
}

void GetTranslation(double &tx, double &ty, double &tz)
{
    tx = g_TransX;
    ty = g_TransY;
    tz = g_TransZ;
}

void ApplyRotation(double &x, double &y, double &z)
{
    if(!g_IsRotated)
        return;

    double xn = g_R[0][0] * x + g_R[0][1] * y + g_R[0][2] * z;
    double yn = g_R[1][0] * x + g_R[1][1] * y + g_R[1][2] * z;
    double zn = g_R[2][0] * x + g_R[2][1] * y + g_R[2][2] * z;

    x = xn;
    y = yn;
    z = zn;
}

void ApplyInverseRotation(double &x, double &y, double &z)
{
    if(!g_IsRotated)
        return;

    // Transpose multiplication
    double xn = g_R[0][0] * x + g_R[1][0] * y + g_R[2][0] * z;
    double yn = g_R[0][1] * x + g_R[1][1] * y + g_R[2][1] * z;
    double zn = g_R[0][2] * x + g_R[1][2] * y + g_R[2][2] * z;

    x = xn;
    y = yn;
    z = zn;
}

void ApplyTransformation(double &x, double &y, double &z)
{
    // 1. Rotation (Local -> Aligned Global)
    ApplyRotation(x, y, z);

    // 2. Translation (Aligned Global -> Global)
    x += g_TransX;
    y += g_TransY;
    z += g_TransZ;
}

void ApplyInverseTransformation(double &x, double &y, double &z)
{
    // 1. Inverse Translation (Global -> Aligned Global)
    x -= g_TransX;
    y -= g_TransY;
    z -= g_TransZ;

    // 2. Inverse Rotation (Aligned Global -> Local)
    ApplyInverseRotation(x, y, z);
}

int GetBoardGlobalOffset(int board_id)
{
    auto it = BOARD_OFFSETS.find(board_id);
    return (it != BOARD_OFFSETS.end()) ? it->second : -1;
}

const std::vector<CylinderConfig> &GetCylinders()
{
    if(g_GeometryDirty || g_Cylinders.empty())
    {
        g_Cylinders = { {
                            0, // Cylinder 1
                            { 45, 17.0, 4.294 + g_DELTA1 + g_OFFSET_EXP, -1, 632 }, // Inner (Red)
                            { 49, 17.0, 3.829 + g_DELTA1 + g_OFFSET_EXP, 1, 807 } // Outer (Orange)
                        },
            {
                1, // Cylinder 2
                { 59, 21.0, 2.609 + g_DELTA2 + g_OFFSET_EXP, -1, 600 }, // Inner (Blue)
                { 60, 21.0, 3.351 + g_DELTA2 + g_OFFSET_EXP, 1, 432 } // Outer (Cyan)
            } };
        g_GeometryDirty = false;
    }
    return g_Cylinders;
}

double wrap0_2pi(double angle)
{
    double wrapped = std::fmod(angle, 2.0 * M_PI);
    if(wrapped < 0)
        wrapped += 2.0 * M_PI;
    return wrapped;
}

FiberProp GetFiberProp(int b_id)
{
    const auto &cylinders = GetCylinders();
    int offset = 0;

    for(const auto &cyl : cylinders)
    {
        // Check Inner
        if(b_id < offset + cyl.inner.nBundles)
        {
            int b_loc = b_id - offset;
            return { cyl.id, 0, cyl.inner.radius,
                wrap0_2pi(2.0 * M_PI * (-b_loc) / cyl.inner.nBundles + cyl.inner.phiOffset),
                cyl.inner.direction, cyl.inner.color };
        }
        offset += cyl.inner.nBundles;

        // Check Outer
        if(b_id < offset + cyl.outer.nBundles)
        {
            int b_loc = b_id - offset;
            return { cyl.id, 1, cyl.outer.radius,
                wrap0_2pi(2.0 * M_PI * (b_loc) / cyl.outer.nBundles + cyl.outer.phiOffset),
                cyl.outer.direction, cyl.outer.color };
        }
        offset += cyl.outer.nBundles;
    }
    // Return invalid prop
    return { -1, -1, 0, 0, 0, 0 };
}

int GetGlobalBundleId(int board_id, int channel_id)
{
    auto boardMapIt = BOARD_MAPS.find(board_id);
    if(boardMapIt == BOARD_MAPS.end())
        return -1;

    auto chanIt = boardMapIt->second.find(channel_id);
    if(chanIt == boardMapIt->second.end())
        return -1;

    int local_val = chanIt->second;

    // Apply offset
    int boardOffset = GetBoardGlobalOffset(board_id);
    if(boardOffset == -1)
        return -1;

    return local_val + boardOffset;
}

std::vector<BundlesIntersection> FindIntersections(const std::vector<int> &hit_ids)
{
    std::vector<BundlesIntersection> out;
    // Optimization: reserve memory if hits are many, though tricky to guess how
    // many intersections

    for(size_t i = 0; i < hit_ids.size(); ++i)
    {
        for(size_t j = i + 1; j < hit_ids.size(); ++j)
        {
            FiberProp p1 = GetFiberProp(hit_ids[i]);
            FiberProp p2 = GetFiberProp(hit_ids[j]);

            // Intersection logic: Same cylinder, different layers
            if(p1.cylinderId == p2.cylinderId && p1.cylinderId != -1 && p1.layerId != p2.layerId)
            {
                double dphi0 = p2.phi0 - p1.phi0;
                double ddir = p1.dir - p2.dir;

                // Typical case: +1 vs -1 -> diff is 2 or -2.
                // Loop over possible periodic intersections
                for(int k = -2; k <= 2; ++k)
                {
                    // Avoid division by zero
                    if(std::abs(ddir) < 1e-9)
                        continue;

                    double alpha = (dphi0 + 2.0 * k * M_PI) / (ddir * M_PI);

                    if(alpha >= 0.0 && alpha <= 1.0)
                    {
                        double zi = alpha * (2.0 * L_HALF) - L_HALF;
                        double ph = p1.phi0 + p1.dir * alpha * M_PI;

                        // Z Thickness Calculation
                        // DeltaZ depends on Bundle Width and Helix slope.
                        double deltaZ = (BUNDLE_WIDTH / p1.r) * (L_HALF / M_PI);

                        // Calculate local coordinates
                        double x_loc = p1.r * std::cos(ph);
                        double y_loc = p1.r * std::sin(ph);
                        double z_loc = zi;

                        // Apply global rotation (to copies)
                        double x_rot = x_loc;
                        double y_rot = y_loc;
                        double z_rot = z_loc;
                        ApplyTransformation(x_rot, y_rot, z_rot);

                        out.push_back(
                            { z_rot, x_rot, y_rot, z_loc, x_loc, y_loc, p1.cylinderId, deltaZ });
                    }
                }
            }
        }
    }
    return out;
}

// --- Debug / Helper Functions Implementation ---

int GetGlobalIdFromGeometry(int cyl_id, int layer_id, int layer_idx)
{
    const auto &cylinders = GetCylinders();
    int global_offset = 0;

    for(const auto &cyl : cylinders)
    {
        if(cyl.id == cyl_id)
        {
            // Check limits
            int limit = (layer_id == 0) ? cyl.inner.nBundles : cyl.outer.nBundles;

            if(layer_idx < 0 || layer_idx >= limit)
            {
                std::cerr << " [ERROR] Index " << layer_idx << " out of range for Cylinder "
                          << cyl_id << " Layer " << layer_id << "!\n";
                return -1;
            }

            // Calculate ID: Offset + Local Base + Index
            int local_base = (layer_id == 0) ? 0 : cyl.inner.nBundles;
            return global_offset + local_base + layer_idx;
        }
        // Accumulate offset if not current cylinder
        global_offset += (cyl.inner.nBundles + cyl.outer.nBundles);
    }
    return -1; // Cylinder not found
}

void PrintBundleMapping(int global_id)
{
    if(global_id < 0)
        return;

    std::cout << "\n >>> INFO BUNDLE ID: " << global_id << " <<<\n";

    // A. Check Hardware Mapping
    bool found_hw = false;
    for(auto const &[board, map] : BOARD_MAPS)
    {
        int off = GetBoardGlobalOffset(board);
        if(off < 0)
            continue;

        // Reverse lookup: value == global_id - offset
        int target_local = global_id - off;

        for(auto const &[chan, val] : map)
        {
            if(val == target_local)
            {
                std::cout << " [HARDWARE] Board: " << board << " | Channel: " << chan << "\n";
                found_hw = true;
                break;
            }
        }
        if(found_hw)
            break;
    }
    if(!found_hw)
        std::cout << " [HARDWARE] Not mapped on any board.\n";

    // B. Check Geometry Mapping
    const auto &cylinders = GetCylinders();
    int current_offset = 0;
    bool found_geo = false;

    for(const auto &cyl : cylinders)
    {
        int tot = cyl.inner.nBundles + cyl.outer.nBundles;
        if(global_id >= current_offset && global_id < current_offset + tot)
        {
            int local = global_id - current_offset;
            bool isInner = (local < cyl.inner.nBundles);
            int layerIdx = isInner ? local : (local - cyl.inner.nBundles);

            std::cout << " [GEOMETRY] Cylinder: " << cyl.id
                      << " | Layer: " << (isInner ? "INNER" : "OUTER") << " | Index: " << layerIdx
                      << "\n";
            found_geo = true;
            break;
        }
        current_offset += tot;
    }

    if(!found_geo)
        std::cout << " [GEOMETRY] ID out of geometric range.\n";
}

void MapExplorer()
{
    int mode;
    while(true)
    {
        std::cout << "\n-----------------------------------\n";
        std::cout << " 1. Info from Global ID\n";
        std::cout << " 2. Info from Board/Channel\n";
        std::cout << " 3. Info from Geometry (Cyl, Lay, Idx)\n";
        std::cout << " 4. First bundle of Cylinder (N/A in this version)\n"; // Placeholder
                                                                             // based
                                                                             // on
                                                                             // original
                                                                             // menu
                                                                             // logic
        std::cout << " 5. First bundle of Layer\n";
        std::cout << " 0. Quit\n";
        std::cout << "Choice: ";
        std::cin >> mode;

        if(!std::cin)
        { // Handle non-integer input to avoid infinite loop
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }

        if(mode == 0)
            break;

        int target_id = -1;

        if(mode == 1)
        {
            std::cout << "Global ID: ";
            std::cin >> target_id;
        }
        else if(mode == 2)
        {
            int b, c;
            std::cout << "Board: ";
            std::cin >> b;
            std::cout << "Chan: ";
            std::cin >> c;
            target_id = GetGlobalBundleId(b, c);
        }
        else if(mode >= 3 && mode <= 5)
        {
            int c, l = 0, i = 0;
            std::cout << "Cylinder ID: ";
            std::cin >> c;
            if(mode == 3 || mode == 5)
            {
                std::cout << "Layer (0=In, 1=Out): ";
                std::cin >> l;
            }
            if(mode == 3)
            {
                std::cout << "Index: ";
                std::cin >> i;
            }
            // Use the geometry helper
            target_id = GetGlobalIdFromGeometry(c, l, i);
        }

        // Print result
        if(target_id != -1)
            PrintBundleMapping(target_id);
        else
            std::cout << " -> Mapping not found or wrong parameters.\n";
    }
}

} // namespace Config
} // namespace CHeT
