#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

#include <TMath.h>
#include <TRandom3.h>

#include "CHeT/CHeTGlobalSettings.hh"
#include "ToyGenerator.hh"

using namespace std;
using namespace CHeT;

namespace ToyMC
{

CosmicTrack GenerateCosmic()
{
    static TRandom3 rnd(0);
    CosmicTrack track;
    Double_t halfX_up = 90.;
    Double_t halfZ_up = 90.;
    Double_t yUp = 386.0;
    Double_t halfX_down = 90.;
    Double_t halfZ_down = 90.;
    Double_t yDown = -476.0;
    Bool_t accepted = false;

    while(!accepted)
    {
        Double_t cosTheta = rnd.Uniform(-1.0, 0.0);
        Double_t cosThetaSq_test = rnd.Uniform(0.0, 1.0);
        if(cosThetaSq_test > cosTheta * cosTheta)
            continue;

        Double_t phi = rnd.Uniform(0.0, 2.0 * M_PI);
        Double_t sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        track.ux = sinTheta * cos(phi);
        track.uy = cosTheta;
        track.uz = sinTheta * sin(phi);

        track.x0 = rnd.Uniform(-halfX_up, halfX_up);
        track.y0 = yUp;
        track.z0 = rnd.Uniform(-halfZ_up, halfZ_up);

        Double_t t = (yDown - track.y0) / track.uy;
        Double_t x_proj = track.x0 + track.ux * t;
        Double_t z_proj = track.z0 + track.uz * t;

        if(abs(x_proj) <= halfX_down && abs(z_proj) <= halfZ_down)
        {
            accepted = true;
        }
    }
    track.active = true;
    return track;
}

HitResult FindCosmicHits(const CosmicTrack &track, double efficiency)
{
    std::map<int, std::tuple<double, double, double>> hit_map;
    auto cylinders = Config::GetCylinders();
    int current_global_offset = 0;

    for(const auto &cyl : cylinders)
    {
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int i_lay = 0; i_lay < 2; ++i_lay)
        {
            const auto &lay = *layers[i_lay];
            Double_t R = lay.radius;
            Double_t A = track.ux * track.ux + track.uy * track.uy;
            Double_t B = 2.0 * (track.x0 * track.ux + track.y0 * track.uy);
            Double_t C = track.x0 * track.x0 + track.y0 * track.y0 - R * R;
            Double_t delta = B * B - 4 * A * C;

            if(delta >= 0)
            {
                Double_t t_sol[2]
                    = { (-B + sqrt(delta)) / (2.0 * A), (-B - sqrt(delta)) / (2.0 * A) };
                for(int i = 0; i < 2; ++i)
                {
                    Double_t t = t_sol[i];
                    Double_t zi = track.z0 + track.uz * t;
                    if(abs(zi) <= Config::L_HALF)
                    {
                        Double_t xi = track.x0 + track.ux * t;
                        Double_t yi = track.y0 + track.uy * t;
                        Double_t phi_track = atan2(yi, xi);

                        for(int b = 0; b < lay.nBundles; ++b)
                        {
                            int b_id = current_global_offset + b;
                            Config::FiberProp p = Config::GetFiberProp(b_id);
                            Double_t alpha = (zi + Config::L_HALF) / (2.0 * Config::L_HALF);
                            Double_t phi_f = p.phi0 + p.dir * alpha * M_PI;
                            Double_t dphi
                                = abs(Config::wrap0_2pi(phi_track) - Config::wrap0_2pi(phi_f));
                            if(dphi > M_PI)
                                dphi = 2.0 * M_PI - dphi;

                            if(dphi < (M_PI / lay.nBundles) && gRandom->Rndm() <= efficiency)
                            {
                                hit_map[b_id] = { xi, yi, zi };
                            }
                        }
                    }
                }
            }
            current_global_offset += lay.nBundles;
        }
    }
    HitResult res;
    for(const auto &kv : hit_map)
    {
        res.bundles.push_back(kv.first);
        res.x.push_back(std::get<0>(kv.second));
        res.y.push_back(std::get<1>(kv.second));
        res.z.push_back(std::get<2>(kv.second));
    }
    return res;
}

double GenerateMichelEnergy()
{
    const double E_max = 52.8;

    while(true)
    {
        double x = gRandom->Uniform(0.0, E_max);
        double y = gRandom->Uniform(0.0, 1.0);
        double w = x / E_max;
        double f = 2.0 * w * w * (3.0 - 2.0 * w); // Standard Michel spectrum for unpolarized muons
        if(y <= f)
            return x;
    }
}

MichelTrack GenerateMichelTrack(bool requireHitsBoundary)
{
    const double B_tesla = 2.89;
    const double m_e = 0.511;
    const double R_start = 30.0;
    const double R_limit = 199.0 / 2.0;

    MichelTrack tr;
    bool accepted = false;

    while(!accepted)
    {
        tr.E_kin = GenerateMichelEnergy();
        double gamma = gRandom->Uniform(0.0, 2.0 * M_PI);
        double cosTheta = gRandom->Uniform(-0.99, 0.99);
        tr.theta_rad = std::acos(cosTheta);
        tr.phi_dir = gRandom->Uniform(0.0, 2.0 * M_PI);

        tr.x0 = R_start * std::cos(gamma);
        tr.y0 = R_start * std::sin(gamma);
        double p = std::sqrt(tr.E_kin * tr.E_kin + 2.0 * m_e * tr.E_kin);
        double pt = p * std::sin(tr.theta_rad);
        tr.radius = pt / (0.29979 * std::abs(B_tesla));

        double omega0 = (B_tesla > 0) ? -1.0 : 1.0;
        tr.t_min = tr.phi_dir - omega0 * (M_PI / 2.0);
        while(tr.t_min > M_PI)
            tr.t_min -= 2.0 * M_PI;
        while(tr.t_min < -M_PI)
            tr.t_min += 2.0 * M_PI;

        tr.cx = tr.x0 - tr.radius * std::cos(tr.t_min);
        tr.cy = tr.y0 - tr.radius * std::sin(tr.t_min);
        double rho_c = std::sqrt(tr.cx * tr.cx + tr.cy * tr.cy);
        double phi_c = std::atan2(tr.cy, tr.cx);

        tr.dz_dt = omega0 * tr.radius / std::tan(tr.theta_rad);
        tr.z0 = -tr.dz_dt * tr.t_min;
        double delta_t_z
            = (std::abs(tr.dz_dt) > 1e-9) ? (Config::L_HALF / std::abs(tr.dz_dt)) : 1e9;
        double t_end_z = tr.t_min + omega0 * delta_t_z;
        tr.t_max = t_end_z;
        tr.hits_boundary = false;

        double arg = (R_limit * R_limit - rho_c * rho_c - tr.radius * tr.radius)
            / (2.0 * tr.radius * rho_c);
        if(std::abs(arg) <= 1.0)
        {
            double dt_acos = std::acos(arg);
            double closest_dt = 1e9;
            bool found_radial = false;
            for(double ang : { phi_c + dt_acos, phi_c - dt_acos })
            {
                double diff = ang - tr.t_min;
                if(omega0 < 0)
                {
                    while(diff > 0)
                        diff -= 2.0 * M_PI;
                    while(diff < -2.0 * M_PI)
                        diff += 2.0 * M_PI;
                }
                else
                {
                    while(diff < 0)
                        diff += 2.0 * M_PI;
                    while(diff > 2.0 * M_PI)
                        diff += 2.0 * M_PI;
                }
                if(std::abs(diff) < std::abs(closest_dt))
                {
                    closest_dt = diff;
                    found_radial = true;
                }
            }
            if(found_radial)
            {
                double t_radial = tr.t_min + closest_dt;
                if(std::abs(t_radial - tr.t_min) < std::abs(t_end_z - tr.t_min))
                {
                    tr.t_max = t_radial;
                    tr.hits_boundary = true;
                }
            }
        }

        if(!requireHitsBoundary || tr.hits_boundary)
        {
            accepted = true;
        }
    }
    return tr;
}

HitResult FindMichelHits(const MichelTrack &tr, double efficiency)
{
    std::map<int, std::tuple<double, double, double>> hit_map;
    auto cylinders = Config::GetCylinders();
    int current_global_offset = 0;

    double rho_c = std::sqrt(tr.cx * tr.cx + tr.cy * tr.cy);
    double phi_c = std::atan2(tr.cy, tr.cx);
    double t_start_search = std::min(tr.t_min, tr.t_max);
    double t_end_search = std::max(tr.t_min, tr.t_max);

    for(const auto &cyl : cylinders)
    {
        const Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
        for(int l = 0; l < 2; ++l)
        {
            const auto &lay = *layers[l];
            double arg_layer = (lay.radius * lay.radius - rho_c * rho_c - tr.radius * tr.radius)
                / (2.0 * tr.radius * rho_c);

            if(std::abs(arg_layer) <= 1.0)
            {
                double dt_layer = std::acos(arg_layer);
                double t_base[2] = { phi_c + dt_layer, phi_c - dt_layer };

                for(int k = -25; k <= 25; ++k)
                {
                    for(double t_b : t_base)
                    {
                        double t = t_b + 2.0 * M_PI * k;
                        if(t >= t_start_search && t <= t_end_search)
                        {
                            double z = tr.z0 + tr.dz_dt * t;
                            if(std::abs(z) <= Config::L_HALF)
                            {
                                double x = tr.cx + tr.radius * std::cos(t);
                                double y = tr.cy + tr.radius * std::sin(t);
                                double phi_hit = std::atan2(y, x);

                                for(int b = 0; b < lay.nBundles; ++b)
                                {
                                    int b_id = current_global_offset + b;
                                    auto prop = Config::GetFiberProp(b_id);
                                    double alpha = (z + Config::L_HALF) / (2.0 * Config::L_HALF);
                                    double phi_f = prop.phi0 + prop.dir * alpha * M_PI;

                                    double dphi = std::abs(
                                        Config::wrap0_2pi(phi_hit) - Config::wrap0_2pi(phi_f));
                                    if(dphi > M_PI)
                                        dphi = 2.0 * M_PI - dphi;

                                    if(dphi < (M_PI / lay.nBundles)
                                        && gRandom->Rndm() <= efficiency)
                                    {
                                        hit_map[b_id] = { x, y, z };
                                    }
                                }
                            }
                        }
                    }
                }
            }
            current_global_offset += lay.nBundles;
        }
    }
    HitResult res;
    for(const auto &kv : hit_map)
    {
        res.bundles.push_back(kv.first);
        res.x.push_back(std::get<0>(kv.second));
        res.y.push_back(std::get<1>(kv.second));
        res.z.push_back(std::get<2>(kv.second));
    }
    return res;
}

} // namespace ToyMC
