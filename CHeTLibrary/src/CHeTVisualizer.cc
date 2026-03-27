#include "CHeT/CHeTGlobalSettings.hh"
#include "CHeT/CHeTVisualizer.hh"

// --- Internal Implementation Details ---
namespace
{

// Helper: Clipping algorithm to keep lines inside the detector box.
// Moved to anonymous namespace to keep the header clean.
bool ClipLineToBox(double x0, double y0, double z0, double ux, double uy, double uz, double xmin,
    double xmax, double ymin, double ymax, double zmin, double zmax, double &tmin, double &tmax)
{
    // Use bounds provided by caller (e.g. -inf, +inf for line, or 0, L for segment)

    auto clip = [&](double p0, double dp, double pmin, double pmax) -> bool
    {
        if(std::abs(dp) < 1e-12)
            return (p0 >= pmin && p0 <= pmax);

        double t1 = (pmin - p0) / dp;
        double t2 = (pmax - p0) / dp;
        if(t1 > t2)
            std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        return tmin <= tmax;
    };

    if(!clip(x0, ux, xmin, xmax))
        return false;
    if(!clip(y0, uy, ymin, ymax))
        return false;
    if(!clip(z0, uz, zmin, zmax))
        return false;

    return true;
}

} // anonymous namespace

namespace CHeT
{
namespace Vis
{

// --- Overloads for different track types in 2D ---

void RenderTrack2D(const VisLineTrack &tr, TMultiGraph *mg_xy)
{
    TGraph *gxy = new TGraph();
    // Use a large enough range so the line covers the detector view.
    // The view limits are set by Draw2DCore based on geometry.
    double range = 2000.0;
    double pts[2] = { -range, range };

    for(int i = 0; i < 2; i++)
    {
        if(std::abs(tr.ux) > 1e-6)
        {
            double t = (pts[i] - tr.x0) / tr.ux;
            gxy->SetPoint(i, pts[i], tr.y0 + tr.uy * t);
        }
        else
        {
            gxy->SetPoint(0, tr.x0, -range);
            gxy->SetPoint(1, tr.x0, range);
        }
    }
    gxy->SetLineColor(tr.color);
    gxy->SetLineWidth(tr.width);
    gxy->SetLineStyle(tr.style);
    mg_xy->Add(gxy, "L");
}

void RenderTrack2D(const VisHelixTrack &tr, TMultiGraph *mg_xy)
{
    int n_pts = 100;
    TGraph *gxy = new TGraph(n_pts);
    double dt = (tr.t_max - tr.t_min) / (n_pts - 1);
    for(int i = 0; i < n_pts; ++i)
    {
        double t = tr.t_min + i * dt;
        gxy->SetPoint(i, tr.cx + tr.radius * std::cos(t), tr.cy + tr.radius * std::sin(t));
    }
    gxy->SetLineColor(tr.color);
    gxy->SetLineWidth(tr.width);
    gxy->SetLineStyle(tr.style);
    mg_xy->Add(gxy, "L");
}

void RenderTrack2D(const VisGenericTrack &tr, TMultiGraph *mg_xy)
{
    if(tr.points.empty())
        return;
    TGraph *gxy = new TGraph(tr.points.size());
    for(size_t i = 0; i < tr.points.size(); ++i)
    {
        gxy->SetPoint(i, tr.points[i].x, tr.points[i].y);
    }
    gxy->SetLineColor(tr.color);
    gxy->SetLineWidth(tr.width);
    gxy->SetLineStyle(tr.style);
    mg_xy->Add(gxy, "L");
}

// --- ZX Projection (X vs Z) ---
void RenderTrackZX(const VisLineTrack &tr, TMultiGraph *mg_zx)
{
    TGraph *g = new TGraph();
    double range = 2000.0;
    double pts[2] = { -range, range };

    for(int i = 0; i < 2; i++)
    {
        if(std::abs(tr.ux) > 1e-6)
        {
            double t = (pts[i] - tr.x0) / tr.ux;
            g->SetPoint(i, pts[i], tr.z0 + tr.uz * t);
        }
        else
        {
            g->SetPoint(0, tr.x0, -range);
            g->SetPoint(1, tr.x0, range);
        }
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zx->Add(g, "L");
}

void RenderTrackZX(const VisHelixTrack &tr, TMultiGraph *mg_zx)
{
    int n_pts = 100;
    TGraph *g = new TGraph(n_pts);
    double dt = (tr.t_max - tr.t_min) / (n_pts - 1);
    for(int i = 0; i < n_pts; ++i)
    {
        double t = tr.t_min + i * dt;
        double x = tr.cx + tr.radius * std::cos(t);
        double z = tr.z0 + tr.dz_dt * t;
        g->SetPoint(i, x, z);
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zx->Add(g, "L");
}

void RenderTrackZX(const VisGenericTrack &tr, TMultiGraph *mg_zx)
{
    if(tr.points.empty())
        return;
    TGraph *g = new TGraph(tr.points.size());
    for(size_t i = 0; i < tr.points.size(); ++i)
        g->SetPoint(i, tr.points[i].x, tr.points[i].z);
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zx->Add(g, "L");
}

// --- ZY Projection (Z vs Y) ---
void RenderTrackZY(const VisLineTrack &tr, TMultiGraph *mg_zy)
{
    TGraph *g = new TGraph();
    double range = 2000.0;
    double pts[2] = { -range, range }; // Z is horizontal

    for(int i = 0; i < 2; i++)
    {
        if(std::abs(tr.uz) > 1e-6)
        {
            double t = (pts[i] - tr.z0) / tr.uz;
            g->SetPoint(i, pts[i], tr.y0 + tr.uy * t);
        }
        else
        {
            g->SetPoint(0, tr.z0, -range);
            g->SetPoint(1, tr.z0, range);
        }
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zy->Add(g, "L");
}

void RenderTrackZY(const VisHelixTrack &tr, TMultiGraph *mg_zy)
{
    int n_pts = 100;
    TGraph *g = new TGraph(n_pts);
    double dt = (tr.t_max - tr.t_min) / (n_pts - 1);
    for(int i = 0; i < n_pts; ++i)
    {
        double t = tr.t_min + i * dt;
        double y = tr.cy + tr.radius * std::sin(t);
        double z = tr.z0 + tr.dz_dt * t;
        g->SetPoint(i, z, y);
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zy->Add(g, "L");
}

void RenderTrackZY(const VisGenericTrack &tr, TMultiGraph *mg_zy)
{
    if(tr.points.empty())
        return;
    TGraph *g = new TGraph(tr.points.size());
    for(size_t i = 0; i < tr.points.size(); ++i)
        g->SetPoint(i, tr.points[i].z, tr.points[i].y);
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_zy->Add(g, "L");
}

// --- PhiZ Projection (Phi vs Z) ---
void RenderTrackPhiZ(const VisLineTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi)
{
    int n_pts = 200;
    TGraph *g = new TGraph();
    double range = 1000.0;
    int cnt = 0;
    for(int i = 0; i < n_pts; ++i)
    {
        double t = -range + i * (2 * range) / (n_pts - 1);
        double x = tr.x0 + tr.ux * t;
        double y = tr.y0 + tr.uy * t;
        double z = tr.z0 + tr.uz * t;

        if(std::abs(z) > 400)
            continue;

        double phi = std::atan2(y, x);
        if(wrap_phi && phi < 0)
            phi += 2 * M_PI;
        g->SetPoint(cnt++, phi, z);
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_phiz->Add(g, "L");
}

void RenderTrackPhiZ(const VisHelixTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi)
{
    int n_pts = 100;
    TGraph *g = new TGraph(n_pts);
    double dt = (tr.t_max - tr.t_min) / (n_pts - 1);
    for(int i = 0; i < n_pts; ++i)
    {
        double t = tr.t_min + i * dt;
        double x = tr.cx + tr.radius * std::cos(t);
        double y = tr.cy + tr.radius * std::sin(t);
        double z = tr.z0 + tr.dz_dt * t;

        double phi = std::atan2(y, x);
        if(wrap_phi && phi < 0)
            phi += 2 * M_PI;

        g->SetPoint(i, phi, z);
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_phiz->Add(g, "L");
}

void RenderTrackPhiZ(const VisGenericTrack &tr, TMultiGraph *mg_phiz, bool wrap_phi)
{
    if(tr.points.empty())
        return;
    TGraph *g = new TGraph(tr.points.size());
    for(size_t i = 0; i < tr.points.size(); ++i)
    {
        double phi = std::atan2(tr.points[i].y, tr.points[i].x);
        if(wrap_phi && phi < 0)
            phi += 2 * M_PI;
        g->SetPoint(i, phi, tr.points[i].z);
    }
    g->SetLineColor(tr.color);
    g->SetLineWidth(tr.width);
    g->SetLineStyle(tr.style);
    mg_phiz->Add(g, "L");
}

// --- Overloads for different track types in 3D ---

void RenderTrack3D(const VisLineTrack &tr_in, double bounds[6])
{
    VisLineTrack tr = tr_in;
    if(tr.isLocalFrame)
    {
        CHeT::Config::ApplyTransformation(tr.x0, tr.y0, tr.z0);
        CHeT::Config::ApplyRotation(tr.ux, tr.uy, tr.uz);
    }

    double tmin = -1e30;
    double tmax = 1e30;
    bool ok = ClipLineToBox(tr.x0, tr.y0, tr.z0, tr.ux, tr.uy, tr.uz, bounds[0], bounds[1],
        bounds[2], bounds[3], bounds[4], bounds[5], tmin, tmax);
    if(ok)
    {
        TPolyLine3D *lt = new TPolyLine3D(2);
        lt->SetPoint(0, tr.z0 + tr.uz * tmin, tr.x0 + tr.ux * tmin, tr.y0 + tr.uy * tmin);
        lt->SetPoint(1, tr.z0 + tr.uz * tmax, tr.x0 + tr.ux * tmax, tr.y0 + tr.uy * tmax);
        lt->SetLineColor(tr.color);
        lt->SetLineWidth(tr.width);
        lt->SetLineStyle(tr.style);
        lt->Draw("same");
    }
}

void RenderTrack3D(const VisHelixTrack &tr, double bounds[6])
{
    int n_pts = 200; // High resolution for 3D helix
    TPolyLine3D *lt = new TPolyLine3D(n_pts);
    double dt = (tr.t_max - tr.t_min) / (n_pts - 1);

    for(int i = 0; i < n_pts; ++i)
    {
        double t = tr.t_min + i * dt;
        double x = tr.cx + tr.radius * std::cos(t);
        double y = tr.cy + tr.radius * std::sin(t);
        double z = tr.z0 + tr.dz_dt * t;

        if(tr.isLocalFrame)
        {
            CHeT::Config::ApplyTransformation(x, y, z);
        }
        // ROOT Frame: Z_root=Y_phys, etc. (Respecting original code mapping)
        lt->SetPoint(i, z, x, y);
    }
    lt->SetLineColor(tr.color);
    lt->SetLineWidth(tr.width);
    lt->SetLineStyle(tr.style);
    lt->Draw("same");
}

void RenderTrack3D(const VisGenericTrack &tr, double bounds[6])
{
    if(tr.points.empty())
        return;
    TPolyLine3D *lt = new TPolyLine3D(tr.points.size());

    for(size_t i = 0; i < tr.points.size(); ++i)
    {
        double x = tr.points[i].x;
        double y = tr.points[i].y;
        double z = tr.points[i].z;

        if(tr.isLocalFrame)
        {
            CHeT::Config::ApplyTransformation(x, y, z);
        }
        lt->SetPoint(i, z, x, y);
    }
    lt->SetLineColor(tr.color);
    lt->SetLineWidth(tr.width);
    lt->SetLineStyle(tr.style);
    lt->Draw("same");
}

// --- Core Functions ---

void Draw2DCore(const std::vector<int> &bundle_ids, const std::vector<VisPoint2D> &extraPoints,
    const std::function<void(TMultiGraph *, TMultiGraph *, TMultiGraph *, TMultiGraph *)>
        &trackDrawer,
    bool wrap_phi)
{
    // --- 1. Canvas Phi-Z (Unrolled) ---
    TCanvas *c_phiz = (TCanvas *)gROOT->FindObject("c_phiz");
    if(!c_phiz)
        c_phiz = new TCanvas("c_phiz", "Map Phi-Z", 750, 850);
    else
        c_phiz->Clear();

    c_phiz->SetLeftMargin(0.15);
    c_phiz->SetGrid();

    TMultiGraph *mg_phiz = new TMultiGraph();

    // --- 2. Canvas Lateral Projections ---
    TCanvas *c_lat = (TCanvas *)gROOT->FindObject("c_lat");
    if(!c_lat)
        c_lat = new TCanvas("c_lat", "Projections ZX / ZY", 1200, 750);
    else
        c_lat->Clear();

    TPad *pad_zx = new TPad("pad_zx", "ZX", 0.01, 0.01, 0.42, 1);
    TPad *pad_zy = new TPad("pad_zy", "ZY", 0.4, 0.15, 1, 0.8);
    pad_zx->SetLeftMargin(0.15);
    pad_zy->SetBottomMargin(0.15);
    pad_zx->Draw();
    pad_zy->Draw();

    TMultiGraph *mg_zx = new TMultiGraph();
    TMultiGraph *mg_zy = new TMultiGraph();

    // --- 3. Canvas XY View ---
    TCanvas *c_xy = (TCanvas *)gROOT->FindObject("c_xy");
    if(!c_xy)
        c_xy = new TCanvas("c_xy", "Transverse View XY", 700, 700);
    else
        c_xy->Clear();

    c_xy->SetGrid();
    TMultiGraph *mg_xy = new TMultiGraph();

    // --- Draw Bundles ---
    for(int b_id : bundle_ids)
    {
        CHeT::Config::FiberProp p = CHeT::Config::GetFiberProp(b_id);
        double d_phi_w = (CHeT::Config::BUNDLE_WIDTH / p.r) / 2.0;

        std::vector<double> vz, vphi, vx, vy;
        vz.reserve(500);
        vphi.reserve(500);
        vx.reserve(500);
        vy.reserve(500);

        // Generate points along the fiber
        for(int i = 0; i < 500; ++i)
        {
            double z = -CHeT::Config::L_HALF + i * (2.0 * CHeT::Config::L_HALF / 499.0);
            double alpha = (z + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
            double ph = p.phi0 + p.dir * alpha * M_PI;
            vz.push_back(z);
            vphi.push_back(ph);
            vx.push_back(p.r * std::cos(ph));
            vy.push_back(p.r * std::sin(ph));
        }

        // Linear Graphs (ZX, ZY, XY)
        TGraph *l_zx = new TGraph(500, &vx[0], &vz[0]);
        l_zx->SetLineColor(p.color);
        l_zx->SetLineWidth(2);
        mg_zx->Add(l_zx, "L");

        TGraph *l_zy = new TGraph(500, &vz[0], &vy[0]);
        l_zy->SetLineColor(p.color);
        l_zy->SetLineWidth(2);
        mg_zy->Add(l_zy, "L");

        TGraph *l_xy = new TGraph(500, &vx[0], &vy[0]);
        l_xy->SetLineColor(p.color);
        l_xy->SetLineWidth(1);
        mg_xy->Add(l_xy, "L");

        // Phi-Z Graphs (Handle 2pi wrapping)
        std::vector<double> sz, sf;
        for(size_t i = 0; i < 500; ++i)
        {
            double cf = wrap_phi ? CHeT::Config::wrap0_2pi(vphi[i]) : vphi[i];

            // Check for wrap-around jump
            if(wrap_phi && i > 0 && std::abs(cf - CHeT::Config::wrap0_2pi(vphi[i - 1])) > M_PI)
            {
                // Draw current segment
                TGraph *h = new TGraph();
                int nn = sz.size();
                for(int j = 0; j < nn; ++j)
                    h->SetPoint(j, sf[j] + d_phi_w, sz[j]);
                for(int j = 0; j < nn; ++j)
                    h->SetPoint(nn + j, sf[nn - 1 - j] - d_phi_w, sz[nn - 1 - j]);

                h->SetFillColorAlpha(p.color, 0.3);
                h->SetLineWidth(0);
                mg_phiz->Add(h, "F");

                TGraph *l = new TGraph(nn, &sf[0], &sz[0]);
                l->SetLineColor(p.color);
                l->SetLineWidth(2);
                mg_phiz->Add(l, "L");

                sz.clear();
                sf.clear();
            }
            sz.push_back(vz[i]);
            sf.push_back(cf);
        }

        // Draw final segment
        if(!sz.empty())
        {
            TGraph *hf = new TGraph();
            int nn = sz.size();
            for(int j = 0; j < nn; ++j)
                hf->SetPoint(j, sf[j] + d_phi_w, sz[j]);
            for(int j = 0; j < nn; ++j)
                hf->SetPoint(nn + j, sf[nn - 1 - j] - d_phi_w, sz[nn - 1 - j]);

            hf->SetFillColorAlpha(p.color, 0.3);
            hf->SetLineWidth(0);
            mg_phiz->Add(hf, "F");

            TGraph *lf = new TGraph(nn, &sf[0], &sz[0]);
            lf->SetLineColor(p.color);
            lf->SetLineWidth(2);
            mg_phiz->Add(lf, "L");
        }
    }

    // --- Draw Tracks (All Projections) ---
    if(trackDrawer)
        trackDrawer(mg_xy, mg_zx, mg_zy, mg_phiz);

    // --- Draw Extra Points (XY) ---
    for(const auto &pt : extraPoints)
    {
        TGraph *gp = new TGraph(1);
        gp->SetPoint(0, pt.x, pt.y);
        gp->SetMarkerColor(pt.color);
        gp->SetMarkerStyle(pt.markerStyle);
        gp->SetMarkerSize(pt.size);
        mg_xy->Add(gp, "P");
    }

    // --- Draw Intersections ---
    auto inters = Config::FindIntersections(bundle_ids);
    if(!inters.empty())
    {
        // XY
        TGraph *gi_xy = new TGraph();
        for(size_t i = 0; i < inters.size(); ++i)
            gi_xy->SetPoint(i, inters[i].x_loc, inters[i].y_loc);
        gi_xy->SetMarkerStyle(20);
        gi_xy->SetMarkerSize(0.8);
        gi_xy->SetMarkerColor(kBlack);
        mg_xy->Add(gi_xy, "P");

        // ZX (x vs z)
        TGraph *gi_zx = new TGraph();
        for(size_t i = 0; i < inters.size(); ++i)
            gi_zx->SetPoint(i, inters[i].x_loc, inters[i].z_loc);
        gi_zx->SetMarkerStyle(20);
        gi_zx->SetMarkerSize(0.8);
        gi_zx->SetMarkerColor(kBlack);
        mg_zx->Add(gi_zx, "P");

        // ZY (z vs y)
        TGraph *gi_zy = new TGraph();
        for(size_t i = 0; i < inters.size(); ++i)
            gi_zy->SetPoint(i, inters[i].z_loc, inters[i].y_loc);
        gi_zy->SetMarkerStyle(20);
        gi_zy->SetMarkerSize(0.8);
        gi_zy->SetMarkerColor(kBlack);
        mg_zy->Add(gi_zy, "P");

        // PhiZ (phi vs z)
        TGraph *gi_phiz = new TGraph();
        for(size_t i = 0; i < inters.size(); ++i)
        {
            double phi = std::atan2(inters[i].y_loc, inters[i].x_loc);
            if(wrap_phi && phi < 0)
                phi += 2 * M_PI;
            gi_phiz->SetPoint(i, phi, inters[i].z_loc);
        }
        gi_phiz->SetMarkerStyle(20);
        gi_phiz->SetMarkerSize(0.8);
        gi_phiz->SetMarkerColor(kBlack);
        mg_phiz->Add(gi_phiz, "P");
    }

    // --- Final Rendering ---

    // 1. Phi-Z
    c_phiz->cd();
    mg_phiz->Draw("A");
    mg_phiz->SetTitle("Detector Map #phi-z; #phi [rad]; z [mm]");
    if(wrap_phi)
    {
        mg_phiz->GetXaxis()->SetLimits(0, 2 * M_PI);
        mg_phiz->GetXaxis()->SetNdivisions(-504);

        // Custom axis labels
        TAxis *ax = mg_phiz->GetXaxis();
        ax->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
        ax->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
        ax->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
        ax->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/2");
        ax->ChangeLabel(5, -1, -1, -1, -1, -1, "2#pi");
    }

    // 2. ZX and ZY
    auto cyls = CHeT::Config::GetCylinders();

    double max_R = 0;
    for(const auto &cyl : cyls)
    {
        if(cyl.outer.radius > max_R)
            max_R = cyl.outer.radius;
    }
    if(max_R < 30.0)
        max_R = 30.0;

    double margin = 10.0;
    double limit = max_R + margin;

    pad_zx->cd();
    pad_zx->SetGrid();
    mg_zx->Draw("A");
    mg_zx->SetTitle("Top View ZX; x [mm]; z [mm]");
    mg_zx->GetXaxis()->SetLimits(-limit, limit);
    mg_zx->GetYaxis()->SetRangeUser(-CHeT::Config::L_HALF - 20, CHeT::Config::L_HALF + 20);
    for(auto &c : cyls)
    {
        TBox *b = new TBox(
            -c.inner.radius, -CHeT::Config::L_HALF, c.inner.radius, CHeT::Config::L_HALF);
        b->SetFillStyle(0);
        b->SetLineColor(kGray + 1);
        b->Draw("same");
    }

    pad_zy->cd();
    pad_zy->SetGrid();
    mg_zy->Draw("A");
    mg_zy->SetTitle("Lateral View ZY; z [mm]; y [mm]");
    mg_zy->GetXaxis()->SetLimits(-CHeT::Config::L_HALF - 20, CHeT::Config::L_HALF + 20);
    mg_zy->GetYaxis()->SetRangeUser(-limit, limit);
    for(auto &c : cyls)
    {
        TBox *b = new TBox(
            -CHeT::Config::L_HALF, -c.inner.radius, CHeT::Config::L_HALF, c.inner.radius);
        b->SetFillStyle(0);
        b->SetLineColor(kGray + 1);
        b->Draw("same");
    }

    // 3. XY
    c_xy->cd();
    mg_xy->Draw("A");
    mg_xy->SetTitle("XY View; x [mm]; y [mm]");
    mg_xy->GetXaxis()->SetLimits(-limit, limit);
    mg_xy->GetYaxis()->SetRangeUser(-limit, limit);
    for(auto &c : cyls)
    {
        TEllipse *e = new TEllipse(0, 0, c.inner.radius);
        e->SetFillStyle(0);
        e->Draw();
    }

    c_phiz->Update();
    c_lat->Update();
    c_xy->Update();
}

void Draw3DCore(const std::vector<int> &hit_ids, const std::vector<VisPoint3D> &points,
    bool drawSkeleton, const std::function<void(double[6])> &trackDrawer)
{
    // gStyle->SetCanvasPreferGL(kTRUE);

    TCanvas *c_3d = (TCanvas *)gROOT->FindObject("c_3d");
    if(!c_3d)
        c_3d = new TCanvas("c_3d", "Detector 3D View", 1200, 800);
    else
        c_3d->Clear();

    // Calculate dynamic bounding box based on geometry and transformation
    double min_x_phys = 1e9, max_x_phys = -1e9;
    double min_y_phys = 1e9, max_y_phys = -1e9;
    double min_z_phys = 1e9, max_z_phys = -1e9;

    auto cylinders = CHeT::Config::GetCylinders();
    double max_R = 0;
    for(const auto &cyl : cylinders)
    {
        if(cyl.outer.radius > max_R)
            max_R = cyl.outer.radius;
    }
    if(max_R < 1.0)
        max_R = 100.0; // Fallback

    double corners[8][3]
        = { { max_R, max_R, CHeT::Config::L_HALF }, { max_R, -max_R, CHeT::Config::L_HALF },
              { -max_R, max_R, CHeT::Config::L_HALF }, { -max_R, -max_R, CHeT::Config::L_HALF },
              { max_R, max_R, -CHeT::Config::L_HALF }, { max_R, -max_R, -CHeT::Config::L_HALF },
              { -max_R, max_R, -CHeT::Config::L_HALF }, { -max_R, -max_R, -CHeT::Config::L_HALF } };

    for(int i = 0; i < 8; ++i)
    {
        double x = corners[i][0];
        double y = corners[i][1];
        double z = corners[i][2];
        CHeT::Config::ApplyTransformation(x, y, z);

        if(x < min_x_phys)
            min_x_phys = x;
        if(x > max_x_phys)
            max_x_phys = x;
        if(y < min_y_phys)
            min_y_phys = y;
        if(y > max_y_phys)
            max_y_phys = y;
        if(z < min_z_phys)
            min_z_phys = z;
        if(z > max_z_phys)
            max_z_phys = z;
    }

    double margin = 20.0;
    min_x_phys -= margin;
    max_x_phys += margin;
    min_y_phys -= margin;
    max_y_phys += margin;
    min_z_phys -= margin;
    max_z_phys += margin;

    // Force aspect ratio to be isotropic (Cube) only if rotated
    double rx, ry, rz;
    CHeT::Config::GetRotation(rx, ry, rz);
    bool isRotated = (std::abs(rx) > 1e-9 || std::abs(ry) > 1e-9 || std::abs(rz) > 1e-9);

    if(isRotated)
    {
        double dx = max_x_phys - min_x_phys;
        double dy = max_y_phys - min_y_phys;
        double dz = max_z_phys - min_z_phys;
        double max_dim = std::max({ dx, dy, dz });

        double cx = (min_x_phys + max_x_phys) / 2.0;
        double cy = (min_y_phys + max_y_phys) / 2.0;
        double cz = (min_z_phys + max_z_phys) / 2.0;

        min_x_phys = cx - max_dim / 2.0;
        max_x_phys = cx + max_dim / 2.0;
        min_y_phys = cy - max_dim / 2.0;
        max_y_phys = cy + max_dim / 2.0;
        min_z_phys = cz - max_dim / 2.0;
        max_z_phys = cz + max_dim / 2.0;
    }

    double bounds[6] = { min_x_phys, max_x_phys, min_y_phys, max_y_phys, min_z_phys, max_z_phys };

    // ROOT Frame: X->Z_phys, Y->X_phys, Z->Y_phys
    TH3F *h_frame = new TH3F("h_frame", "; Z [mm]; X [mm]; Y [mm]", 1, min_z_phys, max_z_phys, 1,
        min_x_phys, max_x_phys, 1, min_y_phys, max_y_phys);
    h_frame->SetDirectory(0); // Detach from directory to prevent ROOT auto-deletion issues
    h_frame->SetStats(0);
    h_frame->GetXaxis()->SetTitleOffset(1.5);
    h_frame->GetYaxis()->SetTitleOffset(1.5);
    h_frame->GetZaxis()->SetTitleOffset(1.);
    h_frame->Draw();

    // --- Draw Skeleton (Inactive fibers) ---
    if(drawSkeleton)
    {
        auto cylinders = CHeT::Config::GetCylinders();
        for(const auto &cyl : cylinders)
        {
            const CHeT::Config::LayerConfig *layers[2] = { &cyl.inner, &cyl.outer };
            for(int l = 0; l < 2; ++l)
            {
                for(int b = 0; b < layers[l]->nBundles; ++b)
                {
                    int global_id = CHeT::Config::GetGlobalIdFromGeometry(cyl.id, l, b);
                    CHeT::Config::FiberProp p = CHeT::Config::GetFiberProp(global_id);

                    // Ghost fiber
                    TPolyLine3D *bg_f = new TPolyLine3D(10);
                    for(int i = 0; i < 10; ++i)
                    {
                        double z = -CHeT::Config::L_HALF + i * (2.0 * CHeT::Config::L_HALF / 9.0);
                        double a = (z + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
                        double ph = p.phi0 + p.dir * a * M_PI;

                        double x3 = p.r * std::cos(ph);
                        double y3 = p.r * std::sin(ph);
                        double z3 = z;
                        CHeT::Config::ApplyTransformation(x3, y3, z3);

                        // Map to ROOT: X_root=z, Y_root=x, Z_root=y
                        bg_f->SetPoint(i, z3, x3, y3);
                    }
                    bg_f->SetLineColorAlpha(p.color, 0.1); // Very transparent
                    bg_f->Draw("same");
                }
            }
        }
    }

    // --- Draw Hits (Active fibers) ---
    for(int id : hit_ids)
    {
        if(id < 0)
            continue;
        CHeT::Config::FiberProp p = CHeT::Config::GetFiberProp(id);

        TPolyLine3D *fl = new TPolyLine3D(50); // More detailed for active ones
        for(int i = 0; i < 50; ++i)
        {
            double z = -CHeT::Config::L_HALF + i * (2.0 * CHeT::Config::L_HALF / 49.0);
            double a = (z + CHeT::Config::L_HALF) / (2.0 * CHeT::Config::L_HALF);
            double ph = p.phi0 + p.dir * a * M_PI;

            double x3 = p.r * std::cos(ph);
            double y3 = p.r * std::sin(ph);
            double z3 = z;
            CHeT::Config::ApplyTransformation(x3, y3, z3);

            // Map to ROOT: X_root=z, Y_root=x, Z_root=y
            fl->SetPoint(i, z3, x3, y3);
        }
        fl->SetLineColor(p.color);
        fl->SetLineWidth(3);
        fl->Draw("same");
    }

    // --- Draw Tracks ---
    if(trackDrawer)
        trackDrawer(bounds);

    // --- Draw 3D Points ---
    for(const auto &pt_in : points)
    {
        VisPoint3D pt = pt_in;
        if(pt.isLocalFrame)
        {
            CHeT::Config::ApplyTransformation(pt.x, pt.y, pt.z);
        }

        TPolyMarker3D *pm = new TPolyMarker3D(1);
        // Note: Mapping Physics (x,y,z) to ROOT 3D Frame (Z, X, Y)
        // ROOT X axis = Physics Z
        // ROOT Y axis = Physics X
        // ROOT Z axis = Physics Y
        pm->SetPoint(0, pt.z, pt.x, pt.y);

        pm->SetMarkerColor(pt.color);
        pm->SetMarkerStyle(pt.markerStyle);
        pm->SetMarkerSize(pt.size);
        pm->Draw("same");
    }

    c_3d->Update();
}

} // namespace Vis
} // namespace CHeT
