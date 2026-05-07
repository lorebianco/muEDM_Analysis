#include <cmath>
#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TRandom.h>
#include <TString.h>
#include <TTree.h>

#include "CHeT/CHeTGlobalSettings.hh"
#include "ToyGenerator.hh"

using namespace ftxui;

struct AppConfig
{
    AppConfig()
    {
        CHeT::Config::GetTranslation(tx, ty, tz);
        CHeT::Config::GetRotation(rx, ry, rz);
    }

    std::string mode = "cosmic";
    int nEvents = 1000;
    std::string outFile = "toy_data.root";
    double efficiency = 1.0;

    double tx = 0.0, ty = 0.0, tz = 0.0;
    double rx = 0.0, ry = 0.0, rz = 0.0;
    double offset_exp = CHeT::Config::GetOffsetExp();
    std::vector<double> deltas = CHeT::Config::GetDeltas();

    bool active_cyls[6] = { true, true, true, true, true, true };

    ToyMC::CosmicDetConfig cosmicDet;
};

void PrintUsage(const char *progName)
{
    std::cout << "Usage: " << progName << " <mode> <n_events> <output_file> [efficiency]\n"
              << "  mode:        cosmic | michel\n"
              << "  n_events:    number of events to generate\n"
              << "  output_file: output ROOT file (e.g., toy_data.root)\n"
              << "  efficiency:  (optional) detector efficiency [0.0 - 1.0], default 1.0\n"
              << "\nRun without arguments to launch the Interactive GUI!\n";
}

template <typename T> std::string format_val(T val)
{
    std::ostringstream out;
    out << val;
    return out.str();
}

void RunInteractiveGUI(AppConfig &config)
{
    auto screen = ScreenInteractive::TerminalOutput();

    std::vector<std::string> mode_entries = { "cosmic", "michel" };
    int mode_selected = (config.mode == "michel") ? 1 : 0;

    std::string ev_str = format_val(config.nEvents);
    std::string out_str = config.outFile;
    std::string eff_str = format_val(config.efficiency);

    std::string tx_str = format_val(config.tx);
    std::string ty_str = format_val(config.ty);
    std::string tz_str = format_val(config.tz);

    std::string rx_str = format_val(config.rx);
    std::string ry_str = format_val(config.ry);
    std::string rz_str = format_val(config.rz);

    std::string off_str = format_val(config.offset_exp);
    std::vector<std::string> d_str(6);
    for(int i = 0; i < 6; ++i)
        d_str[i] = format_val(config.deltas.size() > i ? config.deltas[i] : 0.0);

    std::string xup_str = format_val(config.cosmicDet.xUp);
    std::string yup_str = format_val(config.cosmicDet.yUp);
    std::string zup_str = format_val(config.cosmicDet.zUp);
    std::string hxup_str = format_val(config.cosmicDet.halfX_up);
    std::string hzup_str = format_val(config.cosmicDet.halfZ_up);

    std::string xdn_str = format_val(config.cosmicDet.xDown);
    std::string ydn_str = format_val(config.cosmicDet.yDown);
    std::string zdn_str = format_val(config.cosmicDet.zDown);
    std::string hxdn_str = format_val(config.cosmicDet.halfX_down);
    std::string hzdn_str = format_val(config.cosmicDet.halfZ_down);

    InputOption input_opt;
    input_opt.transform = [](InputState state)
    {
        state.element |= state.focused ? bgcolor(Color::Blue) | color(Color::White) : nothing;
        return state.element;
    };

    Component radio_mode = Radiobox(&mode_entries, &mode_selected);
    Component inp_ev = Input(&ev_str, "1000", input_opt);
    Component inp_out = Input(&out_str, "toy_data.root", input_opt);
    Component inp_eff = Input(&eff_str, "1.0", input_opt);

    Component inp_tx = Input(&tx_str, "0.0", input_opt);
    Component inp_ty = Input(&ty_str, "0.0", input_opt);
    Component inp_tz = Input(&tz_str, "0.0", input_opt);
    Component inp_rx = Input(&rx_str, "0.0", input_opt);
    Component inp_ry = Input(&ry_str, "0.0", input_opt);
    Component inp_rz = Input(&rz_str, "0.0", input_opt);
    Component inp_off = Input(&off_str, "0.0", input_opt);
    std::vector<Component> inp_d(6);
    for(int i = 0; i < 6; ++i)
        inp_d[i] = Input(&d_str[i], "0.0", input_opt);

    Component inp_xup = Input(&xup_str, "0.0", input_opt);
    Component inp_yup = Input(&yup_str, "386.0", input_opt);
    Component inp_zup = Input(&zup_str, "0.0", input_opt);
    Component inp_hxup = Input(&hxup_str, "90.0", input_opt);
    Component inp_hzup = Input(&hzup_str, "90.0", input_opt);

    Component inp_xdn = Input(&xdn_str, "0.0", input_opt);
    Component inp_ydn = Input(&ydn_str, "-476.0", input_opt);
    Component inp_zdn = Input(&zdn_str, "0.0", input_opt);
    Component inp_hxdn = Input(&hxdn_str, "90.0", input_opt);
    Component inp_hzdn = Input(&hzdn_str, "90.0", input_opt);

    Component cb_cyl0 = Checkbox("0", &config.active_cyls[0]);
    Component cb_cyl1 = Checkbox("1", &config.active_cyls[1]);
    Component cb_cyl2 = Checkbox("2", &config.active_cyls[2]);
    Component cb_cyl3 = Checkbox("3", &config.active_cyls[3]);
    Component cb_cyl4 = Checkbox("4", &config.active_cyls[4]);
    Component cb_cyl5 = Checkbox("5", &config.active_cyls[5]);

    bool proceed = false;
    Component btn_run = Button(
        "Run Simulation",
        [&]
        {
            proceed = true;
            screen.ExitLoopClosure()();
        },
        ButtonOption::Animated());
    Component btn_exit = Button("Exit", screen.ExitLoopClosure(), ButtonOption::Animated());

    // --- LOGICA DEI CONTAINER (Riorganizzata) ---

    // 1. Gruppo di sinistra (Generation Settings)
    auto container_settings_left = Container::Vertical({
        radio_mode,
        inp_ev,
        inp_out,
        inp_eff,
    });

    // 2. Gruppo di destra (Detector Geometry)
    // Creiamo una riga orizzontale per i checkbox
    auto container_cylinders
        = Container::Horizontal({ cb_cyl0, cb_cyl1, cb_cyl2, cb_cyl3, cb_cyl4, cb_cyl5 });

    auto container_geometry_right = Container::Vertical({
        container_cylinders, // Ora i checkbox sono un sottogruppo orizzontale
        Container::Horizontal({ inp_tx, inp_ty, inp_tz }),
        Container::Horizontal({ inp_rx, inp_ry, inp_rz }),
        Container::Horizontal(
            { inp_off, inp_d[0], inp_d[1], inp_d[2], inp_d[3], inp_d[4], inp_d[5] }),
    });

    // 3. Uniamo i due blocchi superiori in una riga
    auto container_top_row = Container::Horizontal({
        container_settings_left,
        container_geometry_right,
    });

    // 4. Blocco inferiore (Cosmic Detectors)
    auto container_cosmic_bottom = Container::Vertical({
        Container::Horizontal({ inp_xup, inp_yup, inp_zup, inp_hxup, inp_hzup }),
        Container::Horizontal({ inp_xdn, inp_ydn, inp_zdn, inp_hxdn, inp_hzdn }),
    });

    // 5. Bottoni
    auto container_buttons = Container::Horizontal({ btn_run, btn_exit });

    // 6. Container PRINCIPALE (Verticale)
    auto main_container = Container::Vertical({
        container_top_row,
        container_cosmic_bottom,
        container_buttons,
    });

    // --- RENDERER ---
    auto renderer = Renderer(main_container,
        [&]
        {
            return vbox({
                text(" muEDM ToyMC Generator ") | bold | center,
                separator(),
                hbox({
                    // Colonna Sinistra
                    vbox({
                        text("--- Generation Settings ---") | bold,
                        hbox(text("Mode:       "), radio_mode->Render()),
                        hbox(text("Events:     "), inp_ev->Render() | size(WIDTH, EQUAL, 15)),
                        hbox(text("Out File:   "), inp_out->Render() | size(WIDTH, EQUAL, 25)),
                        hbox(text("Efficiency: "), inp_eff->Render() | size(WIDTH, EQUAL, 10)),
                    }) | border
                        | flex,

                    // Colonna Destra
                    vbox({
                        text("--- Detector Geometry (Misalignment) ---") | bold,
                        hbox({ text("Active Cyls: "), cb_cyl0->Render(), text(" "),
                            cb_cyl1->Render(), text(" "), cb_cyl2->Render(), text(" "),
                            cb_cyl3->Render(), text(" "), cb_cyl4->Render(), text(" "),
                            cb_cyl5->Render() }),
                        hbox(text("Trans (X,Y,Z): "), inp_tx->Render() | size(WIDTH, EQUAL, 6),
                            text(" "), inp_ty->Render() | size(WIDTH, EQUAL, 6), text(" "),
                            inp_tz->Render() | size(WIDTH, EQUAL, 6)),
                        hbox(text("Rot (X,Y,Z):   "), inp_rx->Render() | size(WIDTH, EQUAL, 6),
                            text(" "), inp_ry->Render() | size(WIDTH, EQUAL, 6), text(" "),
                            inp_rz->Render() | size(WIDTH, EQUAL, 6)),
                        hbox(text("Offset (Exp):  "), inp_off->Render() | size(WIDTH, EQUAL, 6)),
                        hbox(text("Deltas (0,1,2):"), inp_d[0]->Render() | size(WIDTH, EQUAL, 6),
                            text(" "), inp_d[1]->Render() | size(WIDTH, EQUAL, 6), text(" "),
                            inp_d[2]->Render() | size(WIDTH, EQUAL, 6)),
                        hbox(text("Deltas (3,4,5):"), inp_d[3]->Render() | size(WIDTH, EQUAL, 6),
                            text(" "), inp_d[4]->Render() | size(WIDTH, EQUAL, 6), text(" "),
                            inp_d[5]->Render() | size(WIDTH, EQUAL, 6)),
                    }) | border
                        | flex,
                }),

                // Blocco Cosmic
                vbox({
                    text("--- Cosmic Detectors (Only for cosmic mode) ---") | bold,
                    hbox(text("Up Pos(X,Y,Z): "), inp_xup->Render() | size(WIDTH, EQUAL, 7),
                        text(" "), inp_yup->Render() | size(WIDTH, EQUAL, 7), text(" "),
                        inp_zup->Render() | size(WIDTH, EQUAL, 7), text("  Half(X,Z): "),
                        inp_hxup->Render() | size(WIDTH, EQUAL, 7), text(" "),
                        inp_hzup->Render() | size(WIDTH, EQUAL, 7)),
                    hbox(text("Dn Pos(X,Y,Z): "), inp_xdn->Render() | size(WIDTH, EQUAL, 7),
                        text(" "), inp_ydn->Render() | size(WIDTH, EQUAL, 7), text(" "),
                        inp_zdn->Render() | size(WIDTH, EQUAL, 7), text("  Half(X,Z): "),
                        hxdn_str == "" ? inp_hxdn->Render()
                                       : inp_hxdn->Render() | size(WIDTH, EQUAL, 7),
                        text(" "), inp_hzdn->Render() | size(WIDTH, EQUAL, 7)),
                }) | border,

                separator(),
                hbox({ btn_run->Render(), text("   "), btn_exit->Render() }) | center,
            });
        });

    screen.Loop(renderer);

    std::cout << "\x1b[?25h\x1b[0 q\x1b[0m\x1b[2J\x1b[H" << std::flush;

    if(!proceed)
    {
        exit(0);
    }

    // Parse outputs back
    config.mode = mode_entries[mode_selected];

    auto parse_double = [](const std::string &s, double def)
    {
        try
        {
            return std::stod(s);
        }
        catch(...)
        {
            return def;
        }
    };
    auto parse_int = [](const std::string &s, int def)
    {
        try
        {
            return std::stoi(s);
        }
        catch(...)
        {
            return def;
        }
    };

    config.nEvents = parse_int(ev_str, 1000);
    config.outFile = out_str;
    config.efficiency = parse_double(eff_str, 1.0);

    config.tx = parse_double(tx_str, 0.0);
    config.ty = parse_double(ty_str, 0.0);
    config.tz = parse_double(tz_str, 0.0);
    config.rx = parse_double(rx_str, 0.0);
    config.ry = parse_double(ry_str, 0.0);
    config.rz = parse_double(rz_str, 0.0);
    config.offset_exp = parse_double(off_str, 0.0);
    config.deltas.resize(6);
    for(int i = 0; i < 6; ++i)
        config.deltas[i] = parse_double(d_str[i], 0.0);

    config.cosmicDet.xUp = parse_double(xup_str, 0.0);
    config.cosmicDet.yUp = parse_double(yup_str, 386.0);
    config.cosmicDet.zUp = parse_double(zup_str, 0.0);
    config.cosmicDet.halfX_up = parse_double(hxup_str, 90.0);
    config.cosmicDet.halfZ_up = parse_double(hzup_str, 90.0);

    config.cosmicDet.xDown = parse_double(xdn_str, 0.0);
    config.cosmicDet.yDown = parse_double(ydn_str, -476.0);
    config.cosmicDet.zDown = parse_double(zdn_str, 0.0);
    config.cosmicDet.halfX_down = parse_double(hxdn_str, 90.0);
    config.cosmicDet.halfZ_down = parse_double(hzdn_str, 90.0);
}

int main(int argc, char **argv)
{
    AppConfig config;

    if(argc == 1)
    {
        RunInteractiveGUI(config);
    }
    else if(argc >= 4)
    {
        config.mode = argv[1];
        config.nEvents = std::stoi(argv[2]);
        config.outFile = argv[3];
        if(argc >= 5)
        {
            config.efficiency = std::stod(argv[4]);
        }
    }
    else
    {
        PrintUsage(argv[0]);
        return 1;
    }

    if(config.mode != "cosmic" && config.mode != "michel")
    {
        std::cerr << "Error: Unknown mode '" << config.mode << "'. Must be 'cosmic' or 'michel'.\n";
        return 1;
    }

    // Apply geometry misalignments to the global CHeT settings
    std::vector<int> act_cyls;
    for(int i = 0; i < 6; ++i)
    {
        if(config.active_cyls[i])
            act_cyls.push_back(i);
    }
    CHeT::Config::SetActiveCylinders(act_cyls);

    CHeT::Config::SetTranslation(config.tx, config.ty, config.tz);
    CHeT::Config::SetRotation(config.rx, config.ry, config.rz);
    CHeT::Config::SetOffsetExp(config.offset_exp);
    CHeT::Config::SetDeltas(config.deltas);

    std::string outFile = config.outFile;
    if(outFile.find('/') == std::string::npos)
    {
        outFile = "../../data/input/toy/" + outFile;
    }

    std::cout << "--- muEDM ToyMC Generator ---\n"
              << "Mode:       " << config.mode << "\n"
              << "Events:     " << config.nEvents << "\n"
              << "Output:     " << outFile << "\n"
              << "Efficiency: " << config.efficiency << "\n"
              << "-----------------------------\n";

    // Create the output file
    TFile *fOut = TFile::Open(outFile.c_str(), "RECREATE");
    if(!fOut || fOut->IsZombie())
    {
        std::cerr << "Error: Could not open output file " << outFile << "\n";
        return 1;
    }

    // =========================================================================
    // 1. CHeT Tree (Detector Data - High Level)
    // =========================================================================
    TTree *treeCHeT = new TTree("chet", "Single Event Detector Data");

    int eventID = 0;
    std::vector<int> all_bundle;

    treeCHeT->Branch("EventID", &eventID);
    treeCHeT->Branch("All_Bundle", &all_bundle);

    // =========================================================================
    // 2. SIM Tree (MonteCarlo Truth Data)
    // =========================================================================
    TTree *treeSIM = new TTree("sim", "MonteCarlo Truth Data");

    treeSIM->Branch("EventID", &eventID);

    std::vector<double> true_hit_x, true_hit_y, true_hit_z;
    treeSIM->Branch("mc_hits_x", &true_hit_x);
    treeSIM->Branch("mc_hits_y", &true_hit_y);
    treeSIM->Branch("mc_hits_z", &true_hit_z);

    int trackID = 0;
    int particleID = 0; // e.g., 11 for electron, 0 for cosmic
    treeSIM->Branch("TrackID", &trackID);
    treeSIM->Branch("ParticleID", &particleID);

    // Common Truth Variables (can be padded with 0 if not used in a specific mode)
    // Track Variables (Mathematical parametrization for fit comparison)
    double trk_x0 = 0, trk_y0 = 0, trk_z0 = 0;
    double trk_ux = 0, trk_uy = 0, trk_uz = 0;
    double trk_sx = 0, trk_sz = 0;
    double trk_R = 0, trk_cx = 0, trk_cy = 0, trk_tmin = 0, trk_tmax = 0;

    // MC Truth Variables (Physical generation points and sanity checks)
    double mc_x = 0, mc_y = 0, mc_z = 0;
    double mc_px = 0, mc_py = 0, mc_pz = 0;
    double mc_E = 0;

    treeSIM->Branch("mc_x", &mc_x);
    treeSIM->Branch("mc_y", &mc_y);
    treeSIM->Branch("mc_z", &mc_z);
    treeSIM->Branch("mc_px", &mc_px);
    treeSIM->Branch("mc_py", &mc_py);
    treeSIM->Branch("mc_pz", &mc_pz);

    if(config.mode == "cosmic")
    {
        // 50% positive (13 is muon-) and 50% negative (13 is muon-) -> wait, PDG for mu- is 13, mu+
        // is -13
        particleID = (gRandom->Rndm() > 0.5) ? 13 : -13;
        treeSIM->Branch("mc_E", &mc_E);
        treeSIM->Branch("trk_x0", &trk_x0);
        treeSIM->Branch("trk_z0", &trk_z0);
        treeSIM->Branch("trk_sx", &trk_sx);
        treeSIM->Branch("trk_sz", &trk_sz);
    }
    else if(config.mode == "michel")
    {
        particleID = -11; // Michel Positron
        treeSIM->Branch("mc_E", &mc_E);
        treeSIM->Branch("trk_R", &trk_R);
        treeSIM->Branch("trk_cx", &trk_cx);
        treeSIM->Branch("trk_cy", &trk_cy);
        treeSIM->Branch("trk_z0", &trk_z0);
        treeSIM->Branch("trk_uz", &trk_uz);
        treeSIM->Branch("trk_tmin", &trk_tmin);
        treeSIM->Branch("trk_tmax", &trk_tmax);
    }

    std::cout << "Generating events...\n";
    for(eventID = 0; eventID < config.nEvents; ++eventID)
    {
        trackID = eventID;
        if(config.mode == "cosmic")
        {
            particleID = (gRandom->Rndm() > 0.5) ? 13 : -13; // 50% mu-, 50% mu+
        }

        all_bundle.clear();
        true_hit_x.clear();
        true_hit_y.clear();
        true_hit_z.clear();

        if(config.mode == "cosmic")
        {
            ToyMC::CosmicTrack tr = ToyMC::GenerateCosmic(config.cosmicDet);
            mc_x = tr.x0;
            mc_y = tr.y0;
            mc_z = tr.z0;
            mc_E = tr.E_kin;

            // Calculate momentum using relativistic formula
            double p = std::sqrt(tr.E_kin * tr.E_kin + 2.0 * 105.658 * tr.E_kin);
            mc_px = p * tr.ux;
            mc_py = p * tr.uy;
            mc_pz = p * tr.uz;

            double t_to_y0 = -tr.y0 / tr.uy;
            trk_x0 = tr.x0 + tr.ux * t_to_y0;
            trk_z0 = tr.z0 + tr.uz * t_to_y0;
            trk_sx = tr.ux / tr.uy;
            trk_sz = tr.uz / tr.uy;

            ToyMC::HitResult res = ToyMC::FindCosmicHits(tr, config.efficiency);
            all_bundle = res.bundles;
            true_hit_x = res.x;
            true_hit_y = res.y;
            true_hit_z = res.z;
        }
        else if(config.mode == "michel")
        {
            ToyMC::MichelTrack tr = ToyMC::GenerateMichelTrack(false);
            mc_x = tr.x0;
            mc_y = tr.y0;
            mc_z = 0.0;
            mc_E = tr.E_kin;

            double p = std::sqrt(tr.E_kin * tr.E_kin + 2.0 * 0.511 * tr.E_kin);
            double pt = p * std::sin(tr.theta_rad);
            mc_px = pt * std::cos(tr.phi_dir);
            mc_py = pt * std::sin(tr.phi_dir);
            mc_pz = p * std::cos(tr.theta_rad);

            trk_R = tr.radius;
            trk_cx = tr.cx;
            trk_cy = tr.cy;
            trk_z0 = tr.z0;
            trk_uz = tr.dz_dt;
            trk_tmin = tr.t_min;
            trk_tmax = tr.t_max;
            ToyMC::HitResult res = ToyMC::FindMichelHits(tr, config.efficiency);
            all_bundle = res.bundles;
            true_hit_x = res.x;
            true_hit_y = res.y;
            true_hit_z = res.z;
        }

        // Save parallel entries to both trees
        treeCHeT->Fill();
        treeSIM->Fill();

        if((eventID + 1) % 1000 == 0)
        {
            std::cout << "Processed " << (eventID + 1) << " / " << config.nEvents << " events\r"
                      << std::flush;
        }
    }
    std::cout << "\nGeneration complete. Writing to disk...\n";

    fOut->cd();

    // =========================================================================
    // Save Geometry configuration to UserInfo for propagation
    // =========================================================================
    auto save_metadata = [&](TTree *t)
    {
        t->GetUserInfo()->Add(new TNamed("Geometry_Version", "Nominal_With_Misalignments"));
        t->GetUserInfo()->Add(new TParameter<double>("Efficiency", config.efficiency));

        t->GetUserInfo()->Add(new TParameter<double>("Geom_Tx", config.tx));
        t->GetUserInfo()->Add(new TParameter<double>("Geom_Ty", config.ty));
        t->GetUserInfo()->Add(new TParameter<double>("Geom_Tz", config.tz));
        t->GetUserInfo()->Add(new TParameter<double>("Geom_Rx", config.rx));
        t->GetUserInfo()->Add(new TParameter<double>("Geom_Ry", config.ry));
        t->GetUserInfo()->Add(new TParameter<double>("Geom_Rz", config.rz));

        t->GetUserInfo()->Add(new TParameter<double>("Geom_OffsetExp", config.offset_exp));
        for(size_t i = 0; i < config.deltas.size(); ++i)
        {
            t->GetUserInfo()->Add(
                new TParameter<double>(Form("Geom_Delta_%zu", i), config.deltas[i]));
        }

        auto active_cyls = CHeT::Config::GetActiveCylinders();
        t->GetUserInfo()->Add(new TParameter<int>("ActiveCyl_Count", active_cyls.size()));
        for(size_t i = 0; i < active_cyls.size(); ++i)
        {
            t->GetUserInfo()->Add(new TParameter<int>(Form("ActiveCyl_%zu", i), active_cyls[i]));
        }
    };

    save_metadata(treeCHeT);
    save_metadata(treeSIM);

    treeCHeT->Write();
    treeSIM->Write();
    fOut->Close();

    std::cout << "Saved successfully to " << outFile << "\n";

    return 0;
}
