#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#include <TApplication.h>
#include <TEveManager.h>
#include <TROOT.h>

#include "CLI11.hpp"
#include "auxiliaryalgorithms.hh"
#include "config.hh"
#include "fitteralgorithms.hh"

using namespace std;

template <typename T> std::string format_val(T val)
{
    std::ostringstream out;
    out << val;
    return out.str();
}

void RunInteractiveGUI(Config &config)
{
    using namespace ftxui;

    // --- 1. MAPPATURA ENUM ---
    std::vector<FitterType> fitter_values
        = { FitterType::CosmicFitter, FitterType::SpacepointFitter, FitterType::HelixFitter };
    std::vector<muEDM::Fields::FieldBType> bfield_values = { muEDM::Fields::FieldBType::MultiMap,
        muEDM::Fields::FieldBType::SingleMap, muEDM::Fields::FieldBType::Constant };

    // --- 2. DICHIARAZIONE VARIABILI ---
    auto screen = ScreenInteractive::TerminalOutput();

    std::string input_dir = config.inputDir;
    std::string input_file = config.inputFileName;
    std::string output_dir = config.outputDir;
    std::string geom_dir = config.geomDir;
    std::string geom_file = config.geomFileName;
    std::string b_dir = config.bFieldDir;
    std::string b_path = config.bFieldMapName;
    std::string b_scale_str = format_val(config.bScale);
    std::string b_const_str = format_val(config.bConstZ);
    std::string entry_max_str = format_val(config.rangeLoop.second);
    std::string turns_str = format_val(config.turnID);
    std::string toa_min_str = format_val(config.cutToAMin);
    std::string toa_max_str = format_val(config.cutToAMax);
    std::string tot_min_str = format_val(config.cutToTMin);
    std::string tot_max_str = format_val(config.cutToTMax);
    std::string tx_str = format_val(config.geom_tx);
    std::string ty_str = format_val(config.geom_ty);
    std::string tz_str = format_val(config.geom_tz);
    std::string rx_str = format_val(config.geom_rx);
    std::string ry_str = format_val(config.geom_ry);
    std::string rz_str = format_val(config.geom_rz);
    std::string off_str = format_val(config.geom_offset);

    std::vector<std::string> d_str(6);
    for(int i = 0; i < 6; ++i)
        d_str[i] = format_val(config.geom_deltas.size() > i ? config.geom_deltas[i] : 0.0);

    std::string event_str = format_val(config.event);
    std::string start_str = format_val(config.rangeLoop.first);

    int selected_runmode = config.processSingle ? (config.event == -1 ? 1 : 0) : 2;
    int selected_fitter = 0;
    for(size_t i = 0; i < fitter_values.size(); ++i)
        if(config.fitter == fitter_values[i])
            selected_fitter = (int)i;
    int selected_btype = 0;
    for(size_t i = 0; i < bfield_values.size(); ++i)
        if(config.bFieldType == bfield_values[i])
            selected_btype = (int)i;

    // --- 3. CREAZIONE COMPONENTI ---
    std::vector<std::string> fitter_entries = { "Cosmic", "Spacepoint", "Helix" };
    std::vector<std::string> bfield_entries = { "MultiMap", "SingleMap", "Constant" };
    std::vector<std::string> runmode_entries = { "Single Event", "Random Event", "Event Loop" };

    InputOption input_opt;
    input_opt.transform = [](InputState state)
    {
        state.element |= state.focused ? bgcolor(Color::Blue) | color(Color::White) : nothing;
        return state.element;
    };

    Component inp_input_dir = Input(&input_dir, "path/to/input/", input_opt);
    Component inp_input = Input(&input_file, "input.root", input_opt);
    Component inp_output_dir = Input(&output_dir, "path/to/output/", input_opt);
    Component inp_geom_dir = Input(&geom_dir, "path/to/geom/", input_opt);
    Component inp_geom = Input(&geom_file, "geom.gdml", input_opt);
    Component menu_btype = Radiobox(&bfield_entries, &selected_btype);
    Component inp_bdir = Input(&b_dir, "path/to/maps/", input_opt);
    Component inp_bpath = Input(&b_path, "map_folder", input_opt);
    Component inp_bscale = Input(&b_scale_str, "-1.0", input_opt);
    Component inp_bconst = Input(&b_const_str, "28.9", input_opt);
    Component menu_fitter = Radiobox(&fitter_entries, &selected_fitter);
    Component cb_prefitter = Checkbox("Use Pre-fitter", &config.usePrefitter);
    Component cb_quiet = Checkbox("Quiet Mode", &config.quietMode);
    Component cb_smearing = Checkbox("Apply Smearing", &config.useSmearing);
    Component cb_pttrec = Checkbox("Pattern Recognition", &config.pttrecMode);
    Component cb_over_geom = Checkbox("Override Nominal Geometry", &config.overrideGeom);
    Component cb_mc_geom = Checkbox("Use True MC Geometry", &config.useTrueMCGeom);
    Component cb_mc_hits = Checkbox("Use True MC Hits", &config.useTrueMCHits);

    // FIX SPAZIATURA CHECKBOX
    Component cb_cyls = Container::Horizontal(
        { Checkbox("0 ", &config.active_cyls[0]), Checkbox("1 ", &config.active_cyls[1]),
            Checkbox("2 ", &config.active_cyls[2]), Checkbox("3 ", &config.active_cyls[3]),
            Checkbox("4 ", &config.active_cyls[4]), Checkbox("5 ", &config.active_cyls[5]) });

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

    Component menu_runmode = Radiobox(&runmode_entries, &selected_runmode);
    Component inp_event = Input(&event_str, "0", input_opt);
    Component inp_start = Input(&start_str, "0", input_opt);
    Component inp_entrymax = Input(&entry_max_str, "1000", input_opt);
    Component inp_turns = Input(&turns_str, "0", input_opt);
    Component inp_toa_min = Input(&toa_min_str, "215", input_opt);
    Component inp_toa_max = Input(&toa_max_str, "250", input_opt);
    Component inp_tot_min = Input(&tot_min_str, "40", input_opt);
    Component inp_tot_max = Input(&tot_max_str, "220", input_opt);

    bool should_exit = false;
    Component btn_run = Button(
        "  RUN ANALYSIS  ", [&] { screen.Exit(); }, ButtonOption::Animated());
    Component btn_exit = Button(
        "  EXIT  ",
        [&]
        {
            should_exit = true;
            screen.Exit();
        },
        ButtonOption::Animated());

    // --- 4. ORGANIZZAZIONE LAYOUT ---
    auto geom_inputs = Container::Vertical({
        cb_cyls,
        Container::Horizontal({ inp_tx, inp_ty, inp_tz }),
        Container::Horizontal({ inp_rx, inp_ry, inp_rz }),
        Container::Horizontal({ inp_off, inp_d[0], inp_d[1], inp_d[2] }),
        Container::Horizontal({ inp_d[3], inp_d[4], inp_d[5] }),
    });
    auto maybe_geom = Maybe(geom_inputs, &config.overrideGeom);

    auto col1_layout
        = Container::Vertical({ inp_input_dir, inp_input, inp_output_dir, inp_geom_dir, inp_geom,
            menu_btype, inp_bdir, inp_bpath, Container::Horizontal({ inp_bscale, inp_bconst }),
            menu_fitter, Container::Horizontal({ cb_prefitter, cb_quiet }) });

    auto col2_layout = Container::Vertical({ Container::Horizontal({ inp_toa_min, inp_toa_max }),
        Container::Horizontal({ inp_tot_min, inp_tot_max }), inp_turns, menu_runmode, inp_event,
        Container::Horizontal({ inp_start, inp_entrymax }), cb_over_geom, maybe_geom, cb_pttrec,
        cb_smearing, cb_mc_geom, cb_mc_hits });

    auto main_layout = Container::Vertical({ Container::Horizontal({ col1_layout, col2_layout }),
        Container::Horizontal({ btn_run, btn_exit }) | center });

    auto renderer = Renderer(main_layout,
        [&]
        {
            return vbox({ text(" CHeTFitter - Control Panel ") | bold | center | color(Color::Cyan)
                    | border,
                hbox({
                    // COLONNA SINISTRA
                    vbox({ window(text(" I/O & Geometry "),
                               vbox({
                                   hbox(text("InDir:  "), inp_input_dir->Render()),
                                   hbox(text("Input:  "), inp_input->Render()),
                                   hbox(text("OutDir: "), inp_output_dir->Render()),
                                   hbox(text("GeoDir: "), inp_geom_dir->Render()),
                                   hbox(text("Geom:   "), inp_geom->Render()),
                               })),
                        window(text(" B-Field "),
                            vbox({
                                hbox(text("Type: "), menu_btype->Render()),
                                hbox(text("B-Dir: "), inp_bdir->Render()),
                                hbox(text("Path:  "), inp_bpath->Render()),
                                hbox(text("Scale/Const: "),
                                    inp_bscale->Render() | size(WIDTH, EQUAL, 8), text(" / "),
                                    inp_bconst->Render() | size(WIDTH, EQUAL, 8)),
                            })),
                        window(text(" Finder/Fitter Settings "),
                            vbox({ hbox(text("Algorithm: "), menu_fitter->Render()),
                                hbox({ cb_prefitter->Render(), separator(),
                                    cb_quiet->Render() }) })) })
                        | flex,

                    // COLONNA DESTRA
                    vbox({ window(text(" Analysis Settings "),
                               vbox({ hbox(text("ToA Min/Max: "),
                                          inp_toa_min->Render() | size(WIDTH, EQUAL, 8),
                                          text(" / "),
                                          inp_toa_max->Render() | size(WIDTH, EQUAL, 8)),
                                   hbox(text("ToT Min/Max: "),
                                       inp_tot_min->Render() | size(WIDTH, EQUAL, 8), text(" / "),
                                       inp_tot_max->Render() | size(WIDTH, EQUAL, 8)),
                                   hbox(text("Turns:       "),
                                       inp_turns->Render() | size(WIDTH, EQUAL, 8)) })),
                        window(text(" Event Loop "),
                            vbox({ hbox(text("Run Mode:  "), menu_runmode->Render()),
                                hbox(text("Event ID:  "),
                                    inp_event->Render() | size(WIDTH, EQUAL, 10)),
                                hbox(text("Start/End: "),
                                    inp_start->Render() | size(WIDTH, EQUAL, 8), text(" / "),
                                    inp_entrymax->Render() | size(WIDTH, EQUAL, 8)) })),
                        window(text(" Assumed Recon Geometry "),
                            vbox({ cb_over_geom->Render(),
                                config.overrideGeom ? vbox({ hbox(text("  Active:  "),
                                                                 cb_cyls->Render()),
                                    hbox(text("  Trans:   "),
                                        inp_tx->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_ty->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_tz->Render() | size(WIDTH, EQUAL, 6)),
                                    hbox(text("  Rot:     "),
                                        inp_rx->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_ry->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_rz->Render() | size(WIDTH, EQUAL, 6)),
                                    hbox(text("  Off/D0-2:"),
                                        inp_off->Render() | size(WIDTH, EQUAL, 6), text(" |"),
                                        inp_d[0]->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_d[1]->Render() | size(WIDTH, EQUAL, 6), text(","),
                                        inp_d[2]->Render() | size(WIDTH, EQUAL, 6)),
                                    hbox(text("  D3-5:    "), filler() | size(WIDTH, EQUAL, 8),
                                        text(" |"), inp_d[3]->Render() | size(WIDTH, EQUAL, 6),
                                        text(","), inp_d[4]->Render() | size(WIDTH, EQUAL, 6),
                                        text(","), inp_d[5]->Render() | size(WIDTH, EQUAL, 6)) })
                                                    : text("  (Using Nominal Geometry) ") | dim })),
                        window(text(" MC Only Options "),
                            vbox({ hbox({ cb_pttrec->Render(), cb_smearing->Render() }),
                                hbox({ cb_mc_geom->Render(), cb_mc_hits->Render() }) })) })
                        | flex,
                }),
                filler(), hbox({ btn_run->Render(), text("   "), btn_exit->Render() }) | center });
        });

    // --- 5. ESECUZIONE ---
    screen.Loop(renderer);

    std::cout << "\x1b[?25h\x1b[0 q\x1b[0m\x1b[2J\x1b[H" << std::flush;
    if(should_exit)
        exit(0);

    // --- 6. RIVERSAMENTO DATI ---
    config.inputDir = input_dir;
    config.inputFileName = input_file;
    if(!input_file.empty())
    {
        config.inputDataFiles.clear();
        config.inputDataFiles.push_back(input_dir + input_file);
    }
    config.outputDir = output_dir;
    config.outputFile = output_dir + "fitted_" + input_file;
    config.geomDir = geom_dir;
    config.geomFileName = geom_file;
    config.geometryFile = geom_dir + geom_file;
    config.bFieldDir = b_dir;
    config.bFieldMapName = b_path;
    config.bFieldPath = b_dir + b_path;

    auto safe_stod = [](const std::string &s, double def)
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
    auto safe_stol = [](const std::string &s, long def)
    {
        try
        {
            return std::stol(s);
        }
        catch(...)
        {
            return def;
        }
    };
    auto safe_stoi = [](const std::string &s, int def)
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

    config.bScale = safe_stod(b_scale_str, config.bScale);
    config.bConstZ = safe_stod(b_const_str, config.bConstZ);
    config.rangeLoop.first = safe_stol(start_str, 0);
    config.rangeLoop.second = safe_stol(entry_max_str, 1000);
    config.event = safe_stoi(event_str, 0);
    config.turnID = safe_stoi(turns_str, 0);
    config.cutToAMin = safe_stod(toa_min_str, 215.0);
    config.cutToAMax = safe_stod(toa_max_str, 250.0);
    config.cutToTMin = safe_stod(tot_min_str, 40.0);
    config.cutToTMax = safe_stod(tot_max_str, 220.0);

    config.geom_tx = safe_stod(tx_str, 0.0);
    config.geom_ty = safe_stod(ty_str, 0.0);
    config.geom_tz = safe_stod(tz_str, 0.0);
    config.geom_rx = safe_stod(rx_str, 0.0);
    config.geom_ry = safe_stod(ry_str, 0.0);
    config.geom_rz = safe_stod(rz_str, 0.0);
    config.geom_offset = safe_stod(off_str, 0.0);

    config.geom_deltas.resize(6);
    for(int i = 0; i < 6; ++i)
        config.geom_deltas[i] = safe_stod(d_str[i], 0.0);

    if(config.turnID > 0)
        config.turnMode = true;
    if(selected_runmode == 0)
    {
        config.processSingle = true;
        config.rangeLoop = { 0, 99999999 };
    }
    else if(selected_runmode == 1)
    {
        config.processSingle = true;
        config.event = -1;
        config.rangeLoop = { 0, 99999999 };
    }
    else
    {
        config.processSingle = false;
    }

    config.fitter = fitter_values[selected_fitter];
    config.bFieldType = bfield_values[selected_btype];
}

Int_t main(Int_t argc, char **argv)
{
    if(!gApplication)
    {
        int rootArgc = 0;
        char **rootArgv = nullptr;
        new TApplication("CHeTFitter", &rootArgc, rootArgv);

        // Load rootlogon if present
        if(!gSystem->AccessPathName(".rootlogon.C"))
        {
            gROOT->ProcessLine(".x .rootlogon.C");
        }
    }

    auto &config = Config::get();
    CLI::App app { "CHeTFitter Analysis Tool" };
    app.set_config("--config", "macros/edm_tracker.mac", "Read a configuration file", true);

    // --- MAPS FOR ENUM ---
    // Fitter
    map<string, FitterType> mapFitter { { "cosmic", FitterType::CosmicFitter },
        { "spacepoint", FitterType::SpacepointFitter }, { "helix", FitterType::HelixFitter } };

    // B-Field Type
    map<string, muEDM::Fields::FieldBType> mapBField { { "multi",
                                                           muEDM::Fields::FieldBType::MultiMap },
        { "single", muEDM::Fields::FieldBType::SingleMap },
        { "constant", muEDM::Fields::FieldBType::Constant } };

    // --- DEFINE OPTIONS ---

    // I/O Data Files
    app.add_option("--inDir", config.inputDir, "Input directory");
    app.add_option("-i,--input", config.inputFileName, "Input ROOT data file name");
    app.add_option("--outDir", config.outputDir, "Output directory");
    // outputFile will be set in POST-PARSING LOGIC

    // Event loop
    auto optEvent = app.add_option("-e,--event", config.event, "Single event mode");
    app.add_option("-M,--entryMax", config.rangeLoop.second, "End event loop at this entry");
    app.add_option(
        "--rangeLoop", config.rangeLoop, "Start and stop entries for event loop [start, stop)");

    // Geometry / Bfields Files
    app.add_option("--geomDir", config.geomDir, "Detector geometry directory");
    app.add_option(
        "-G,--geometry", config.geomFileName, "Detector geometry file name (.gdml/.root)");
    app.add_option("--b-type", config.bFieldType, "Magnetic Field Type (multi, single, constant)")
        ->transform(CLI::CheckedTransformer(mapBField, CLI::ignore_case));
    app.add_option("--b-dir", config.bFieldDir, "Directory containing field maps");
    app.add_option("--b-path", config.bFieldMapName, "Name of the field map folder/file");
    app.add_option("--b-scale", config.bScale, "Global B scaling factor");
    app.add_option("--b-const", config.bConstZ, "Bz value for constant field mode (kGauss)");

    // -F [type]
    app.add_option("-F,--fitter", config.fitter, "Fitter type: cosmic, spacepoint, helix")
        ->transform(CLI::CheckedTransformer(mapFitter, CLI::ignore_case));
    // -P
    app.add_flag("-P,--prefitter", config.usePrefitter, "Use pre-fitter");
    // -t [nTurns]
    auto optTurn = app.add_option("-t,--turns", config.turnID, "Apply turn analysis (set Turn ID)");
    // -s
    app.add_flag("-s,--smearing", config.useSmearing, "Apply smearing analysis");
    // -u
    app.add_flag("-u,--use-true-hits", config.useTrueMCHits, "Use true MC hits");
    // -p
    app.add_flag("-p,--pattern-rec", config.pttrecMode, "Apply pattern recognition analysis");
    // -q
    app.add_flag("-q,--quiet", config.quietMode, "Quiet mode");

    app.add_option("--toa-min", config.cutToAMin, "Minimum ToA cut");
    app.add_option("--toa-max", config.cutToAMax, "Maximum ToA cut");
    app.add_option("--tot-min", config.cutToTMin, "Minimum ToT cut");
    app.add_option("--tot-max", config.cutToTMax, "Maximum ToT cut");

    // --- PARSING / INTERFACE ---
    try
    {
        app.parse(argc, argv);
    }
    catch(const CLI::ParseError &e)
    {
        return app.exit(e);
    }

    // If -e, enable processSingle
    if(optEvent->count() > 0)
    {
        config.processSingle = true;
        config.rangeLoop = { 0, numeric_limits<Long_t>::max() };
    }

    // If -t, enable turnMode
    if(optTurn->count() > 0)
        config.turnMode = true;

    if(argc == 1)
    {
        // Se non passo argomenti, entro in modalità interattiva
        RunInteractiveGUI(config);
    }
    else
    {
        config.inputDataFiles.clear();
        config.inputDataFiles.push_back(config.inputDir + config.inputFileName);
        config.outputFile = config.outputDir + "fitted_" + config.inputFileName;
        config.geometryFile = config.geomDir + config.geomFileName;
        config.bFieldPath = config.bFieldDir + config.bFieldMapName;
    }

    // --- POST-PARSING LOGIC ---

    if(config.quietMode)
    {
        gROOT->SetBatch(kTRUE);
    }
    else
    {
        if(!gEve)
            TEveManager::Create(kTRUE, "FVI");
    }

    // --- EXECUTION ---

    // Select the fitter
    switch(config.fitter)
    {
        case FitterType::CosmicFitter:
            cout << ">>> Running Cosmic Fitter...\n" << endl;
            FITALG::CosmicFitter();
            break;
        case FitterType::SpacepointFitter:
            cout << ">>> Running Spacepoint Fitter...\n" << endl;
            FITALG::SpacepointFitter();
            break;
        case FitterType::HelixFitter:
            cout << ">>> Running Helix Fitter...\n" << endl;
            FITALG::HelixFitter();
            break;
        default:
            cerr << ">>> Unknown fitter type or not specified!" << endl;
            cout << app.help() << endl;
            // CLI11 può anche rendere obbligatorio -F se vuoi:
            // app.add_option(...)->required();
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
