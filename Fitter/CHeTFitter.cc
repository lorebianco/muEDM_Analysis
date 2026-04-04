#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <iostream>
#include <map>
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

void RunInteractiveGUI(Config &config)
{
    using namespace ftxui;

    // --- 1. MAPPATURA ENUM (Per risolvere il problema del Fitter che non cambia)
    // ---
    std::vector<FitterType> fitter_values
        = { FitterType::CosmicFitter, FitterType::SpacepointFitter, FitterType::HelixFitter };

    std::vector<muEDM::Fields::FieldBType> bfield_values = { muEDM::Fields::FieldBType::MultiMap,
        muEDM::Fields::FieldBType::SingleMap, muEDM::Fields::FieldBType::Constant };

    // --- 2. DICHIARAZIONE VARIABILI (Prima dei componenti!) ---
    auto screen = ScreenInteractive::TerminalOutput();

    // Stringhe di appoggio per gli input di testo
    std::string input_dir = config.inputDir;
    std::string input_file = config.inputFileName;
    std::string output_dir = config.outputDir;
    std::string geom_dir = config.geomDir;
    std::string geom_file = config.geomFileName;
    std::string b_dir = config.bFieldDir;
    std::string b_path = config.bFieldMapName;
    std::string b_scale_str = std::to_string(config.bScale);
    std::string b_const_str = std::to_string(config.bConstZ);
    std::string entry_max_str = std::to_string(config.rangeLoop.second);
    std::string turns_str = std::to_string(config.turnID);

    std::string toa_min_str = std::to_string(config.cutToAMin);
    std::string toa_max_str = std::to_string(config.cutToAMax);
    std::string tot_min_str = std::to_string(config.cutToTMin);
    std::string tot_max_str = std::to_string(config.cutToTMax);

    // Run Mode
    int selected_runmode = 0;
    if(config.processSingle)
    {
        if(config.event == -1)
            selected_runmode = 1; // Random Single
        else
            selected_runmode = 0; // Specific Single
    }
    else
    {
        selected_runmode = 2; // Loop
    }
    std::string event_str = std::to_string(config.event);
    std::string start_str = std::to_string(config.rangeLoop.first);

    // Indici per i menu (Radiobox o Dropdown)
    int selected_fitter = 0;
    for(size_t i = 0; i < fitter_values.size(); ++i)
    {
        if(config.fitter == fitter_values[i])
            selected_fitter = (int)i;
    }

    int selected_btype = 0;
    for(size_t i = 0; i < bfield_values.size(); ++i)
    {
        if(config.bFieldType == bfield_values[i])
            selected_btype = (int)i;
    }

    // --- 3. CREAZIONE COMPONENTI ---
    std::vector<std::string> fitter_entries = { "Cosmic", "Spacepoint", "Helix" };
    std::vector<std::string> bfield_entries = { "MultiMap", "SingleMap", "Constant" };

    Component inp_input_dir = Input(&input_dir, "path/to/input/");
    Component inp_input = Input(&input_file, "input.root");
    Component inp_output_dir = Input(&output_dir, "path/to/output/");
    Component inp_geom_dir = Input(&geom_dir, "path/to/geom/");
    Component inp_geom = Input(&geom_file, "geom.gdml");

    Component menu_btype = Radiobox(&bfield_entries, &selected_btype);
    Component inp_bdir = Input(&b_dir, "path/to/maps/");
    Component inp_bpath = Input(&b_path, "map_folder");
    Component inp_bscale = Input(&b_scale_str, "-1.0");
    Component inp_bconst = Input(&b_const_str, "28.9");

    Component menu_fitter = Radiobox(&fitter_entries, &selected_fitter);
    Component cb_prefitter = Checkbox("Use Pre-fitter", &config.usePrefitter);
    Component cb_quiet = Checkbox("Quiet Mode", &config.quietMode);
    Component cb_smearing = Checkbox("Apply Smearing", &config.useSmearing);
    Component cb_pttrec = Checkbox("Pattern Recognition", &config.pttrecMode);

    std::vector<std::string> runmode_entries = { "Single Event", "Random Event", "Event Loop" };
    Component menu_runmode = Radiobox(&runmode_entries, &selected_runmode);
    Component inp_event = Input(&event_str, "Event ID");
    Component inp_start = Input(&start_str, "Start Entry");
    Component inp_entrymax = Input(&entry_max_str, "1000");
    Component inp_turns = Input(&turns_str, "0");

    Component inp_toa_min = Input(&toa_min_str, "215.0");
    Component inp_toa_max = Input(&toa_max_str, "250.0");
    Component inp_tot_min = Input(&tot_min_str, "40.0");
    Component inp_tot_max = Input(&tot_max_str, "220.0");

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

    // --- 4. LAYOUT E RENDERER ---
    auto layout = Container::Vertical(
        { inp_input_dir, inp_input, inp_output_dir, inp_geom_dir, inp_geom, menu_btype, inp_bdir,
            inp_bpath, inp_bscale, inp_bconst, menu_fitter, cb_prefitter, cb_quiet, menu_runmode,
            inp_event, inp_start, inp_entrymax, inp_turns, inp_toa_min, inp_toa_max, inp_tot_min,
            inp_tot_max, cb_smearing, cb_pttrec, Container::Horizontal({ btn_run, btn_exit }) });

    auto renderer = Renderer(layout,
        [&]
        {
            return vbox({ text(" CHeTFitter - Control Panel ") | bold | center | color(Color::Cyan)
                    | border,
                hbox({
                    vbox({
                        window(text(" I/O & Geometry "),
                            vbox({
                                hbox(text("InDir:  "), inp_input_dir->Render()),
                                hbox(text("Input:  "), inp_input->Render()),
                                hbox(text("OutDir: "), inp_output_dir->Render()),
                                hbox(text("GeoDir: "), inp_geom_dir->Render()),
                                hbox(text("Geom:   "), inp_geom->Render()),
                            })),
                        window(text(" Fitter "),
                            vbox({
                                menu_fitter->Render(),
                                hbox({ cb_prefitter->Render(), separator(), cb_quiet->Render() }),
                            })),
                    }) | flex,
                    vbox({
                        window(text(" B-Field "),
                            vbox({
                                hbox(text("Type: "), menu_btype->Render()),
                                hbox(text("B-Dir: "), inp_bdir->Render()),
                                hbox(text("Path: "), inp_bpath->Render()),
                                hbox(text("Scale/Const: "), inp_bscale->Render(), text(" / "),
                                    inp_bconst->Render()),
                            })),
                        window(text(" Events & Sim "),
                            vbox({
                                hbox(text("Run Mode: "), menu_runmode->Render()),
                                hbox(text("Event ID: "), inp_event->Render()),
                                hbox(text("Start/End Entries: "), inp_start->Render(), text(" / "),
                                    inp_entrymax->Render()),
                                hbox(text("Turns:   "), inp_turns->Render()),
                                hbox(text("ToA Min/Max: "), inp_toa_min->Render(), text(" / "),
                                    inp_toa_max->Render()),
                                hbox(text("ToT Min/Max: "), inp_tot_min->Render(), text(" / "),
                                    inp_tot_max->Render()),
                                cb_smearing->Render(),
                                cb_pttrec->Render(),
                            })),
                    }) | flex,
                }),
                filler(),
                hbox({ btn_run->Render() | flex | color(Color::Green),
                    btn_exit->Render() | flex | color(Color::Red) }) });
        });

    // --- 5. ESECUZIONE ---
    screen.Loop(renderer);

    std::cout << "\x1b[?25h\x1b[0 q\x1b[0m\x1b[2J\x1b[H" << std::flush;

    if(should_exit)
        exit(0);

    // --- 6. RIVERSAMENTO DATI (Dalla GUI al config) ---
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

    // Conversione stringhe -> numeri
    try
    {
        config.bScale = std::stod(b_scale_str);
    }
    catch(...)
    {
    }
    try
    {
        config.bConstZ = std::stod(b_const_str);
    }
    catch(...)
    {
    }
    try
    {
        config.rangeLoop.first = std::stol(start_str);
    }
    catch(...)
    {
    }
    try
    {
        config.rangeLoop.second = std::stol(entry_max_str);
    }
    catch(...)
    {
    }
    try
    {
        config.event = std::stoi(event_str);
    }
    catch(...)
    {
    }
    try
    {
        config.turnID = std::stoi(turns_str);
    }
    catch(...)
    {
    }

    try
    {
        config.cutToAMin = std::stod(toa_min_str);
    }
    catch(...)
    {
    }
    try
    {
        config.cutToAMax = std::stod(toa_max_str);
    }
    catch(...)
    {
    }
    try
    {
        config.cutToTMin = std::stod(tot_min_str);
    }
    catch(...)
    {
    }
    try
    {
        config.cutToTMax = std::stod(tot_max_str);
    }
    catch(...)
    {
    }

    if(config.turnID > 0)
        config.turnMode = true;

    if(selected_runmode == 0)
    {
        config.processSingle = true;
        config.rangeLoop = { 0, std::numeric_limits<Long_t>::max() };
    }
    else if(selected_runmode == 1)
    {
        config.processSingle = true;
        config.event = -1;
        config.rangeLoop = { 0, std::numeric_limits<Long_t>::max() };
    }
    else
    {
        config.processSingle = false;
    }

    // Mappatura Enum corretta
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
