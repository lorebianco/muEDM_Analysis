// rootlogon.C
// Lorenzo Bianco
//

{
    // ==============================================
    // 0. CHeT Library Auto-Load
    // ==============================================
    // Definiamo i percorsi relativi alla cartella corrente
    TString libPath = "CHeTLibrary/install/lib/libCHeTLib.so";
    TString incPath = "CHeTLibrary/install/include";

    // gSystem->AccessPathName restituisce FALSE (0) se il file ESISTE (sì, è al
    // contrario)
    if(!gSystem->AccessPathName(libPath))
    {

        // 1. Carica la libreria condivisa
        gSystem->Load(libPath);

        // 2. Aggiunge il path degli header per l'interprete (Cling)
        //    Così puoi fare #include "CHeT/GlobalSettings.hh" nelle macro
        gInterpreter->AddIncludePath(incPath);
        gInterpreter->ProcessLine("#include \"CHeT/CHeTVisualizer.hh\"");

        cout << " [CHeTLib] Library loaded successfully!" << endl;
        cout << "           Path: " << libPath << endl;
    }
    else
    {
        cout << " [CHeTLib] WARNING: Library NOT found!" << endl;
        cout << "           Expected at: " << libPath << endl;
        cout << "           Did you run 'make install' in CHeTLibrary/build?" << endl;
    }

    // ==============================================
    // Lorenzo's Style Definitions
    // ==============================================

    TStyle *lbStyle = new TStyle("lbStyle", "Lorenzo's Root Styles");

    // 1. OpenGL Force
    lbStyle->SetCanvasPreferGL(kTRUE);

    // 2. Colori di Sfondo
    lbStyle->SetCanvasColor(kWhite);
    lbStyle->SetPadColor(kWhite);
    lbStyle->SetFrameFillColor(kWhite);
    lbStyle->SetFrameBorderMode(0);
    lbStyle->SetPadBorderMode(0);
    lbStyle->SetCanvasBorderMode(0);

    // 3. Fonts & Text
    Int_t fontId = 22;
    Double_t titleSize = 0.05;
    Double_t labelSize = 0.045;

    lbStyle->SetTextFont(fontId);
    lbStyle->SetTextSize(labelSize);
    lbStyle->SetLabelFont(fontId, "xyz");
    lbStyle->SetLabelSize(labelSize, "xyz");
    lbStyle->SetLabelOffset(0.01, "xyz");
    lbStyle->SetTitleFont(fontId, "xyz");
    lbStyle->SetTitleSize(titleSize, "xyz");
    lbStyle->SetTitleOffset(1.1, "x");
    lbStyle->SetTitleOffset(1.2, "y");
    lbStyle->SetTitleFont(fontId, "");
    lbStyle->SetTitleBorderSize(0);
    lbStyle->SetTitleStyle(0);
    lbStyle->SetTitleX(0.5);
    lbStyle->SetTitleAlign(23);
    lbStyle->SetTitleH(0.050);

    // Data/Ora
    lbStyle->SetOptDate(0);
    lbStyle->GetAttDate()->SetTextFont(fontId);
    lbStyle->GetAttDate()->SetTextSize(0.03);

    // 4. Color Palette
    lbStyle->SetPalette(kBird);
    lbStyle->SetNumberContours(255);

    // 5. Statistics Box
    lbStyle->SetOptStat(1110);
    lbStyle->SetOptFit(111);
    lbStyle->SetStatX(0.90);
    lbStyle->SetStatY(0.90);
    lbStyle->SetStatW(0.20);
    lbStyle->SetStatH(0.15);
    lbStyle->SetStatColor(kWhite);
    lbStyle->SetStatStyle(1001);
    lbStyle->SetStatBorderSize(1);
    lbStyle->SetStatFont(fontId);
    lbStyle->SetStatFontSize(0.03);
    lbStyle->SetFitFormat("5.4g");

    // 6. Canvas & Margins
    lbStyle->SetCanvasDefW(800);
    lbStyle->SetCanvasDefH(600);
    lbStyle->SetPadBottomMargin(0.13);
    lbStyle->SetPadTopMargin(0.08);
    lbStyle->SetPadLeftMargin(0.14);
    lbStyle->SetPadRightMargin(0.12);
    lbStyle->SetFrameLineWidth(1);

    // 7. Histogram Styling
    lbStyle->SetHistLineWidth(2);
    lbStyle->SetHistLineColor(kBlue + 1);
    lbStyle->SetFuncColor(kRed);
    lbStyle->SetFuncWidth(2);

    // 8. Apply Style
    gROOT->SetStyle("lbStyle");
    gROOT->ForceStyle();

    cout << " [Style]   Lorenzo's Style Loaded!" << endl;

    return;
}

// Funzioni Helper definite fuori dallo scope anonimo per essere visibili
// globalmente

void AddWatermark(const char *text = "PRELIMINARY")
{
    if(!gPad)
    {
        printf("Error: No Canvas/Pad\n");
        return;
    }
    TLatex *watermark = new TLatex();
    watermark->SetNDC();
    watermark->SetTextFont(42);
    watermark->SetTextSize(0.12);
    watermark->SetTextAlign(22);
    watermark->SetTextAngle(45);
    watermark->SetTextColorAlpha(kGray + 1, 0.3);
    watermark->DrawLatex(0.5, 0.5, text);
    gPad->Update();
}

TPaveText *print_acquisition_date(const char *start, const char *stop, Bool_t fromto)
{
    TString str_0 = "Acquisition of ";
    TString str_1 = "Acquisition from ";
    TString str_2 = " to ";
    TString str = str_1 + start + str_2 + stop;
    TString str_s = str_0 + start;
    TPaveText *pt = new TPaveText(0.02, 0.95, 0.32, 1, "NB NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.035);
    if(fromto)
        pt->AddText(str);
    else
        pt->AddText(str_s);
    return pt;
}
