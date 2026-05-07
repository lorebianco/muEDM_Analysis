#ifndef LBROOTSTYLE_HH
#define LBROOTSTYLE_HH

#include <TColor.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

// ============================================================================
// Lorenzo's Publication Style for Experimental Particle Physics
// ============================================================================

inline void SetLBStyle()
{

    TStyle *lbStyle = new TStyle("lbStyle", "Lorenzo's Publication Style");

    // ---------------------------------------------------------
    // 1. Colori e Sfondi (Tutto bianco, niente bordi grigi)
    // ---------------------------------------------------------
    lbStyle->SetCanvasColor(kWhite);
    lbStyle->SetPadColor(kWhite);
    lbStyle->SetFrameFillColor(kWhite);
    lbStyle->SetStatColor(kWhite);
    lbStyle->SetTitleFillColor(kWhite);

    lbStyle->SetCanvasBorderMode(0);
    lbStyle->SetPadBorderMode(0);
    lbStyle->SetFrameBorderMode(0);
    lbStyle->SetTitleBorderSize(0);

    // ---------------------------------------------------------
    // 2. Dimensioni Canvas e Margini
    // ---------------------------------------------------------
    // Margini ampi per gli assi, ridotti in alto e a destra.
    lbStyle->SetCanvasDefW(800);
    lbStyle->SetCanvasDefH(600);

    lbStyle->SetPadTopMargin(0.06); // Poco spazio sopra (niente titolo del plot)
    lbStyle->SetPadBottomMargin(0.14); // Spazio per il titolo asse X
    lbStyle->SetPadLeftMargin(0.16); // Spazio per il titolo asse Y (evita sovrapposizioni)
    lbStyle->SetPadRightMargin(0.05); // Poco spazio a destra
    // NOTA: se fai un plot 2D (COLZ), cambia manualmente il margine destro nel
    // tuo script: gPad->SetRightMargin(0.15);

    // ---------------------------------------------------------
    // 3. Font (La regola d'oro: Font 42)
    // ---------------------------------------------------------
    Int_t fontId = 42; // Helvetica (Standard CERN/PRL/Nature)
    Double_t tsize = 0.05; // Dimensione testo standard

    lbStyle->SetTextFont(fontId);
    lbStyle->SetTextSize(tsize);

    lbStyle->SetLabelFont(fontId, "xyz");
    lbStyle->SetLabelSize(tsize, "xyz");
    lbStyle->SetTitleFont(fontId, "xyz");
    lbStyle->SetTitleSize(tsize, "xyz");

    // ---------------------------------------------------------
    // 4. Impostazioni degli Assi
    // ---------------------------------------------------------
    lbStyle->SetTitleOffset(1.2, "X");
    lbStyle->SetTitleOffset(1.4, "Y"); // Allontana la Y per non coprire i numeri
    lbStyle->SetTitleOffset(1.2, "Z");

    lbStyle->SetLabelOffset(0.01, "xyz");

    // Ottimizza i ticks degli assi (fondamentale per le pubblicazioni)
    lbStyle->SetNdivisions(505, "xyz"); // Evita sovraffollamento dei numeri
    lbStyle->SetStripDecimals(kFALSE); // Mantiene l'allineamento (es. 1.0, 1.5, 2.0)

    // Disegna le tacche anche sui lati superiore e destro del frame (Standard
    // HEP)
    lbStyle->SetPadTickX(1);
    lbStyle->SetPadTickY(1);

    // ---------------------------------------------------------
    // 5. Istogrammi, Grafici (TGraph) e Marker
    // ---------------------------------------------------------
    lbStyle->SetHistLineWidth(2);
    lbStyle->SetHistLineColor(kBlack);

    // TGraph (Punti sperimentali)
    lbStyle->SetMarkerStyle(20); // Pallino pieno
    lbStyle->SetMarkerSize(0.9);
    lbStyle->SetMarkerColor(kBlack);
    lbStyle->SetEndErrorSize(0); // Niente sbarrette orizzontali agli estremi degli errori
    lbStyle->SetErrorX(0.5); // Dimensione standard per le barre su X

    // Forza l'asse Y a partire da zero quando si riempie l'istogramma
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 16, 0)
    lbStyle->SetHistMinimumZero(kTRUE);
#endif

    // ---------------------------------------------------------
    // 6. Legenda e Statistiche
    // ---------------------------------------------------------
    // Nelle pubblicazioni, Stat e Titoli dei plot sono disattivati (si usano le
    // caption testuali)
    lbStyle->SetOptStat(0);
    lbStyle->SetOptFit(0);
    lbStyle->SetOptTitle(0);

    // Se proprio devi usare la stat box per analisi interne, rendila pulita:
    lbStyle->SetStatFont(fontId);
    lbStyle->SetStatFontSize(0.04);
    lbStyle->SetStatBorderSize(1);
    lbStyle->SetStatX(0.92);
    lbStyle->SetStatY(0.92);

    // Legenda trasparente e senza bordo
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 16, 0)
    lbStyle->SetLegendBorderSize(0);
    // lbStyle->SetLegendFillStyle(0);
    lbStyle->SetLegendFont(fontId);
#endif

    // ---------------------------------------------------------
    // 7. Funzioni (Fit) e Palette 2D
    // ---------------------------------------------------------
    lbStyle->SetFuncColor(kRed + 1); // Un rosso leggermente più scuro e visibile
    lbStyle->SetFuncWidth(2);

    lbStyle->SetPalette(kBird); // Palette perceptually uniform (ottima scelta)
    lbStyle->SetNumberContours(255);

    // ---------------------------------------------------------
    // 8. Applicazione dello stile
    // ---------------------------------------------------------
    // lbStyle->SetCanvasPreferGL(kTRUE); // (Sconsigliato per l'export PDF dei
    // paper, lascialo off)

    gROOT->SetStyle("lbStyle");
    gROOT->ForceStyle();

    // Forza ROOT a usare l'esponenziale (10^x) standard di TGaxis se i numeri
    // sono grandi
    TGaxis::SetMaxDigits(4);

    // Aggiorna la GUI
    gSystem->ProcessEvents();

    std::cout << " [lbStyle] Lorenzo's Publication Style successfully loaded!" << std::endl;
}

// --- Helper functions ---

inline void AddWatermark(const char *text = "PRELIMINARY")
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

inline TPaveText *AddAcquisitionDate(const char *start, const char *stop, Bool_t fromto)
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

#endif // LBROOTSTYLE_HH
