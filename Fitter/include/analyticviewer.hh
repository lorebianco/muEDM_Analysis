#ifndef ANALYTICVIEWER_HH
#define ANALYTICVIEWER_HH

#include <TEveBrowser.h>
#include <TEveEventManager.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveScene.h>
#include <TEveStraightLineSet.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGWidget.h>
#include <TGeoManager.h>

#include "CHeT/CHeTGlobalSettings.hh"
#include "CHeT/CHeTVisualizer.hh"
#include "fitteralgorithms.hh"

namespace FITALG
{

class AnalyticViewer
{
  private:
    struct EventData
    {
        int eventID;
        FITALG::FitOutput result;
    };

    std::vector<EventData> events;
    int current_idx = 0;

    // GUI
    TGNumberEntry *event_num_display = nullptr;
    TGCheckButton *cb_geom = nullptr;
    TGCheckButton *cb_hits = nullptr;
    TGCheckButton *cb_track = nullptr;

    bool drawGeom = true;
    bool drawHits = true;
    bool drawTrack = true;

    class GUIFrame : public TGMainFrame
    {
      public:
        AnalyticViewer *viewer;
        GUIFrame(AnalyticViewer *v)
            : TGMainFrame(gClient->GetRoot(), 300, 400)
            , viewer(v)
        {
        }
        Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override
        {
            return viewer->ProcessMessage(msg, parm1, parm2);
        }
    };

    GUIFrame *frmMain = nullptr;

  public:
    AnalyticViewer()
        : current_idx(0)
    {
        // Inizializza l'ambiente EVE se non esiste
        if(!gEve)
            TEveManager::Create();
    }

    // Carica la geometria ed imposta lo stile
    void LoadGeometry(const char *filepath)
    {
        if(!gGeoManager)
            TGeoManager::Import(filepath);

        if(gGeoManager && gGeoManager->GetTopNode())
        {
            TObjArray *volumes = gGeoManager->GetListOfVolumes();
            if(volumes)
            {
                for(int v = 0; v < volumes->GetEntries(); ++v)
                {
                    TGeoVolume *vol = (TGeoVolume *)volumes->At(v);
                    if(vol)
                    {
                        vol->SetLineColor(12);
                        vol->SetTransparency(50);
                    }
                }
            }
            TEveGeoTopNode *eveNode = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
            eveNode->IncDenyDestroy();
            eveNode->SetMainTransparency(50);
            gEve->AddGlobalElement(eveNode);
        }
    }

    // Carica i risultati del tuo fit analitico
    void AddEvent(int id, const FITALG::FitOutput &res)
    {
        events.push_back({ id, res });
    }

    void Show()
    {
        if(events.empty())
        {
            printf(">>> Viewer: Nessun evento da mostrare.\n");
            return;
        }
        BuildGUI();
        DrawEvent(0);
    }

    void DrawEvent(int idx)
    {
        if(idx < 0 || idx >= (int)events.size())
            return;
        current_idx = idx;

        // Pulisce la visualizzazione precedente
        if(!gEve->GetCurrentEvent())
            gEve->AddEvent(new TEveEventManager("Event", "Analytic Event"));
        else
            gEve->GetCurrentEvent()->DestroyElements();

        auto &ev = events[idx];
        auto &tr = ev.result.track;

        // 1. Disegna i punti (Hits fittate sulle fibre) usando TEveStraightLineSet per il depth
        // bypass
        if(drawHits)
        {
            TEveStraightLineSet *ps
                = new TEveStraightLineSet(Form("Fitted Hits (Ev %d)", ev.eventID));
            ps->SetMarkerColor(kRed);
            ps->SetMarkerStyle(20);
            ps->SetMarkerSize(1.5);
            ps->SetRnrMarkers(kTRUE);
            ps->SetRnrLines(kFALSE);
            ps->SetDepthTest(kFALSE); // Il trucco magico per i marker on top

            for(const auto &pt : ev.result.fittedPoints)
            {
                double x = pt.x, y = pt.y, z = pt.z;
                CHeT::Config::ApplyTransformation(x, y, z);
                ps->AddMarker(x / 10.0, y / 10.0, z / 10.0);
            }
            gEve->AddElement(ps);
        }

        // 2. Disegna la retta del FIT analitico
        if(drawTrack && tr.converged)
        {
            TEveStraightLineSet *ls = new TEveStraightLineSet("Analytical Track");
            ls->SetLineColor(kCyan);
            ls->SetLineWidth(3);
            ls->SetDepthTest(kFALSE); // Linea sempre visibile

            double y_start = -500.0, y_end = 500.0;
            double x_start = tr.x0 + tr.sx * y_start, z_start = tr.z0 + tr.sz * y_start;
            double x_end = tr.x0 + tr.sx * y_end, z_end = tr.z0 + tr.sz * y_end;

            CHeT::Config::ApplyTransformation(x_start, y_start, z_start);
            CHeT::Config::ApplyTransformation(x_end, y_end, z_end);

            ls->AddLine(x_start / 10.0, y_start / 10.0, z_start / 10.0, x_end / 10.0, y_end / 10.0,
                z_end / 10.0);
            gEve->AddElement(ls);
        }

        // Handle geometry visibility
        if(gEve->GetGlobalScene())
        {
            for(auto el : gEve->GetGlobalScene()->RefChildren())
            {
                el->SetRnrState(drawGeom);
            }
        }

        // Aggiorna il numero dell'evento nella GUI
        if(event_num_display)
            event_num_display->SetIntNumber(ev.eventID);

        gEve->Redraw3D(kFALSE);
    }

    void BuildGUI()
    {
        TEveBrowser *browser = gEve->GetBrowser();
        browser->StartEmbedding(TRootBrowser::kLeft);

        frmMain = new GUIFrame(this);
        frmMain->SetWindowName("Analytic Viewer");
        frmMain->SetCleanup(kDeepCleanup); // Evita doppie chiusure che causano BadWindow

        TGVerticalFrame *vf = new TGVerticalFrame(frmMain, 300, 400);

        // Controllo Eventi
        TGHorizontalFrame *hf = new TGHorizontalFrame(vf);
        TGTextButton *bPrev = new TGTextButton(hf, " < Prev ", 1);
        bPrev->Associate(frmMain);
        event_num_display = new TGNumberEntry(hf, 0, 6, -1, TGNumberFormat::kNESInteger);
        TGTextButton *bNext = new TGTextButton(hf, " Next > ", 2);
        bNext->Associate(frmMain);
        hf->AddFrame(bPrev, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 5, 5, 5));
        hf->AddFrame(
            event_num_display, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 5, 5));
        hf->AddFrame(bNext, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 5, 5, 5, 5));
        vf->AddFrame(hf, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

        // Checkbox opzioni di visualizzazione
        TGGroupFrame *gf = new TGGroupFrame(vf, "Display Options");
        cb_geom = new TGCheckButton(gf, "Draw Geometry", 10);
        cb_geom->SetState(drawGeom ? kButtonDown : kButtonUp);
        cb_geom->Associate(frmMain);
        gf->AddFrame(cb_geom, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

        cb_hits = new TGCheckButton(gf, "Draw Hits", 11);
        cb_hits->SetState(drawHits ? kButtonDown : kButtonUp);
        cb_hits->Associate(frmMain);
        gf->AddFrame(cb_hits, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

        cb_track = new TGCheckButton(gf, "Draw Track", 12);
        cb_track->SetState(drawTrack ? kButtonDown : kButtonUp);
        cb_track->Associate(frmMain);
        gf->AddFrame(cb_track, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

        vf->AddFrame(gf, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

        frmMain->AddFrame(vf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        frmMain->MapSubwindows();
        frmMain->Resize();
        frmMain->MapWindow();

        browser->StopEmbedding();
        browser->SetTabTitle("Analytic Viewer", 0);
    }

    Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
    {
        (void)parm2; // Suppress unused warning
        switch(GET_MSG(msg))
        {
            case kC_COMMAND:
                switch(GET_SUBMSG(msg))
                {
                    case kCM_BUTTON:
                        if(parm1 == 1)
                            Prev();
                        else if(parm1 == 2)
                            Next();
                        break;
                    case kCM_CHECKBUTTON:
                        if(parm1 == 10)
                        {
                            drawGeom = cb_geom->IsOn();
                            DrawEvent(current_idx);
                        }
                        if(parm1 == 11)
                        {
                            drawHits = cb_hits->IsOn();
                            DrawEvent(current_idx);
                        }
                        if(parm1 == 12)
                        {
                            drawTrack = cb_track->IsOn();
                            DrawEvent(current_idx);
                        }
                        break;
                }
                break;
        }
        return kTRUE;
    }

    void Next()
    {
        if(current_idx < (int)events.size() - 1)
            DrawEvent(++current_idx);
    }
    void Prev()
    {
        if(current_idx > 0)
            DrawEvent(--current_idx);
    }
};

} // namespace FITALG

#endif // ANALYTICVIEWER_HH
