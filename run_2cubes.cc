// include C++ headers
#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <tuple>
#include <algorithm>
#include <vector>


// include C headers
#include <sys/stat.h>
#include <time.h>

// include ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TPaveStats.h>
#include <TError.h>
#include <TGraph.h>

using namespace std;

// include my headers
#include "const.h"
#include "ToScintiCh.h"


// EASIROCのチャンネルとホドスコープの位置・チャンネルの対応
tuple<EHodoscope, int> ToHodoscopeCh(int easirocCh, int shiftXup=0, int shiftYup=0, int shiftYDown=0, int shiftXDown=0)
{
    if (easirocCh >= 0 && easirocCh <= 15 && easirocCh - shiftXup >= 0 && easirocCh - shiftXup <= 15)
    {
        return forward_as_tuple(EHodoscope::HSX1, 16 - (easirocCh - shiftXup));
    }
    else if (easirocCh >= 16 && easirocCh <= 31 && easirocCh - shiftYup >= 16 && easirocCh - shiftYup <= 31)
    {
        return forward_as_tuple(EHodoscope::HSY1, 32 - (easirocCh - shiftYup));
    }
    else if (easirocCh >= 31 && easirocCh <= 47 && easirocCh + shiftYDown >= 32 && easirocCh + shiftYDown <= 47)
    {
        return forward_as_tuple(EHodoscope::HSY2, easirocCh + shiftYDown - 31);
    }
    else if (easirocCh >= 48 && easirocCh <= 63 && easirocCh + shiftXDown >= 48 && easirocCh + shiftXDown <= 63)
    {
        return forward_as_tuple(EHodoscope::HSX2, easirocCh + shiftXDown - 47);
    }

    #ifdef DEBUG
        cerr << "EASIROC Ch conversion error! (EASIROC -> Hodoscope)" << endl;
    #endif

    return forward_as_tuple(EHodoscope::None, 0);
}

bool isGap(int hodoX, int hodoY, int protoX, int shiftX=0, int shiftY=0)
{
    if (protoX == 3 && !(hodoX + shiftX >= 6 && hodoX + shiftX <= 11 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }
    if (protoX == 2 && !(hodoX + shiftX >= 1 && hodoX + shiftX <= 5 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }

    return false;
}

// ヒストグラム出力関数
// クラスとして作ったほうが良い：Divideする場合の処理が複雑・引数が複数パターン取れそう（1個のヒストをDrawするときはTH1*、1個以上のときはTH1**を渡さなければならないなど）
// void OutputTH1(TH1* hist, TString figName, int nHistHorizontal, int nHistVertical, int int histWidth = 800, int histHeight = 600)
// {
//     // error
//     if(nHistVertical <= 0 || nHistHorizontal <= 0) return;
//
//     TCanvas canvas = new TCanvas("canvas", "", histWidth * nHistHorizontal, histHeight * nHistVertical);
//     if(nHistHorizontal != 1 && nHistVertical != 1)
//     {
//         canvas->Divide(nHistHorizontal, nHistVertical);
//     }
//     canvas->cd(1);
//     hCrossTalkXY->Draw();
//     canvas->cd(2);
//     hCrossTalkXZ->Draw();
//     figName = TString::Format("%sCrossTalk_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
//     canvas->SaveAs(figName);
//     canvas->Clear();
// }

void changestatsBoxSize(TH1* hist, double x1, double x2, double y1, double y2)
{
    gPad->Update();
    TPaveStats *st = (TPaveStats*)hist->FindObject("stats");
    st->SetX1NDC(x1);
    st->SetX2NDC(x2);
    st->SetY1NDC(y1);
    st->SetY2NDC(y2);
}

void run_proto(int runnum, int fileCount, int shiftHSX1=0, int shiftHSY1=0, int shiftHSY2=0, int shiftHSX2=0, int min_evt=-1, int max_evt=-1, string outputFileType="png", int gap_point_0=0, int gap_pt2s_0=0, int gap_pt1s_0=0, int gap_hs_0=0, int gap_point_1=0, int gap_pt2s_1=0, int gap_pt1s_1=0, int gap_hs_1=0)
{
    // constants
    gErrorIgnoreLevel = kError;

    const array<double, 2> FitRangeCT = {0, 0.07};
    const array<double, 2> FitRangeCTDarkCut = {0, 0.1};

    const double MinPEScatterCTCenter = -10;
    const double MaxPEScatterCTCenter = 100;
    const double MinPEScatterCTLeft = -10;
    const double MaxPEScatterCTLeft = 100;

    const int hodoLowXForCT = 7;
    const int hodoHighXForCT = 10;
    const int hodoLowYForCT = 7;
    const int hodoHighYForCT = 10;

    const double DarkCutPEForCT = 0.5;

    int omit_lowx1 = 0;
    int omit_highx1 = 0;
    int omit_lowy1 = 0;
    int omit_highy1 = 0;

    int omit_lowx2 = 0;
    int omit_highx2 = 0;
    int omit_lowy2 = 0;
    int omit_highy2 = 0;

    int subrun = 4;

    // const array<int, NHodo> shiftHodo = { 0, 0, -1, 0 };
    // const int shiftHSX1 = 0;
    // const int shiftHSY1 = 0;
    // const int shiftHSY2 = -1;
    // const int shiftHSX2 = 0;


    const double HodoPEThreshold = 2.5;
    const double ProtoPEThreshold = 9.5;
    const double ProtoPEThresholdBeamAxis = 39.5;

    const TString FitOption = "Q";

    // const int NCh2Cubes = 5;
    const int NChOneSide = 5;


    // recycled variables
    TCanvas* canvas;
    TString figName;
    string histAxis, histName;

    // directory settings
    const string rootfile_dir = "../tree_root/";
    const string calibfile_dir = "calibfile/";
    string resultDirTemp = "2cubes/";

    time_t rawTime = time(nullptr);
    struct tm timeInfo = *localtime(&rawTime);
    string dateTimeStr = to_string(timeInfo.tm_sec + 100 * timeInfo.tm_min + 10000 * timeInfo.tm_hour + 1000000 * timeInfo.tm_mday + 100000000 * (timeInfo.tm_mon + 1) + 10000000000 * (1900 + timeInfo.tm_year));
    // resultDirTemp += dateTimeStr + "/";

    const string ResultDir = resultDirTemp;
    const string ScintiPEDir = ResultDir + "scintiPEEachCh/";
    const string EvtDisplayDir = ResultDir + "evtdisplay/";
    const string HodoPEDir = ResultDir + "hodoscopePEEachCh/";

    mkdir(ResultDir.c_str(), 0777);
    mkdir(ScintiPEDir.c_str(), 0777);
    mkdir(EvtDisplayDir.c_str(), 0777);
    mkdir(HodoPEDir.c_str(), 0777);


    // Histograms
    TH2D* hHodo[NHodo];
    const Double_t MinHodoMap = 0.5;
    const Double_t MaxHodoMap = 16.5;
    const Double_t MinBin = 2.5;
    const Double_t MaxBin = 30.;
    for (int i = 0; i < NHodo; i++)
    {
        histName = "h" + HodoName[i];
        histAxis = HodoName[i] + ";;;Light yield (p.e.)";
        if (HodoName[i] == "HSX1" || HodoName[i] == "HSX2")
        {
            hHodo[i] = new TH2D(histName.c_str(), histAxis.c_str(), NScifiEachHodo, MinHodoMap, MaxHodoMap, 1, 0.5, 1.5);
        }
        else
        {
            hHodo[i] = new TH2D(histName.c_str(), histAxis.c_str(), 1, 0.5, 1.5, NScifiEachHodo, MinHodoMap, MaxHodoMap);
        }
        hHodo[i]->SetMinimum(MinBin);
        hHodo[i]->SetMaximum(MaxBin);
    }

    TH1D* hHitRateHodo[NHodo];
    for (int i = 0; i < NHodo; i++)
    {
        histName = "hHitRate_" + HodoName[i];
        histAxis = "HitRate  " + HodoName[i] + ";;Number of events";
        hHitRateHodo[i] = new TH1D(histName.c_str(), histAxis.c_str(), NScifiEachHodo, MinHodoMap, MaxHodoMap);
    }

    TH1D* hPEHodo[NHodo];
    const Int_t NBinPEHodo = 70;
    const Double_t MinPEHodo = 0.0;
    const Double_t MaxPEHodo = 70.0;
    for (int i = 0; i < NHodo; i++)
    {
        histName = "hPE_" + HodoName[i];
        histAxis = "PE " + HodoName[i] + ";Light yield (p.e.);Number of events";
        hPEHodo[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEHodo, MinPEHodo, MaxPEHodo);
    }

    TH1D* hPEHodoEach[NScifi];
    for (int i = 0; i < NScifi; i++)
    {
        histName = "hPE_ch" + to_string(i);
        histAxis = "PE hodoscope ch" + to_string(i) + ";Light yield (p.e.);Number of events";
        hPEHodoEach[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEHodo, MinPEHodo, MaxPEHodo);
    }

    TH2D* hProto[NSurfaceScinti];
    const Int_t NScintiOneSide = 5;
    const Double_t MinProtoMap = 0.5;
    const Double_t MaxProtoMap = 5.5;
    const double MaxPEProto2D = 50.;
    const double MinPEProto2D = 2.5;
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        histName = "hProto_" + SurfaceName[i];
        histAxis = "Prototype : " + SurfaceName[i] + " axis;";
        if (SurfaceName[i] == "XY")
        {
            histAxis += "cube# along X; cube# along Y;Light yield (p.e.)";
        }
        else if (SurfaceName[i] == "ZY")
        {
            histAxis += "cube# along Z; cube# along Y;Light yield (p.e.)";
        }
        else if (SurfaceName[i] == "XZ")
        {
            histAxis += "cube# along X; cube# along Z;Light yield (p.e.)";
        }
        hProto[i] = new TH2D(histName.c_str(), histAxis.c_str(),  NScintiOneSide, MinProtoMap, MaxProtoMap, NScintiOneSide, MinProtoMap, MaxProtoMap);

        hProto[i]->SetMinimum(MinPEProto2D);
        hProto[i]->SetMaximum(MaxPEProto2D);
    }

    TH2D* hHitRateProto[NSurfaceScinti];
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        histName = "hHitRate_" + SurfaceName[i];
        histAxis = "HitRate " + SurfaceName[i] + ";";
        if (SurfaceName[i] == "XY")
        {
            histAxis += "cube # along X; cube# along Y;Number of events";
        }
        else if (SurfaceName[i] == "ZY")
        {
            histAxis += "cube # along Z; cube # along Y;Number of events";
        }
        else if (SurfaceName[i] == "XZ")
        {
            histAxis += "cube # along X; cube # along Z;Number of events";
        }
        hHitRateProto[i] = new TH2D(histName.c_str(), histAxis.c_str(), NScintiOneSide, MinProtoMap, MaxProtoMap, NScintiOneSide, MinProtoMap, MaxProtoMap);
    }

    TH1D* hPEProto[NSurfaceScinti];
    const double MinPEProto = 0.;
    const double MaxPEProto = 80.;
    const int NBinPEProto = 80;
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        histName = "hPE" + SurfaceName[i];
        histAxis = "PE " + SurfaceName[i] + ";Light yield (p.e.);Number of events";
        hPEProto[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);
    }

    TH2D* hHodoHitMapWithProtoHitUp = new TH2D("hHodoHitMapWithProtoHitUp", "Upstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapWithProtoHitDown = new TH2D("hHodoHitMapWithProtoHitDown", "Downstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapUp = new TH2D("hHodoHitMapUp", "Upstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapDown = new TH2D("hHodoHitMapDown", "Downstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapWithStraightBeam = new TH2D("hHodoHitMapWithStraightBeam", "Hodoscope hitmap with straight beam event;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH1D* hPEProtoEach[NSurfaceScinti][NChOneSide][NChOneSide];
    vector<tuple<EEasiroc, int, EScintiSurface, int, int>> usingProtoCh;     // 使用チャンネルを格納 <EASIROC ch, プロトタイプ面, cube # horizontal, cube # vertical>
    for (int i = 0; i < NEasiroc - 1; i++)
    {
        EEasiroc easiroc;
        if (i == 0)
            easiroc = EEasiroc::Scinti2;
        else
            easiroc = EEasiroc::Scinti1;

        ToScintiCh* toScintiCh = new ToScintiCh(EScintiType::TwoCubes);


        for (int ch = 0; ch < NChEasiroc; ch++)
        {
            if (toScintiCh->isConnected(easiroc, ch))
            {
                tuple<EScintiSurface, int, int> scintiCh = toScintiCh->getCh(easiroc, ch);

                EScintiSurface surface = get<0> (scintiCh);
                int horizontal = get<1> (scintiCh);
                int vertical = get<2> (scintiCh);

                usingProtoCh.push_back(forward_as_tuple(easiroc, ch, surface, horizontal, vertical));

                histName = "hPE";
                histAxis = "PE ";
                if (surface == EScintiSurface::XY)
                {
                    histName += "X" + to_string(horizontal) + "Y" + to_string(vertical);
                    histAxis += "X=" + to_string(horizontal) + " Y=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
                }
                else if (surface == EScintiSurface::ZY)
                {
                    histName += "Z" + to_string(horizontal) + "Y" + to_string(vertical);
                    histAxis += "Z=" + to_string(horizontal) + " Y=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
                }
                else
                {
                    histName += "X" + to_string(horizontal) + "Z" + to_string(vertical);
                    histAxis += "X=" + to_string(horizontal) + " Z=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
                }

                hPEProtoEach[static_cast<int> (surface)][horizontal][vertical] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);

                #ifdef DEBUG
                    cout << histName << " created." << endl;
                #endif
            }
        }
    }

    const int NBinCT = 100;
    const double MaxCT = 1.0;
    TH1D* hCrossTalkXZ = new TH1D("hCrossTalkXZ", "L.Y. ratio left/center (using Y readout);L.Y. left(ch40)/center(ch41);Number of events", NBinCT, 0.0, MaxCT);
    TH1D* hCrossTalkXY = new TH1D("hCrossTalkXY", "L.Y. ratio left/center (using Z readout);L.Y. left(ch8)/center(ch9);Number of events", NBinCT, 0.0, MaxCT);

    TH1D* hCrossTalkXYDarkCut = new TH1D("hCrossTalkXYDarkCut", "L.Y. ratio left/center (using Z readout, Dark count cut);L.Y. left(ch8)/center(ch9);Number of events", NBinCT, 0.0, MaxCT);
    TH1D* hCrossTalkXZDarkCut = new TH1D("hCrossTalkXZDarkCut", "L.Y. ratio left/center (using Y readout, Dark count cut);L.Y. left(ch40)/center(ch41);Number of events", NBinCT, 0.0, MaxCT);

    const int GapGraphMax = 50000 * 30;
    TH1D* hGapCount1s = new TH1D("hGapCount1s", "Gap count (hodo & proto1s);event # ;gap count", GapGraphMax, 0, GapGraphMax);
    TH1D* hGapCount2s = new TH1D("hGapCount2s", "Gap count (hodo & proto2s);event # ;gap count", GapGraphMax, 0, GapGraphMax);


    TGraph* scatterCTXY = new TGraph();
    TGraph* scatterCTXZ = new TGraph();


    // event loop variables & flags
    array<bool, NHodo> trigHodo;
    array<int, NHodo> hitCountHodo;
    array<int, NHodo> allHitCountHodo;

    array<double, NHodo> maxPEHodo;
    array<int, NHodo> maxChHodo;

    int maxChProtoXOfXY;
    int maxChProtoXOfXZ;
    double maxPEProtoXOfXY;
    double maxPEProtoXOfXZ;

    bool goodEvent;
    bool singleHit;
    bool scintiHit;
    bool isStraightBeam;
    bool goodEventForCT;

    bool goodEventForCTDownOnly;
    int countCTDownOnly = 0;
    bool goodEventForCTUpDown;
    int countCTUpDown = 0;
    bool goodEventForCTStraight;
    int countCTStraight = 0;

    array<bool, NHodo> hodoHit;

    int gapCount1s = 0;
    int gapCount2s = 0;

    int totalEvt = -1;

    string pt2sfile;
    string pt1sfile;
    string hsfile;

    double centerPEXY;
    double centerPEXZ;
    double leftPEXY;
    double leftPEXZ;
    int countCTPoint = 0;

    for (int j = runnum; j < runnum + fileCount; j++)
    {
        // input ADC data
        TString pt2sfile_name = TString::Format("%sproto2s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);
        TString pt1sfile_name = TString::Format("%sproto1s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);
        TString hsfile_name = TString::Format("%shodo_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);

        std::string pt2sfile_name_s = std::string(pt2sfile_name);
        std::string pt1sfile_name_s = std::string(pt1sfile_name);
        std::string hsfile_name_s = std::string(hsfile_name);

        string::size_type pos1 = pt2sfile_name_s.find(".root");
        string::size_type pos2 = pt1sfile_name_s.find(".root");
        string::size_type pos3 = hsfile_name_s.find(".root");
        string pt2sfile_calib(pt2sfile_name);
        pt2sfile_calib.replace(pos1, 6, "_calib");
        pt2sfile_calib = pt2sfile_calib.substr(13);
        string pt1sfile_calib(pt1sfile_name);
        pt1sfile_calib.replace(pos2, 6, "_calib");
        pt1sfile_calib = pt1sfile_calib.substr(13);
        string hsfile_calib(hsfile_name);
        hsfile_calib.replace(pos3, 6, "_calib");
        hsfile_calib = hsfile_calib.substr(13);

        pt2sfile = pt2sfile_name_s;
        pt1sfile = pt1sfile_name_s;
        hsfile = hsfile_name_s;

        Double_t mean_ped[NEasiroc][NChEasiroc];
        // Double_t mean_1pe[NEasiroc][NChEasiroc];
        Double_t gain[NEasiroc][NChEasiroc];
        // Double_t sigma_ped[NEasiroc][NChEasiroc];
        // Double_t sigma_1pe[NEasiroc][NChEasiroc];
        // Double_t chi2_ndf_ped[NEasiroc][NChEasiroc];
        // Double_t chi2_ndf_1pe[NEasiroc][NChEasiroc];

        int inputCh;
        double mp, sp, cnp, m1, s1, cn1, g;

        TString pt2s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt2sfile_calib.c_str());
        ifstream fin1(pt2s_calibname);
        while (fin1 >> inputCh >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g)
        {
            mean_ped[0][inputCh] = mp;
            // sigma_ped[0][inputCh]    = sp;
            // chi2_ndf_ped[0][inputCh] = cnp;
            // mean_1pe[0][inputCh]     = m1;
            // sigma_1pe[0][inputCh]    = s1;
            // chi2_ndf_1pe[0][inputCh] = cn1;
            gain[0][inputCh] = g;
        }

        TString pt1s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt1sfile_calib.c_str());
        ifstream fin2(pt1s_calibname);
        while (fin2 >> inputCh >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g)
        {
            mean_ped[1][inputCh] = mp;
            // sigma_ped[1][inputCh]    = sp;
            // chi2_ndf_ped[1][inputCh] = cnp;
            // mean_1pe[1][inputCh]     = m1;
            // sigma_1pe[1][inputCh]    = s1;
            // chi2_ndf_1pe[1][inputCh] = cn1;
            gain[1][inputCh] = g;
        }

        TString hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
        ifstream fin3(hs_calibname);
        while (fin3 >> inputCh >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g)
        {
            mean_ped[2][inputCh] = mp;
            // sigma_ped[2][inputCh]    = sp;
            // chi2_ndf_ped[2][inputCh] = cnp;
            // mean_1pe[2][inputCh]     = m1;
            // sigma_1pe[2][inputCh]    = s1;
            // chi2_ndf_1pe[2][inputCh] = cn1;
            gain[2][inputCh] = g;
        }

        TChain* tree1 = new TChain("tree");
        TChain* tree2 = new TChain("tree");
        TChain* tree3 = new TChain("tree");

        Int_t adc[NEasiroc][NChEasiroc];
        Double_t pe[NEasiroc][NChEasiroc];

        tree1->Add(pt2sfile.c_str());
        tree1->SetBranchAddress("ADC", &adc[0]);
        tree2->Add(pt1sfile.c_str());
        tree2->SetBranchAddress("ADC", &adc[1]);
        tree3->Add(hsfile.c_str());
        tree3->SetBranchAddress("ADC", &adc[2]);

        cout << "Protype+Hodoscope (Down): number of events = " << tree1->GetEntries() << endl;
        cout << "Hodoscope (Up): number of events = " << tree2->GetEntries() << endl;

        cout << "Select X1:" << "[" << omit_lowx1 + 1 << ", " << 16 - omit_highx1 << "]" << endl;
        cout << "Select Y1:" << "[" << omit_lowy1 + 1 << ", " << 16 - omit_highy1 << "]" << endl;
        cout << "Select Y2:" << "[" << omit_lowy2 + 1 << ", " << 16 - omit_highy2 << "]" << endl;
        cout << "Select X2:" << "[" << omit_lowx2 + 1 << ", " << 16 - omit_highx2 << "]" << endl;

        int less_evt = 0;   // 3つのEASIROCの中で最も少ないイベント数

        if (tree1->GetEntries() >= tree2->GetEntries())
        {
            less_evt = tree2->GetEntries();
        }
        else
        {
            less_evt = tree1->GetEntries();
        }

        if (less_evt >= tree3->GetEntries())
        {
            less_evt = tree3->GetEntries();
        }

        for (int evt = 0; evt < less_evt; evt++)
        {
            totalEvt++;
            #ifdef DEBUG
                if (evt % 10000 == 0)
                {
                    cout << "event: " << evt + (j - runnum) * less_evt << endl;
                }
            #endif

            if (evt < gap_point_0)
            {
                tree1->GetEntry(evt);
                tree2->GetEntry(evt);
                tree3->GetEntry(evt);
            }
            else if (evt >= gap_point_0 && evt < gap_point_1)
            {
                tree1->GetEntry(evt - gap_pt2s_0);
                tree2->GetEntry(evt - gap_pt1s_0);
                tree3->GetEntry(evt - gap_hs_0);
            }
            else if (evt >= gap_point_1)
            {
                tree1->GetEntry(evt - gap_pt2s_1);
                tree2->GetEntry(evt - gap_pt1s_1);
                tree3->GetEntry(evt - gap_hs_1);
            }


            // Initialize
            hitCountHodo = { };
            trigHodo = { };
            goodEvent = false;
            hodoHit = { };
            scintiHit = false;
            isStraightBeam = false;
            maxPEHodo = { };
            maxChHodo = { };
            maxChProtoXOfXY = 0;
            maxChProtoXOfXZ = 0;
            maxPEProtoXOfXY = 0.0;
            maxPEProtoXOfXZ = 0.0;
            goodEventForCT = false;
            goodEventForCTUpDown = false;
            goodEventForCTDownOnly = false;
            goodEventForCTStraight = false;

            for (int ch = 0; ch < NScifi; ch++)
            {
                pe[static_cast<int> (EEasiroc::Hodoscope)][ch] = (double(adc[static_cast<int> (EEasiroc::Hodoscope)][ch]) - mean_ped[2][ch]) / gain[2][ch];
                // pe[static_cast <int>(EEasiroc::Hodoscope)][ch] = double(floor(pe[static_cast <int>(EEasiroc::Hodoscope)][ch] * 10) / 10);

                tuple<EHodoscope, int> hodoCh = ToHodoscopeCh(ch, shiftHSX1, shiftHSY1, shiftHSY2, shiftHSX2);

                EHodoscope hodoPlace = get<0> (hodoCh);
                if (hodoPlace == EHodoscope::None)
                    continue;

                Int_t hodoNumber = get<1> (hodoCh);

                if (hodoPlace == EHodoscope::HSX1 || hodoPlace == EHodoscope::HSX2)
                {
                    hHodo[static_cast<int> (hodoPlace)]->SetBinContent(hodoNumber, 1, pe[static_cast<int> (EEasiroc::Hodoscope)][ch]);
                }
                else
                {
                    hHodo[static_cast<int> (hodoPlace)]->SetBinContent(1, hodoNumber, pe[static_cast<int> (EEasiroc::Hodoscope)][ch]);
                }
                // TODO: hodoscopeをomitできるようにする
                // if(pe[static_cast<int>(EEasiroc::Hodoscope)][i] > HodoPEThreshold && IsOmitted())
                if (pe[static_cast<int> (EEasiroc::Hodoscope)][ch] > HodoPEThreshold)
                {
                    hHitRateHodo[static_cast<int> (hodoPlace)]->Fill(hodoNumber);
                    hitCountHodo[static_cast<int> (hodoPlace)]++;
                    allHitCountHodo[static_cast<int> (hodoPlace)]++;
                    trigHodo[static_cast<int> (hodoPlace)] = true;
                    if (pe[static_cast<int> (EEasiroc::Hodoscope)][ch] > maxPEHodo[static_cast<int> (hodoPlace)])
                    {
                        maxChHodo[static_cast<int> (hodoPlace)] = hodoNumber;
                        maxPEHodo[static_cast<int> (hodoPlace)] = pe[static_cast<int> (EEasiroc::Hodoscope)][ch];
                    }
                    hPEHodo[static_cast<int> (hodoPlace)]->Fill(pe[static_cast<int> (EEasiroc::Hodoscope)][ch]);
                    hPEHodoEach[ch]->Fill(pe[static_cast<int> (EEasiroc::Hodoscope)][ch]);
                }
            }

            // ホドスコープ全てにヒット→good event
            if (count(trigHodo.begin(), trigHodo.end(), true) == NHodo)
            {
                goodEvent = true;
            }

            // すべてのホドスコープに1本ずつだけヒット→single hit
            // TODO: single hit イベントを飛ばすかどうか考える
            if (count(hitCountHodo.begin(), hitCountHodo.end(), 1) == NHodo)
            {
                singleHit = true;
            }
            else
            {
                // continue;
            }

            // ビームがまっすぐ飛んだ→isStraightBeam
            if (goodEvent && singleHit && maxChHodo[static_cast<int> (EHodoscope::HSX1)] == maxChHodo[static_cast<int> (EHodoscope::HSX2)] && maxChHodo[static_cast<int> (EHodoscope::HSY1)] == maxChHodo[static_cast<int> (EHodoscope::HSY2)])
            {
                isStraightBeam = true;
            }


            if (goodEvent && singleHit)
            {
                hHodoHitMapUp->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX1)], maxChHodo[static_cast<int> (EHodoscope::HSY1)]);
                hHodoHitMapDown->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
            }

            if (isStraightBeam)
            {
                hHodoHitMapWithStraightBeam->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
            }

            // Cross talkを出す用のGood eventフラグ
            if (isStraightBeam && maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXForCT && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYForCT)
            {
                goodEventForCTStraight = true;
                countCTStraight++;
            }
            if (goodEvent && singleHit &&
                maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXForCT && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYForCT &&
                maxChHodo[static_cast<int> (EHodoscope::HSX1)] >= hodoLowXForCT && maxChHodo[static_cast<int> (EHodoscope::HSX1)] <= hodoHighXForCT && maxChHodo[static_cast<int> (EHodoscope::HSY1)] >= hodoLowYForCT && maxChHodo[static_cast<int> (EHodoscope::HSY1)] <= hodoHighYForCT)
            {
                goodEventForCTUpDown = true;
                countCTUpDown++;
            }
            if (goodEvent && singleHit && maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXForCT && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYForCT && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYForCT)
            {
                goodEventForCTDownOnly = true;
                countCTDownOnly++;
            }

            goodEventForCT = goodEventForCTStraight;




            // prototype loop
            for (auto itr = usingProtoCh.begin(); itr != usingProtoCh.end(); ++itr)
            {
                EEasiroc easiroc = get<0> (*itr);
                int easirocCh = get<1> (*itr);
                EScintiSurface surface = get<2> (*itr);
                int horizontal = get<3> (*itr);
                int vertical = get<4> (*itr);

                pe[static_cast<int> (easiroc)][easirocCh] = (double(adc[static_cast<int> (easiroc)][easirocCh]) - mean_ped[static_cast<int> (easiroc)][easirocCh]) / gain[static_cast<int> (easiroc)][easirocCh];
                // pe[static_cast <int>(easiroc)][easirocCh] = double(floor(pe[static_cast <int>(easiroc)][easirocCh] * 10) / 10);

                hProto[static_cast<int> (surface)]->SetBinContent(horizontal, vertical, pe[static_cast<int> (easiroc)][easirocCh]);



                if (pe[static_cast<int> (easiroc)][easirocCh] > ProtoPEThreshold)
                {
                    hHitRateProto[static_cast<int> (surface)]->Fill(horizontal, vertical);
                    hPEProto[static_cast<int> (surface)]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
                    hPEProtoEach[static_cast<int> (surface)][horizontal][vertical]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);

                    scintiHit = true;

                    if (goodEvent && singleHit)
                    {
                        // hPEProtoGood[static_cast<int> (surface)]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
                    }

                    if (surface == EScintiSurface::XY && pe[static_cast<int> (easiroc)][easirocCh] > maxPEProtoXOfXY)
                    {
                        maxChProtoXOfXY = horizontal;
                        maxPEProtoXOfXY = pe[static_cast<int> (easiroc)][easirocCh];
                    }
                    if (surface == EScintiSurface::XZ && pe[static_cast<int> (easiroc)][easirocCh] > maxPEProtoXOfXZ)
                    {
                        maxChProtoXOfXZ = horizontal;
                        maxPEProtoXOfXZ = pe[static_cast<int> (easiroc)][easirocCh];
                    }
                }
            }

            centerPEXY = pe[static_cast<int> (EEasiroc::Scinti1)][9];
            leftPEXY = pe[static_cast<int> (EEasiroc::Scinti1)][8];
            centerPEXZ = pe[static_cast<int> (EEasiroc::Scinti2)][41];
            leftPEXZ = pe[static_cast<int> (EEasiroc::Scinti2)][40];

            // ホドスコープ内側4ch×4chヒットをクロストークのヒストグラムに使う
            if (goodEventForCT)
            {
                hCrossTalkXY->Fill(leftPEXY / centerPEXY);
                hCrossTalkXZ->Fill(leftPEXZ / centerPEXZ);

                if (leftPEXY < DarkCutPEForCT)
                {
                    hCrossTalkXYDarkCut->Fill(0);
                }
                else
                {
                    hCrossTalkXYDarkCut->Fill(leftPEXY / centerPEXY);
                }
                if (leftPEXZ < DarkCutPEForCT)
                {
                    hCrossTalkXZDarkCut->Fill(0);
                }
                else
                {
                    hCrossTalkXZDarkCut->Fill(leftPEXZ / centerPEXZ);
                }

                scatterCTXY->SetPoint(countCTPoint, centerPEXY, leftPEXY);
                scatterCTXZ->SetPoint(countCTPoint, centerPEXZ, leftPEXZ);
                countCTPoint++;
            }


            if (isGap(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)], maxChProtoXOfXY, shiftHSX2, shiftHSY2))
            {
                gapCount1s++;
            }
            if (isGap(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)], maxChProtoXOfXZ, shiftHSX2, shiftHSY2))
            {
                gapCount2s++;
            }
            hGapCount1s->SetBinContent(totalEvt, gapCount1s);
            hGapCount2s->SetBinContent(totalEvt, gapCount2s);


            // Draw event display
            // if (evt >= min_evt && evt <= max_evt && runnum == j && goodEvent && scintiHit)
            if (totalEvt >= min_evt && totalEvt <= max_evt && runnum == j && goodEvent)
            {
                const int NHistHori = 4;
                const int NHistVert = 2;
                const int HistWidth = 800;
                const int HistHeight = 600;
                canvas = new TCanvas("canvas", "", HistWidth * NHistHori, HistHeight * NHistVert);
                canvas->Divide(NHistHori, NHistVert);

                gStyle->SetOptStat(0);
                gStyle->SetPaintTextFormat("3.2f");

                const int HodoHistOrder[] = { 5, 1, 4, 8 };
                const double TitleSize = 0.03;
                const double LabelSize = 0.04;
                const double MarkerSize = 2.0;
                for (int i = 0; i < NHodo; i++)
                {
                    canvas->cd(HodoHistOrder[i]);
                    if (HodoName[i] == "HSX1" || HodoName[i] == "HSX2")
                    {
                        hHodo[i]->Draw("text60 colz");
                        hHodo[i]->GetXaxis()->SetNdivisions(NScifiEachHodo);
                        hHodo[i]->GetYaxis()->SetNdivisions(0);
                    }
                    else
                    {
                        hHodo[i]->GetXaxis()->SetNdivisions(0);
                        hHodo[i]->GetYaxis()->SetNdivisions(NScifiEachHodo);
                        hHodo[i]->Draw("text colz");
                    }
                    hHodo[i]->GetYaxis()->SetTitleSize(TitleSize);
                    hHodo[i]->GetXaxis()->SetTitleSize(TitleSize);
                    hHodo[i]->GetYaxis()->SetLabelSize(LabelSize);
                    hHodo[i]->GetXaxis()->SetLabelSize(LabelSize);
                    hHodo[i]->SetMarkerSize(MarkerSize);
                }
                const int ProtoHistOrder[] = { 3, 2, 7 };
                for (int i = 0; i < NSurfaceScinti; i++)
                {
                    canvas->cd(ProtoHistOrder[i]);
                    hProto[i]->GetYaxis()->SetTitleSize(TitleSize);
                    hProto[i]->GetXaxis()->SetTitleSize(TitleSize);
                    hProto[i]->GetYaxis()->SetLabelSize(LabelSize);
                    hProto[i]->GetXaxis()->SetLabelSize(LabelSize);
                    hProto[i]->GetXaxis()->SetNdivisions(NScintiOneSide);
                    hProto[i]->GetYaxis()->SetNdivisions(NScintiOneSide);
                    hProto[i]->SetMarkerSize(MarkerSize);
                    hProto[i]->Draw("text colz");
                }
                figName = TString::Format("%shit_%04d_%04d_evt%d.%s", EvtDisplayDir.c_str(), runnum, subrun, evt, outputFileType.c_str());
                canvas->SaveAs(figName);
                canvas->Clear();
            }
        }
    }


    // Draw Histograms
    int nHistHori = 4;
    int nHistVert = 2;
    int histWidth = 800;
    int histHeight = 600;

    array<int, NHodo> hodoHistOrder = { 5, 1, 4, 8 };
    array<int, NSurfaceScinti> protoHistOrder = { 3, 2, 7 };

    const double TitleSize = 0.03;
    const double LabelSize = 0.04;
    const double MarkerSize = 2.0;

    // hit rate
    gStyle->SetPaintTextFormat("3.0f");
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    for (int i = 0; i < NHodo; i++)
    {
        canvas->cd(hodoHistOrder[i]);
        hHitRateHodo[i]->GetYaxis()->SetTitleSize(TitleSize);
        hHitRateHodo[i]->GetXaxis()->SetTitleSize(TitleSize);
        hHitRateHodo[i]->GetYaxis()->SetLabelSize(LabelSize);
        hHitRateHodo[i]->GetXaxis()->SetLabelSize(LabelSize);
        hHitRateHodo[i]->GetXaxis()->SetNdivisions(16);
        hHitRateHodo[i]->SetMarkerSize(MarkerSize);
        hHitRateHodo[i]->GetYaxis()->SetRangeUser(0, hHitRateHodo[i]->GetMaximum() * 1.1);
        if (HodoName[i] == "HSX1" || HodoName[i] == "HSX2")
        {
            hHitRateHodo[i]->Draw();
        }
        else
        {
            hHitRateHodo[i]->Draw("hbar");
        }
    }
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        canvas->cd(protoHistOrder[i]);
        hHitRateProto[i]->GetYaxis()->SetTitleSize(TitleSize);
        hHitRateProto[i]->GetXaxis()->SetTitleSize(TitleSize);
        hHitRateProto[i]->GetYaxis()->SetLabelSize(LabelSize);
        hHitRateProto[i]->GetXaxis()->SetLabelSize(LabelSize);
        hHitRateProto[i]->GetXaxis()->SetNdivisions(5);
        hHitRateProto[i]->GetYaxis()->SetNdivisions(5);
        hHitRateProto[i]->SetMarkerSize(MarkerSize);
        hHitRateProto[i]->Draw("text colz");
    }

    figName = TString::Format("%shitRate_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    // LY (overall)
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    for (int i = 0; i < NHodo; i++)
    {
        canvas->cd(hodoHistOrder[i]);
        hPEHodo[i]->GetYaxis()->SetTitleSize(TitleSize);
        hPEHodo[i]->GetXaxis()->SetTitleSize(TitleSize);
        hPEHodo[i]->GetYaxis()->SetLabelSize(LabelSize);
        hPEHodo[i]->GetXaxis()->SetLabelSize(LabelSize);
        // hPEHodo[i]->GetXaxis()->SetNdivisions(16);
        hPEHodo[i]->SetMarkerSize(MarkerSize);
        hPEHodo[i]->Draw();
    }
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        canvas->cd(protoHistOrder[i]);
        hPEProto[i]->GetYaxis()->SetTitleSize(TitleSize);
        hPEProto[i]->GetXaxis()->SetTitleSize(TitleSize);
        hPEProto[i]->GetYaxis()->SetLabelSize(LabelSize);
        hPEProto[i]->GetXaxis()->SetLabelSize(LabelSize);
        // hPEProto[i]->GetXaxis()->SetNdivisions(5);
        // hPEProto[i]->GetYaxis()->SetNdivisions(5);
        hPEProto[i]->SetMarkerSize(MarkerSize);
        hPEProto[i]->Draw();
    }

    figName = TString::Format("%sPE_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    // Light yield (each prototype channel)
    for (auto itr = usingProtoCh.begin(); itr != usingProtoCh.end(); ++itr)
    {
        EEasiroc easiroc = get<0> (*itr);
        int easirocCh = get<1> (*itr);
        EScintiSurface surface = get<2> (*itr);
        int horizontal = get<3> (*itr);
        int vertical = get<4> (*itr);

        string chName;
        if (surface == EScintiSurface::XY)
        {
            chName = "XY=(" + to_string(horizontal) + "," + to_string(vertical) + ")";
        }
        else if (surface == EScintiSurface::ZY)
        {
            chName = "ZY=(" + to_string(horizontal) + "," + to_string(vertical) + ")";
        }
        else
        {
            chName = "XZ=(" + to_string(horizontal) + "," + to_string(vertical) + ")";
        }

        canvas = new TCanvas();
        hPEProtoEach[static_cast<int> (surface)][horizontal][vertical]->Draw();

        figName = TString::Format("%sPE_%s_%04d_%04d.%s", ScintiPEDir.c_str(), chName.c_str(), runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }


    // Light yield (each hodoscope channel)
    for (int ch = 0; ch < NScifi; ch++)
    {
        canvas = new TCanvas();

        hPEHodoEach[ch]->Draw();

        figName = TString::Format("%sPE_ch%02d_%04d_%04d.%s", HodoPEDir.c_str(), ch, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }






    // Cross talk


    nHistHori = 2;
    nHistVert = 1;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    canvas->cd(1);
    hCrossTalkXY->Draw();
    TF1* poisson = new TF1("poisson", "[0]*TMath::Poisson(x*[1], [2])", FitRangeCT[0], FitRangeCT[1]);
    poisson->SetParameters(300, 100, 0.09);
    hCrossTalkXY->Fit("poisson", "", "", FitRangeCT[0], FitRangeCT[1]);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);
    changestatsBoxSize(hCrossTalkXY, 0.7, 0.99, 0.65, 0.935);


    canvas->cd(2);
    hCrossTalkXZ->Draw();
    poisson->SetParameters(300, 100, 0.05);
    hCrossTalkXZ->Fit("poisson", "", "", FitRangeCT[0], FitRangeCT[1]);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);
    changestatsBoxSize(hCrossTalkXY, 0.7, 0.99, 0.65, 0.935);


    figName = TString::Format("%sCrossTalk_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);

    canvas->cd(1);
    hCrossTalkXYDarkCut->Draw();
    poisson->SetParameters(300, 100, 0.09);
    hCrossTalkXYDarkCut->Fit("poisson", "", "", FitRangeCT[0], FitRangeCT[1]);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);
    changestatsBoxSize(hCrossTalkXY, 0.7, 0.99, 0.65, 0.935);


    canvas->cd(2);
    hCrossTalkXZDarkCut->Draw();
    poisson->SetParameters(300, 100, 0.05);
    hCrossTalkXZDarkCut->Fit("poisson", "", "", FitRangeCT[0], FitRangeCT[1]);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);
    changestatsBoxSize(hCrossTalkXY, 0.7, 0.99, 0.65, 0.935);

    figName = TString::Format("%sCrossTalkDarkCut_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    // scatter plot of cross talk
    cout << "Straight: " << countCTStraight << endl;
    cout << "UpDown: " << countCTUpDown << endl;
    cout << "DownOnly: " << countCTDownOnly << endl;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    canvas->cd(1);
    gStyle->SetOptFit(0111);
    scatterCTXY->Fit("pol1", FitOption, "", MinPEScatterCTCenter, MaxPEScatterCTCenter);
    scatterCTXY->SetTitle("L.Y. left vs center (using Z readout)");
    scatterCTXY->GetYaxis()->SetTitle("L.Y. left (ch9) (p.e.)");
    scatterCTXY->GetXaxis()->SetTitle("L.Y. center (ch8) (p.e.)");
    scatterCTXY->GetXaxis()->SetRangeUser(MinPEScatterCTCenter, MaxPEScatterCTCenter);
    scatterCTXY->GetYaxis()->SetRangeUser(MinPEScatterCTLeft, MaxPEScatterCTLeft);
    scatterCTXY->Draw("AP");
    canvas->cd(2);
    scatterCTXZ->SetTitle("L.Y. left vs center (using Y readout)");
    scatterCTXZ->GetYaxis()->SetTitle("L.Y. left (ch40) (p.e.)");
    scatterCTXZ->GetXaxis()->SetTitle("L.Y. center (ch41) (p.e.)");
    scatterCTXZ->Fit("pol1", FitOption, "", MinPEScatterCTCenter, MaxPEScatterCTCenter);
    scatterCTXZ->GetXaxis()->SetRangeUser(MinPEScatterCTCenter, MaxPEScatterCTCenter);
    scatterCTXZ->GetYaxis()->SetRangeUser(MinPEScatterCTLeft, MaxPEScatterCTLeft);
    scatterCTXZ->Draw("AP");

    figName = TString::Format("%sCrossTalkScatterPlot_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();





    // gap count
    nHistHori = 2;
    nHistVert = 1;
    histWidth = 800;
    histHeight = 600;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);

    hGapCount1s->GetXaxis()->SetRangeUser(0, totalEvt);
    hGapCount2s->GetXaxis()->SetRangeUser(0, totalEvt);
    hGapCount1s->GetXaxis()->SetNdivisions(fileCount, kFALSE);
    hGapCount2s->GetXaxis()->SetNdivisions(fileCount, kFALSE);
    hGapCount1s->GetXaxis()->SetTickLength(1);
    hGapCount2s->GetXaxis()->SetTickLength(1);

    canvas->cd(1);
    hGapCount1s->Draw("P");
    canvas->cd(2);
    hGapCount2s->Draw("P");
    figName = TString::Format("%sGapCount_%04d_%04d-%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, runnum + fileCount - 1, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();






    // hodoscope hit map
    nHistHori = 3;
    nHistVert = 1;
    histWidth = 600;
    histHeight = 600;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);

    canvas->cd(1);
    // hHodoHitMapUp->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapUp->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapUp->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapUp->GetYaxis()->SetNdivisions(NScifiEachHodo);
    // hHodoHitMapUp->SetMarkerSize(2.0);
    hHodoHitMapUp->Draw("text colz");

    canvas->cd(2);
    // hHodoHitMapDown->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapDown->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapDown->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapDown->GetYaxis()->SetNdivisions(NScifiEachHodo);
    // hHodoHitMapDown->SetMarkerSize(2.0);
    hHodoHitMapDown->Draw("text colz");

    canvas->cd(3);
    hHodoHitMapWithStraightBeam->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->Draw("text colz");

    figName = TString::Format("%sHodoHitMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    //
}

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        cerr << "run_proto <runnum> <# of files> <shiftHSX1> <shiftHSY1> <shiftHSY2> <shiftHSX2> <min_evt> <max_evt>" << endl;

        return -1;
    }

    int runnum = atoi(argv[1]);
    int fileCount = atoi(argv[2]);
    int x1, y1, y2, x2;
    if (argc == 3 + 4)
    {
        x1 = atoi(argv[3]);
        y1 = atoi(argv[4]);
        y2 = atoi(argv[5]);
        x2 = atoi(argv[6]);
    }
    int min_evt, max_evt;
    if (argc == 3 + 4 + 2)
    {
        min_evt = atoi(argv[7]);
        max_evt = atoi(argv[8]);
    }

    // int gap_point_0 = atoi(argv[5]);
    // int gap_pt2s_0  = atoi(argv[6]);
    // int gap_pt1s_0  = atoi(argv[7]);
    // int gap_hs_0    = atoi(argv[8]);
    // int gap_point_1 = atoi(argv[9]);
    // int gap_pt2s_1  = atoi(argv[10]);
    // int gap_pt1s_1  = atoi(argv[11]);
    // int gap_hs_1    = atoi(argv[12]);


    if (argc == 3)
    {
        run_proto(runnum, fileCount);
    }
    else if (argc == 3 + 4)
    {
        run_proto(runnum, fileCount, x1, y1, y2, x2);
    }
    else if (argc == 3 + 4 + 2)
    {
        run_proto(runnum, fileCount, x1, y1, y2, x2, min_evt, max_evt);
    }


    return 0;
}
