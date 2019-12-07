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
#include <boost/format.hpp>


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
#include <TLine.h>
#include <TEllipse.h>

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

void SaveHist(TH1* hist, TString outputFileDir, TString drawOption = "", bool setLogy = false, int histWidth = 0, int histHeight = 0)
{
    TCanvas* canvas;
    if (histWidth == 0 || histHeight == 0)
    {
        canvas = new TCanvas();
    }
    else
    {
        canvas = new TCanvas("canvas", "", histWidth, histHeight);
    }

    if (setLogy)
    {
        canvas->SetLogy();
    }
    hist->Draw(drawOption);
    canvas->SaveAs(outputFileDir);
    canvas->Clear();
}

void SaveHodoMap(TH2* hist, TString outputFileDir, int nCellOneSide)
{
    TCanvas* canvas = new TCanvas("canvas", "", 1280, 1200);
    hist->Draw("text colz");
    hist->GetXaxis()->SetNdivisions(nCellOneSide);
    hist->GetYaxis()->SetNdivisions(nCellOneSide);
    gPad->SetRightMargin(0.17);
    hist->GetZaxis()->SetTitleOffset(1.4);
    hist->SetStats(kFALSE);
    canvas->SaveAs(outputFileDir);
    canvas->Clear();
}

void changestatsBoxSize(TH1* hist, double x1, double x2, double y1, double y2)
{
    gPad->Update();
    TPaveStats* st = (TPaveStats*) hist->FindObject("stats");
    st->SetX1NDC(x1);
    st->SetX2NDC(x2);
    st->SetY1NDC(y1);
    st->SetY2NDC(y2);
}

// lineColor:
// 0:white, 1:black, 2:red, 3:green, 4:blue, 5:yellow, 6:magenta, 7:cyan, 8:dark green, 9:purple
void drawCubeLine(string config="", int lineColor=7)
{
    //     const double binmin = -0.1;
    // const double binmax = 30.;

    // Physical length (mm)
    const double HoleRadius = 0.75;
    const double HoleCenterPosFromCubeEdge = 3.0;
    const double CubeSize = 10.;
    const double SciFiWidth = 1.7;

    const double low_x = 0.1 / SciFiWidth + 5.5;
    const double upp_x = low_x + CubeSize / SciFiWidth;
    const double low_y = low_x;
    const double upp_y = low_y + CubeSize / SciFiWidth;
    const double low_fiber_pos = upp_x - (HoleCenterPosFromCubeEdge - HoleRadius) / SciFiWidth;
    const double upp_fiber_pos = upp_x - (HoleCenterPosFromCubeEdge + HoleRadius) / SciFiWidth;
    const double Zradius = HoleRadius / SciFiWidth;
    const double center_pos = HoleCenterPosFromCubeEdge / SciFiWidth;

    // const int lineColor = 7;
    // 0:white, 1:black, 2:red, 3:green, 4:blue, 5:yellow, 6:magenta, 7:cyan, 8:dark green, 9:purple
    const int LineStyle = 2;
    // 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    const int LineWidth = 8;

    double xShift = 0;
    double yShift = 0;


    if (config == "inj")
    {
        xShift = -0.3;
        yShift = 0.7;
    }
    else if (config == "ent")
    {
        xShift = 0;
        yShift = 0;
    }

    // キューブ左側
    TLine* line1 = new TLine(low_x + xShift, low_y + yShift, low_x + xShift, upp_y + yShift);
    line1->SetLineColor(lineColor);
    line1->SetLineWidth(LineWidth);
    line1->SetLineStyle(LineStyle);
    // キューブ右側
    TLine* line2 = new TLine(upp_x + xShift, low_y + yShift, upp_x + xShift, upp_y + yShift);
    line2->SetLineColor(lineColor);
    line2->SetLineWidth(LineWidth);
    line2->SetLineStyle(LineStyle);
    // キューブ上側
    TLine* line3 = new TLine(low_x + xShift, upp_y + yShift, upp_x + xShift, upp_y + yShift);
    line3->SetLineColor(lineColor);
    line3->SetLineWidth(LineWidth);
    line3->SetLineStyle(LineStyle);
    // キューブ下側
    TLine* line4 = new TLine(low_x + xShift, low_y + yShift, upp_x + xShift, low_y + yShift);
    line4->SetLineColor(lineColor);
    line4->SetLineWidth(LineWidth);
    line4->SetLineStyle(LineStyle);
    // ファイバー縦方向の左側
    TLine* line5 = new TLine(low_fiber_pos + xShift, low_y + yShift, low_fiber_pos + xShift, upp_y + yShift);
    line5->SetLineColor(lineColor);
    line5->SetLineWidth(LineWidth);
    line5->SetLineStyle(LineStyle);
    // ファイバー縦方向の右側
    TLine* line6 = new TLine(upp_fiber_pos + xShift, low_y + yShift, upp_fiber_pos + xShift, upp_y + yShift);
    line6->SetLineColor(lineColor);
    line6->SetLineWidth(LineWidth);
    line6->SetLineStyle(LineStyle);
    // ファイバー横方向の下側
    TLine* line7 = new TLine(low_x + xShift, low_fiber_pos + yShift, upp_x + xShift, low_fiber_pos + yShift);
    line7->SetLineColor(lineColor);
    line7->SetLineWidth(LineWidth);
    line7->SetLineStyle(LineStyle);
    // ファイバー横方向の上側
    TLine* line8 = new TLine(low_x + xShift, upp_fiber_pos + yShift, upp_x + xShift, upp_fiber_pos + yShift);
    line8->SetLineColor(lineColor);
    line8->SetLineWidth(LineWidth);
    line8->SetLineStyle(LineStyle);

    TEllipse* circle = new TEllipse(low_x + center_pos + xShift, low_y + center_pos + yShift, Zradius, Zradius);
    circle->SetLineColor(lineColor);
    circle->SetLineWidth(4);
    circle->SetLineStyle(2);
    // circle->SetLineWidth(LineWidth);
    // circle->SetLineStyle(LineStyle);
    circle->SetFillColorAlpha(0, 0);


    line1->Draw();
    line2->Draw();
    line3->Draw();
    line4->Draw();
    line5->Draw();
    line6->Draw();
    line7->Draw();
    line8->Draw();
    circle->Draw();
}

void run_proto(int runnum, int fileCount, int shiftHSX1=0, int shiftHSY1=0, int shiftHSY2=0, int shiftHSX2=0, int min_evt=-1, int max_evt=-1, string outputFileType="png", int gap_point_0=0, int gap_pt2s_0=0, int gap_pt1s_0=0, int gap_hs_0=0, int gap_point_1=0, int gap_pt2s_1=0, int gap_pt1s_1=0, int gap_hs_1=0)
{
    // constants
    gErrorIgnoreLevel = kError;

    const array<double, 2> FitRangeCT = { 0, 0.16 };
    const array<double, 2> FitRangeCTDarkCut = { 0, 0.1 };
    const array<double, 2> FitRangePECenterXY = { 30, 55 };
    const array<double, 2> FitRangePECenterXZ = { 20, 40 };
    const array<double, 2> FitRangePELeftXY = { -0.5, 5.5 };
    const array<double, 2> FitRangePELeftXZ = { -0.5, 5.5 };

    const double MinPEScatterCTCenter = -10;
    const double MaxPEScatterCTCenter = 100;
    const double MinPEScatterCTAround = -10;
    const double MaxPEScatterCTAround = 100;

    const int hodoLowXForCT = 7;
    const int hodoHighXForCT = 10;
    const int hodoLowYForCT = 7;
    const int hodoHighYForCT = 10;

    // const int hodoLowXForCT = 8;
    // const int hodoHighXForCT = 9;
    // const int hodoLowYForCT = 8;
    // const int hodoHighYForCT = 9;

    const int hodoLowXForCTCell = 1;
    const int hodoHighXForCTCell = 16;
    const int hodoLowYForCTCell = 1;
    const int hodoHighYForCTCell = 16;
    //
    // const int hodoLowXForCTCell = 7;
    // const int hodoHighXForCTCell = 10;
    // const int hodoLowYForCTCell = 7;
    // const int hodoHighYForCTCell = 10;

    const double DarkCutPEForCT = 0.5;

    const double RatioCutForCT = 1.0;

    int omit_lowx1 = 0;
    int omit_highx1 = 0;
    int omit_lowy1 = 0;
    int omit_highy1 = 0;

    int omit_lowx2 = 0;
    int omit_highx2 = 0;
    int omit_lowy2 = 0;
    int omit_highy2 = 0;

    int subrun = 5;


    const double HodoPEThreshold = 2.5;
    const double ProtoPEThreshold = 5.5;
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
    string resultDirTemp = "9cubes/";
    mkdir(resultDirTemp.c_str(), 0777);

    time_t rawTime = time(nullptr);
    struct tm timeInfo = *localtime(&rawTime);
    string dateTimeStr = to_string(timeInfo.tm_sec + 100 * timeInfo.tm_min + 10000 * timeInfo.tm_hour + 1000000 * timeInfo.tm_mday + 100000000 * (timeInfo.tm_mon + 1) + 10000000000 * (1900 + timeInfo.tm_year));
    // resultDirTemp += dateTimeStr + "/";

    const string ResultDir = resultDirTemp;
    const string ScintiPEDir = ResultDir + "scintiPEEachCh/";
    const string EvtDisplayDir = ResultDir + "evtdisplay/";
    const string HodoPEDir = ResultDir + "hodoscopePEEachCh/";
    const string CellPEDir = ResultDir + "scintiPEEachCell/";
    const string CrosstalkDir = ResultDir + "crosstalkEachCell/";

    mkdir(ResultDir.c_str(), 0777);
    mkdir(ScintiPEDir.c_str(), 0777);
    mkdir(EvtDisplayDir.c_str(), 0777);
    mkdir(HodoPEDir.c_str(), 0777);
    mkdir(CellPEDir.c_str(), 0777);
    mkdir(CrosstalkDir.c_str(), 0777);


    // Histograms
    TH2D* hHodo[NHodo];
    const double MinHodoMap = 0.5;
    const double MaxHodoMap = 16.5;
    const double MinBin = 2.5;
    const double MaxBin = 30.;
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
    const double MinPEHodo = 0.0;
    const double MaxPEHodo = 70.0;
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
    const double MinProtoMap = 0.5;
    const double MaxProtoMap = 5.5;
    const double MaxPEProto2D = 50.;
    const double MinPEProto2D = 0;
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
    const double MinPEProto = 2.5;
    const double MaxPEProto = 89.5;
    const int NBinPEProto = 92;
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        histName = "hPE" + SurfaceName[i];
        histAxis = "PE " + SurfaceName[i] + ";Light yield (p.e.);Number of events";
        hPEProto[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);
    }

    const double MinPECenter = -1.5;
    const double MaxPECenter = 98.5;
    const int NBinPECenter = MaxPECenter - MinPECenter;
    TH1D* hPECenterForCTXY = new TH1D("hPECenterForCTXY", "PE center (using Z readout);Light yield (ch9) (p.e.);Number of events", NBinPECenter, MinPECenter, MaxPECenter);


    TH2D* hHodoHitMapWithProtoHitUp = new TH2D("hHodoHitMapWithProtoHitUp", "Upstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapWithProtoHitDown = new TH2D("hHodoHitMapWithProtoHitDown", "Downstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapUp = new TH2D("hHodoHitMapUp", "Upstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapDown = new TH2D("hHodoHitMapDown", "Downstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapWithStraightBeam = new TH2D("hHodoHitMapWithStraightBeam", "Hodoscope hitmap with straight beam event;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH1D* hPEProtoEach[NSurfaceScinti][NChOneSide][NChOneSide];

    vector<tuple<EEasiroc, int, EScintiSurface, int, int>> usingProtoCh;    // 使用チャンネルを格納 <EASIROC ch, プロトタイプ面, cube # horizontal, cube # vertical>
    vector<tuple<EEasiroc, int, EScintiSurface, int, int>> usingProtoChForAllOutput;
    for (int i = 0; i < NEasiroc - 1; i++)
    {
        EEasiroc easiroc;
        if (i == 0)
            easiroc = EEasiroc::Scinti2;
        else
            easiroc = EEasiroc::Scinti1;

        ToScintiCh* toScintiCh = new ToScintiCh(EScintiType::NineCubes);
        ToScintiCh* toScintiCh555 = new ToScintiCh(EScintiType::Proto555);


        for (int ch = 0; ch < NChEasiroc; ch++)
        {
            if(toScintiCh->isConnected(easiroc, ch))
            {
                tuple<EScintiSurface, int, int> scintiCh = toScintiCh->getCh(easiroc, ch);
                EScintiSurface surface = get<0> (scintiCh);
                int horizontal = get<1> (scintiCh);
                int vertical = get<2> (scintiCh);
                usingProtoCh.push_back(forward_as_tuple(easiroc, ch, surface, horizontal, vertical));
            }
            if(toScintiCh555->isConnected(easiroc, ch))
            {
                tuple<EScintiSurface, int, int> scintiCh = toScintiCh555->getCh(easiroc, ch);

                EScintiSurface surface = get<0> (scintiCh);
                int horizontal = get<1> (scintiCh);
                int vertical = get<2> (scintiCh);

                usingProtoChForAllOutput.push_back(forward_as_tuple(easiroc, ch, surface, horizontal, vertical));
                if(horizontal == 0 || vertical == 0)
                    continue;
                // cout << horizontal << ", " << vertical << endl;

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

                hPEProtoEach[static_cast<int> (surface)][horizontal-1][vertical-1] = new TH1D(histName.c_str(), histAxis.c_str(), 10, -1.5, 8.5);

                #ifdef DEBUG
                    cout << histName << " created." << endl;
                #endif
            }
        }
    }


    // Cross talk related
    const int NBinCT = 110;
    const double MinCT = -0.1;
    const double MaxCT = 1.0;
    const int NCubeCT = 8;
    const double MinCTMap = 0;
    const double MaxCTMap = 10;
    const int CenterCubeXOfXY = 3;
    const int CenterCubeYOfXY = 3;
    // For scatter histograms
    const double NBinCTScatter = 100;
    const double MinCTScatter = -4.5;
    const double MaxCTScatter = 95.5;
    const double MinPEAround = -1.5;
    const double MaxPEAround = 18.5;
    const int NBinPEAround = MaxPEAround - MinPEAround;
    // 周辺キューブ用配列添字とChの対応
    // 0 1 2
    // 3 x 4
    // 5 6 7
    const int CubeChMapXY[] = { 11, 12, 28, 8, 26, 5, 6, 24 };
    const string CubeGeometryName[] = { "UpperLeft", "Upper", "UpperRight", "Left", "Right", "LowerLeft", "Lower", "LowerRight" };
    const string CubeGeometryTitle[] = { "upper left", "upper", "upper right", "left", "right", "lower left", "lower", "lower right" };
    array<TH1D*, NCubeCT> hCrosstalkXY;
    array<TH1D*, NCubeCT> hCrosstalkXYDarkCut;

    array<TH1D*, NCubeCT> hPEAroundForCTXY;
    TH1D* hCrosstalkXYDarkCutEachCell[NCubeCT][NScifiEachHodo][NScifiEachHodo];
    array<TH2D*, NCubeCT> hCrosstalkXYDarkCutMap;

    array<TGraph*, NCubeCT> scatterCTXY;
    array<TH2D*, NCubeCT> hCrosstalkScatterXY;
    TH2D* hCrosstalkScatterXYEachCell[NCubeCT][NScifiEachHodo][NScifiEachHodo];
    for (int i = 0; i < NCubeCT; ++i)
    {
        histName = "hCrosstalkXY" + CubeGeometryName[i];
        histAxis = "L.Y. ratio " + CubeGeometryTitle[i] + "/center (using Z readout);L.Y. " + CubeGeometryTitle[i] + "(ch" + to_string(CubeChMapXY[i]) + ")/center(ch9);Number of events";
        hCrosstalkXY[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinCT, MinCT, MaxCT);

        histName += "DarkCut";
        histAxis = "L.Y. ratio " + CubeGeometryTitle[i] + "/center (using Z readout, dark count cut);L.Y. " + CubeGeometryTitle[i] + "(ch" + to_string(CubeChMapXY[i]) + ")/center(ch9);Number of events";
        hCrosstalkXYDarkCut[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinCT, MinCT, MaxCT);


        histName = "hPE" + CubeGeometryName[i] + "ForCTXY";
        // histAxis = "PE left (using Z readout);Light yield (ch8) (p.e.);Number of events";
        histAxis = (boost::format("PE %s (using Z readout);Light yield (ch%d) (p.e.);Number of events") % CubeGeometryTitle[i] % CubeChMapXY[i]).str();
        hPEAroundForCTXY[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEAround, MinPEAround, MaxPEAround);
        for (int j = 0; j < NScifiEachHodo; j++)
        {
            for (int k = 0; k < NScifiEachHodo; k++)
            {
                histName = (boost::format("hCrosstalkXY%sDarkCutX%dY%d") % CubeGeometryName[i] % (j + 1) % (k + 1)).str();
                histAxis = (boost::format("L.Y. ratio %s/center (using Z readout, dark count cut, Cell X=%d Y=%d;L.Y. %s(ch%d)/center(ch9);Number of events") % CubeGeometryTitle[i] % (j + 1) % (k + 1) % CubeGeometryTitle[i] % CubeChMapXY[i]).str();
                hCrosstalkXYDarkCutEachCell[i][j][k] = new TH1D(histName.c_str(), histAxis.c_str(), NBinCT, MinCT, MaxCT);

                histName = (boost::format("hCrosstalk%sScatterXYX%dY%d") % CubeGeometryName[i] % (j + 1) % (k + 1)).str();
                histAxis = (boost::format("L.Y. %s vs center (using Z readout, cell X=%d Y=%d);L.Y. center (ch8) (p.e.);L.Y. %s (ch%d) (p.e.);Number of events") % CubeGeometryTitle[i] % (j + 1) % (k + 1) % CubeGeometryTitle[i] % CubeChMapXY[i]).str();
                hCrosstalkScatterXYEachCell[i][j][k] = new TH2D(histName.c_str(), histAxis.c_str(), NBinCTScatter, MinCTScatter, MaxCTScatter, NBinCTScatter, MinCTScatter, MaxCTScatter);
            }
        }

        histName = (boost::format("hCrosstalkXY%sDarkCutMap") % CubeGeometryName[i]).str();
        histAxis = (boost::format("Crosstalk ratio (%s, using Z readout, dark count cut);cell # along X;cell # along Y;Crosstalk ratio (%%)") % CubeGeometryTitle[i]).str();
        hCrosstalkXYDarkCutMap[i] = new TH2D(histName.c_str(), histAxis.c_str(), NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
        hCrosstalkXYDarkCutMap[i]->SetMinimum(MinCTMap);
        hCrosstalkXYDarkCutMap[i]->SetMaximum(MaxCTMap);

        scatterCTXY[i] = new TGraph();

        histName = "hCrosstalkScatterXY" + CubeGeometryName[i];
        histAxis = (boost::format("L.Y. %s vs center (using Z readout);L.Y. center (ch8) (p.e.);L.Y. %s (ch%d) (p.e.);Number of events") % CubeGeometryTitle[i] % CubeGeometryTitle[i] % CubeChMapXY[i]).str();
        hCrosstalkScatterXY[i] = new TH2D(histName.c_str(), histAxis.c_str(), NBinCTScatter, MinCTScatter, MaxCTScatter, NBinCTScatter, MinCTScatter, MaxCTScatter);
    }


    // 対角のキューブのクロストーク見る用
    const int NCubeCTDiagOpp = 4;
    const int CubeChMapXYDiagOpp[NCubeCTDiagOpp][2] = { { 11, 24 }, { 24, 11 }, { 28, 5 }, { 5, 28 } };
    const string CubeGeometryNameDiagOpp[] = { "LR/UL", "UL/LR", "LL/UR", "UR/LL" };

    const string CubeGeometryTitleDiagOppNumera[] = { "lower right", "upper left", "lower left", "upper right" };
    const string CubeGeometryTitleDiagOppDenomi[] = { "upper left", "lower right", "upper right", "lower left" };
    array<TH1D*, NCubeCTDiagOpp> hCrosstalkXYDarkCutDiagOpp;   // 対角方向の2つのキューブの光量の割合

    array<TH1D*, NCubeCTDiagOpp> hPENoHitForCTXYDiagOpp;
    array<TH1D*, NCubeCTDiagOpp> hPEHitForCTXYDiagOpp;

    TH1D* hCrosstalkXYDarkCutEachCellDiagOpp[NCubeCTDiagOpp][NScifiEachHodo][NScifiEachHodo];
    array<TH2D*, NCubeCTDiagOpp> hCrosstalkXYDarkCutMapDiagOpp;

    array<TGraph*, NCubeCTDiagOpp> scatterCTXYDiagOpp;
    array<TH2D*, NCubeCTDiagOpp> hCrosstalkScatterXYDiagOpp;
    TH2D* hCrosstalkScatterXYEachCellDiagOpp[NCubeCTDiagOpp][NScifiEachHodo][NScifiEachHodo];

    for (int i = 0; i < NCubeCTDiagOpp; i++)
    {
        histName = "hCrosstalkXY" + CubeGeometryNameDiagOpp[i] + "DarkCut";
        histAxis = (boost::format("L.Y. ratio %1%/%2% (using Z readout, dark count cut);L.Y. %1%(ch%3%)/%2%(ch%4%);Number of events") % CubeGeometryTitleDiagOppNumera[i] % CubeGeometryTitleDiagOppDenomi[i] % CubeChMapXYDiagOpp[i][1] % CubeChMapXYDiagOpp[i][0]).str();
        // histAxis = "L.Y. ratio " + CubeGeometryTitleDiagOppNumera[i] + " / " + CubeGeometryTitleDiagOppDenomi[i] + " (using Z readout);L.Y. " + CubeGeometryTitleDiagOppNumera[i] + "(ch" + to_string(CubeChMapXY[i][1]) + ")/" + CubeGeometryTitleDiagOppDenomi[i] + "(ch9);Number of events";
        hCrosstalkXYDarkCutDiagOpp[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinCT, MinCT, MaxCT);

        // クロストーク光量/ビームヒットキューブの光量のうち，クロストークキューブの光量分布
        histName = "hPENumeraOf" + CubeGeometryNameDiagOpp[i];
        histAxis = (boost::format("PE %1% (beam hit to %2% cube) ;L.Y. %1% (ch %3%);Number of events") % CubeGeometryTitleDiagOppNumera[i] % CubeGeometryTitleDiagOppDenomi[i] % CubeChMapXYDiagOpp[i][1]).str();
        hPENoHitForCTXYDiagOpp[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEAround, MinPEAround, MaxPEAround);

        // クロストーク光量/ビームヒットキューブの光量のうち，ビームヒットキューブの光量分布
        histName = "hPEDenomiOf" + CubeGeometryNameDiagOpp[i];
        histAxis = (boost::format("PE %1% (beam hit to %1% cube) ;L.Y. %1% (ch %2%);Number of events") % CubeGeometryTitleDiagOppDenomi[i] % CubeChMapXYDiagOpp[i][0]).str();
        hPEHitForCTXYDiagOpp[i] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPECenter, MinPECenter, MaxPECenter);


        for (int j = 0; j < NScifiEachHodo; j++)
        {
            for (int k = 0; k < NScifiEachHodo; k++)
            {
                histName = (boost::format("hCrosstalkXY%sDarkCutX%dY%d") % CubeGeometryNameDiagOpp[i] % (j + 1) % (k + 1)).str();
                histAxis = (boost::format("L.Y. ratio %1%/%2% (using Z readout, dark count cut, Cell X=%3% Y=%4%;L.Y. %1%(ch%5%)/%2%(ch%6%);Number of events") % CubeGeometryTitleDiagOppNumera[i] % CubeGeometryTitleDiagOppDenomi[i] % (j + 1) % (k + 1) % CubeChMapXYDiagOpp[i][1] % CubeChMapXYDiagOpp[i][0]).str();
                hCrosstalkXYDarkCutEachCellDiagOpp[i][j][k] = new TH1D(histName.c_str(), histAxis.c_str(), NBinCT, MinCT, MaxCT);

                // histName = (boost::format("hCrosstalk%sScatterXYX%dY%d") % CubeGeometryName[i] % (j + 1) % (k + 1)).str();
                // histAxis = (boost::format("L.Y. %s vs center (using Z readout, cell X=%d Y=%d);L.Y. center (ch8) (p.e.);L.Y. %s (ch%d) (p.e.);Number of events") % CubeGeometryTitle[i] % (j + 1) % (k + 1) % CubeGeometryTitle[i] % CubeChMapXY[i]).str();
                // hCrosstalkScatterXYEachCell[i][j][k] = new TH2D(histName.c_str(), histAxis.c_str(), NBinCTScatter, MinCTScatter, MaxCTScatter, NBinCTScatter, MinCTScatter, MaxCTScatter);
            }
        }

        histName = (boost::format("hCrosstalkXY%sDarkCutMap") % CubeGeometryNameDiagOpp[i]).str();
        histAxis = (boost::format("Crosstalk ratio (%1%/%2%, using Z readout, dark count cut);cell # along X;cell # along Y;Crosstalk ratio (%%)") % CubeGeometryTitleDiagOppNumera[i] % CubeGeometryTitleDiagOppDenomi[i]).str();
        hCrosstalkXYDarkCutMapDiagOpp[i] = new TH2D(histName.c_str(), histAxis.c_str(), NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
        hCrosstalkXYDarkCutMapDiagOpp[i]->SetMinimum(MinCTMap);
        hCrosstalkXYDarkCutMapDiagOpp[i]->SetMaximum(MaxCTMap);

        scatterCTXYDiagOpp[i] = new TGraph();

        histName = "hCrosstalkScatterXY" + CubeGeometryNameDiagOpp[i];
        histAxis = (boost::format("L.Y. %1% vs %2% (using Z readout);L.Y. %1% (ch%3%) (p.e.);L.Y. %2% (ch%4%) (p.e.);Number of events") % CubeGeometryTitleDiagOppDenomi[i] % CubeGeometryTitleDiagOppNumera[i] % CubeChMapXYDiagOpp[i][0] % CubeChMapXYDiagOpp[i][1]).str();
        hCrosstalkScatterXYDiagOpp[i] = new TH2D(histName.c_str(), histAxis.c_str(), NBinCTScatter, MinCTScatter, MaxCTScatter, NBinCTScatter, MinCTScatter, MaxCTScatter);
    }

    // Detection efficiency
    TH2D* hDetectionEff[NSurfaceScinti];
    TH2D* hDetectionEffZoom[NSurfaceScinti];
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        histName = "hDetectionEfficiency" + SurfaceName[i];
        histAxis = (boost::format("Detection efficiency (using %s readout);Cell # along X;Cell # along Y;Detection Efficiency") % ReadoutSurfaceName[i]).str();
        hDetectionEff[i] = new TH2D(histName.c_str(), histAxis.c_str(), NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
        hDetectionEff[i]->SetMinimum(0.0);
        hDetectionEff[i]->SetMaximum(1.0);
    }

    // Hit multiplicity
    TH1D* hHitMultiplicity[NHodo];
    for (int i=0; i<NHodo; ++i)
    {
        histName = "hHitMultiplicity" + HodoName[i];
        histAxis = (boost::format("Number of hits per event (%s);Number of hits;Number of events") % HodoName[i]).str();
        hHitMultiplicity[i] = new TH1D(histName.c_str(), histAxis.c_str(), 10, -0.5, 9.5);
    }

    // Gap check
    const int GapGraphMax = 50000 * 30;
    TH1D* hGapCount1s = new TH1D("hGapCount1s", "Gap count (hodo & proto1s);event # ;gap count", GapGraphMax, 0, GapGraphMax);
    TH1D* hGapCount2s = new TH1D("hGapCount2s", "Gap count (hodo & proto2s);event # ;gap count", GapGraphMax, 0, GapGraphMax);



    // event loop variables & flags
    array<bool, NHodo> trigHodo;
    array<int, NHodo> hitCountHodo;
    // array<int, NHodo> allHitCountHodo;

    array<double, NHodo> maxPEHodo;
    array<int, NHodo> maxChHodo;

    int maxChProtoXOfXY;
    int maxChProtoYOfXY;
    int maxChProtoXOfXZ;
    double maxPEProtoXY;
    double maxPEProtoXOfXZ;

    bool goodEvent;
    bool singleHit;
    // bool scintiHit;
    array<bool, 3> scintiHitEachSurface;
    bool isStraightBeam;
    bool goodEventForCT;

    // Cross talkヒストグラムの各選別条件に対するフラグ
    bool goodEventForCTDownOnly;
    int countCTDownOnly = 0;
    bool goodEventForCTUpDown;
    int countCTUpDown = 0;
    bool goodEventForCTStraight;
    int countCTStraight = 0;

    // Cellごとのcross talkヒストグラムに対するフラグ
    bool goodEventForCTCell;
    int countCTCell = 0;

    array<bool, NCubeCTDiagOpp> goodEventForCTDiagOpp;
    array<int, NCubeCTDiagOpp> countCTPointDiagOpp = { };


    // Detection eff.を出すための、読み出し方向ごとのgood event数・すべてのgood event数
    // ここでのgood eventとは、藤田さんの修論P30の定義
    int countHodoHitEachCell[NScifiEachHodo][NScifiEachHodo] = { };
    int countScintiHitEachSurfaceCell[NSurfaceScinti][NScifiEachHodo][NScifiEachHodo] = { };
    int countStraight = 0;
    // array<int, NSurfaceScinti> countScintiHitEachSurface = { };

    array<bool, NHodo> hodoHit;

    int gapCount1s = 0;
    int gapCount2s = 0;

    int totalEvt = -1;

    string pt2sfile;
    string pt1sfile;
    string hsfile;

    double centerPEXY;
    array<double, NCubeCT> aroundPEXY;
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

        double mean_ped[NEasiroc][NChEasiroc];
        // double mean_1pe[NEasiroc][NChEasiroc];
        double gain[NEasiroc][NChEasiroc];
        // double sigma_ped[NEasiroc][NChEasiroc];
        // double sigma_1pe[NEasiroc][NChEasiroc];
        // double chi2_ndf_ped[NEasiroc][NChEasiroc];
        // double chi2_ndf_1pe[NEasiroc][NChEasiroc];

        int inputCh;
        double mp, sp, cnp, m1, s1, cn1, g;

        TString pt2s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt2sfile_calib.c_str());
        ifstream fin1(pt2s_calibname);
        if (fin1.fail())
        {
            cerr << "Calibration file does not exist!!!" << endl;
            return;
        }
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
        double pe[NEasiroc][NChEasiroc];

        tree1->Add(pt2sfile.c_str());
        tree1->SetBranchAddress("ADC", &adc[0]);
        tree2->Add(pt1sfile.c_str());
        tree2->SetBranchAddress("ADC", &adc[1]);
        tree3->Add(hsfile.c_str());
        tree3->SetBranchAddress("ADC", &adc[2]);

        cout << "Protype (XY): number of events = " << tree1->GetEntries() << endl;
        cout << "Protype (XZ & ZY): number of events = " << tree2->GetEntries() << endl;
        cout << "Hodoscope: number of events = " << tree3->GetEntries() << endl;

        if (tree1->GetEntries() == 0 || tree2->GetEntries() == 0 || tree3->GetEntries() == 0)
        {
            cerr << "ROOT file does not exist!!!" << endl;
            return;
        }

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
            #ifdef DEBUG
                if (totalEvt % 10000 == 0)
                {
                    cout << "event: " << totalEvt << endl;
                }
            #endif
            totalEvt++;

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
            hodoHit = { };
            maxPEHodo = { };
            maxChHodo = { };

            maxChProtoXOfXY = 0;
            maxChProtoYOfXY = 0;
            maxChProtoXOfXZ = 0;
            maxPEProtoXY = 0.0;
            maxPEProtoXOfXZ = 0.0;


            goodEvent = false;
            singleHit = false;
            // scintiHit = false;
            scintiHitEachSurface = { };
            isStraightBeam = false;
            goodEventForCT = false;
            goodEventForCTUpDown = false;
            goodEventForCTDownOnly = false;
            goodEventForCTStraight = false;
            goodEventForCTCell = false;
            goodEventForCTDiagOpp = {};

            // Hodoscope loop
            for (int ch = 0; ch < NScifi; ch++)
            {
                pe[static_cast<int> (EEasiroc::Hodoscope)][ch] = (double(adc[static_cast<int> (EEasiroc::Hodoscope)][ch]) - mean_ped[2][ch]) / gain[2][ch];
                // pe[static_cast <int>(EEasiroc::Hodoscope)][ch] = double(floor(pe[static_cast <int>(EEasiroc::Hodoscope)][ch] * 10) / 10);

                tuple<EHodoscope, int> hodoCh = ToHodoscopeCh(ch, shiftHSX1, shiftHSY1, shiftHSY2, shiftHSX2);

                EHodoscope hodoPlace = get<0> (hodoCh);
                if (hodoPlace == EHodoscope::None)
                    continue;

                int hodoNumber = get<1> (hodoCh);

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
                    // allHitCountHodo[static_cast<int> (hodoPlace)]++;
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

            // Hit multiplicity
            for(int i=0; i<NHodo; ++i)
            {
                hHitMultiplicity[i]->Fill(hitCountHodo[i]);
            }


            // 各種Good eventの条件
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
            if (singleHit && maxChHodo[static_cast<int> (EHodoscope::HSX1)] == maxChHodo[static_cast<int> (EHodoscope::HSX2)] && maxChHodo[static_cast<int> (EHodoscope::HSY1)] == maxChHodo[static_cast<int> (EHodoscope::HSY2)])
            {
                hHodoHitMapWithStraightBeam->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
                isStraightBeam = true;
                ++countStraight;
                countHodoHitEachCell[maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]++;
            }


            if (goodEvent && singleHit)
            {
                hHodoHitMapUp->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX1)], maxChHodo[static_cast<int> (EHodoscope::HSY1)]);
                hHodoHitMapDown->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
            }


            // Cross talkを出す用のGood eventフラグ（3つの条件があるが、そのうち1つを下で選択している）
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
            // goodEventForCT = goodEventForCTUpDown;

            // 対角キューブクロストーク用
            // 左上のキューブにヒット
            if (isStraightBeam && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= 4 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= 13)
            {
                goodEventForCTDiagOpp[0] = true;
            }
            // 右上
            else if (isStraightBeam && maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= 13 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= 13)
            {
                goodEventForCTDiagOpp[2] = true;
            }
            // 左下
            else if (isStraightBeam && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= 4 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= 4)
            {
                goodEventForCTDiagOpp[3] = true;
            }
            // 右下
            else if (isStraightBeam && maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= 13 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= 4)
            {
                goodEventForCTDiagOpp[1] = true;
            }

            // CellごとのCross talk分布出す用
            if (isStraightBeam &&
                maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXForCTCell &&
                maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXForCTCell &&
                maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYForCTCell &&
                maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYForCTCell // &&
                // maxChProtoXOfXY == CenterCubeXOfXY &&
                // maxChProtoYOfXY == CenterCubeYOfXY
                )
            {
                goodEventForCTCell = true;
                countCTCell++;
            }

            // Good event ここまで


            // 使用していないチャンネルの光量も出力
            if(isStraightBeam)
            {
                for(auto ch : usingProtoChForAllOutput)
                {
                    EEasiroc easiroc = get<0> (ch);
                    int easirocCh = get<1> (ch);
                    EScintiSurface surface = get<2> (ch);
                    int horizontal = get<3> (ch);
                    int vertical = get<4> (ch);

                    pe[static_cast<int> (easiroc)][easirocCh] = (double(adc[static_cast<int> (easiroc)][easirocCh]) - mean_ped[static_cast<int> (easiroc)][easirocCh]) / gain[static_cast<int> (easiroc)][easirocCh];
                    // cout << pe[static_cast<int> (easiroc)][easirocCh] << endl;

                    hPEProtoEach[static_cast<int> (surface)][horizontal-1][vertical-1]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
                }
            }


            // prototype loop
            for (auto itr = usingProtoCh.begin(); itr != usingProtoCh.end(); ++itr)
            {
                EEasiroc easiroc = get<0> (*itr);
                int easirocCh = get<1> (*itr);
                EScintiSurface surface = get<2> (*itr);
                int horizontal = get<3> (*itr);
                int vertical = get<4> (*itr);
                #ifdef DEBUG
                    cout << "usingProtoCh: " << easirocCh << endl;
                #endif
                // pe[static_cast<int> (easiroc)][easirocCh] = (double(adc[static_cast<int> (easiroc)][easirocCh]) - mean_ped[static_cast<int> (easiroc)][easirocCh]) / gain[static_cast<int> (easiroc)][easirocCh];
                // pe[static_cast <int>(easiroc)][easirocCh] = double(floor(pe[static_cast <int>(easiroc)][easirocCh] * 10) / 10);

                hProto[static_cast<int> (surface)]->SetBinContent(horizontal, vertical, pe[static_cast<int> (easiroc)][easirocCh]);


                if (pe[static_cast<int> (easiroc)][easirocCh] > ProtoPEThreshold)
                {
                    hHitRateProto[static_cast<int> (surface)]->Fill(horizontal, vertical);
                    hPEProto[static_cast<int> (surface)]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
                    // hPEProtoEach[static_cast<int> (surface)][horizontal][vertical]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);

                    // scintiHit = true;

                    if (goodEvent && singleHit)
                    {
                        // hPEProtoGood[static_cast<int> (surface)]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
                    }

                    // Search a channel of max pe
                    if (surface == EScintiSurface::XY && pe[static_cast<int> (easiroc)][easirocCh] > maxPEProtoXY)
                    {
                        maxChProtoXOfXY = horizontal;
                        maxChProtoYOfXY = vertical;
                        maxPEProtoXY = pe[static_cast<int> (easiroc)][easirocCh];
                    }
                    if (surface == EScintiSurface::XZ && pe[static_cast<int> (easiroc)][easirocCh] > maxPEProtoXOfXZ)
                    {
                        maxChProtoXOfXZ = horizontal;
                        maxPEProtoXOfXZ = pe[static_cast<int> (easiroc)][easirocCh];
                    }

                    if (isStraightBeam && horizontal == 3 && vertical == 3)
                    {
                        scintiHitEachSurface[static_cast<int> (surface)] = true;
                    }
                }
            }


            // プロトタイプの読み出し面ごとのgood event数
            for (int i = 0; i < NSurfaceScinti; i++)
            {
                if (scintiHitEachSurface[i])
                    countScintiHitEachSurfaceCell[i][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]++;
            }


            centerPEXY = pe[static_cast<int> (EEasiroc::Scinti1)][9];
            // centerPEXZ = pe[static_cast<int> (EEasiroc::Scinti2)][41];
            // leftPEXZ = pe[static_cast<int> (EEasiroc::Scinti2)][40];
            for (int i = 0; i < NCubeCT; i++)
            {
                aroundPEXY[i] = pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXY[i]];
            }





            // ホドスコープ内側4ch×4chヒットをクロストークのヒストグラムに使う
            if (goodEventForCT)
            {
                for (int i = 0; i < NCubeCT; i++)
                {
                    hCrosstalkXY[i]->Fill(aroundPEXY[i] / centerPEXY);
                    hPEAroundForCTXY[i]->Fill(aroundPEXY[i]);
                    scatterCTXY[i]->SetPoint(countCTPoint, centerPEXY, aroundPEXY[i]);
                    hCrosstalkScatterXY[i]->Fill(centerPEXY, aroundPEXY[i]);
                    #ifdef DEBUG
                        cout << "aroundPEXY" << i << " = " << aroundPEXY[i] << endl;
                    #endif

                    if (aroundPEXY[i] < DarkCutPEForCT)
                    {
                        hCrosstalkXYDarkCut[i]->Fill(0);
                    }
                    else
                    {
                        hCrosstalkXYDarkCut[i]->Fill(aroundPEXY[i] / centerPEXY);
                    }
                }

                hPECenterForCTXY->Fill(centerPEXY);

                countCTPoint++;
            }

            // 対角キューブのクロストーク
            for (int i=0; i<NCubeCTDiagOpp; i++)
            {
                if(goodEventForCTDiagOpp[i])
                {
                    double noHitPEXY = pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][1]];
                    double hitPEXY = pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][0]];
                    hPENoHitForCTXYDiagOpp[i]->Fill(noHitPEXY);
                    hPEHitForCTXYDiagOpp[i]->Fill(hitPEXY);
                    hCrosstalkScatterXYDiagOpp[i]->Fill(hitPEXY, noHitPEXY);
                    if (noHitPEXY < DarkCutPEForCT)
                    {
                        hCrosstalkXYDarkCutDiagOpp[i]->Fill(0);
                    }
                    else
                    {
                        hCrosstalkXYDarkCutDiagOpp[i]->Fill(noHitPEXY / hitPEXY);
                    }
                    break;
                }
            }


            // セルごとのクロストーク
            if (goodEventForCTCell)
            {
                // #ifdef DEBUG
                //
                //     if (leftPEXZ / centerPEXZ > 0.1 && leftPEXZ < centerPEXZ && leftPEXZ >= DarkCutPEForCT && maxChHodo[static_cast<int> (EHodoscope::HSX2)] == 9 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= 5)
                //     {
                //         cout << "selected Hodo Ch: HSX2=" << maxChHodo[static_cast<int> (EHodoscope::HSX2)] << " (" << maxPEHodo[static_cast<int> (EHodoscope::HSX2)] <<  " p.e.), HSY2=" << maxChHodo[static_cast<int> (EHodoscope::HSY2)] << " (" << maxPEHodo[static_cast<int> (EHodoscope::HSY2)] <<  " p.e.)" << endl;
                //         cout << "evt=" << totalEvt << ", XZ, left: " << leftPEXZ << ", center: " << centerPEXZ << endl;
                //     }
                //     if (leftPEXY / centerPEXY > 0.1 && leftPEXY < centerPEXY && leftPEXY >= DarkCutPEForCT && maxChHodo[static_cast<int> (EHodoscope::HSX2)] == 9 && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= 5)
                //     {
                //         cout << "selected Hodo Ch: HSX2=" << maxChHodo[static_cast<int> (EHodoscope::HSX2)] << " (" << maxPEHodo[static_cast<int> (EHodoscope::HSX2)] <<  " p.e.), HSY2=" << maxChHodo[static_cast<int> (EHodoscope::HSY2)] << " (" << maxPEHodo[static_cast<int> (EHodoscope::HSY2)] <<  " p.e.)" << endl;
                //         cout << "evt=" << totalEvt << ", XY, left: " << leftPEXY << ", center: " << centerPEXY << endl;
                //     }
                // #endif

                // Center cubeの光量が一番多いとき
                if (maxChProtoXOfXY == 3 && maxChProtoYOfXY == 3)
                {
                    for (int i = 0; i < NCubeCT; i++)
                    {
                        // 周りのキューブの光量が一定（0.5pe）未満ならば0をFillする
                        if (aroundPEXY[i] < DarkCutPEForCT)
                        {
                            hCrosstalkXYDarkCutEachCell[i][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(0);
                        }
                        else     // if(leftPEXY*RatioCutForCT < centerPEXY)
                        {
                            hCrosstalkXYDarkCutEachCell[i][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(aroundPEXY[i] / centerPEXY);
                        }
                        hCrosstalkScatterXYEachCell[i][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(centerPEXY, aroundPEXY[i]);
                    }
                }

                // 四隅のいずれかのキューブの光量が最大
                if ((maxChProtoXOfXY == 2 || maxChProtoXOfXY == 4) && (maxChProtoYOfXY == 2 || maxChProtoYOfXY == 4))
                {
                    int i;
                    // 左上のキューブの光量が最大
                    if (maxChProtoXOfXY == 2 && maxChProtoYOfXY == 4)
                    {
                        i = 0;
                    }
                    // 右下
                    else if (maxChProtoXOfXY == 4 && maxChProtoYOfXY == 2)
                    {
                        i = 1;
                    }
                    // それ以外ならスキップ
                    else
                    {
                        i = -1;
                    }

                    if (i >= 0)
                    {
                        if (pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][1]] < DarkCutPEForCT)
                        {
                            hCrosstalkXYDarkCutEachCellDiagOpp[0][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(0);
                        }
                        else
                        {
                            hCrosstalkXYDarkCutEachCellDiagOpp[0][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][1]] / pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][0]]);
                        }
                        // hCrosstalkScatterXYEachCellDiagOpp[0][maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][0]], pe[static_cast<int> (EEasiroc::Scinti1)][CubeChMapXYDiagOpp[i][1]]);
                    }
                }
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
            if (totalEvt >= min_evt && totalEvt <= max_evt && goodEventForCTCell)
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
                figName = TString::Format("%shit_%04d_%04d_evt%d.%s", EvtDisplayDir.c_str(), runnum, subrun, totalEvt, outputFileType.c_str());
                canvas->SaveAs(figName);
                canvas->Clear();
            }
        }
    }


    // Draw Histograms
    const double RightMarginForHodoMap = 0.15;

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

    figName = TString::Format("%s/HitRate_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
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
    for (auto itr = usingProtoChForAllOutput.begin(); itr != usingProtoChForAllOutput.end(); ++itr)
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
        // gPad->SetLogy();
        hPEProtoEach[static_cast<int> (surface)][horizontal-1][vertical-1]->Draw();

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


    #ifdef DEBUG
        cout << endl << "--------------event loop end--------------" << endl << "total event: " << totalEvt << endl;
    #endif


    // Hit multiplicity
    nHistHori = 2;
    nHistVert = 2;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    for(int i=0; i<NHodo; ++i)
    {
        canvas->cd(i+1);
        hHitMultiplicity[i]->Draw();
    }
    figName = TString::Format("%sHitMultiplicity_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();




    // Crosstalk
    // normal
    // nHistHori = 2;
    // nHistVert = 1;
    TF1* poissonFix = new TF1("poissonFix", "[0]*TMath::Poisson(x*[1], [1]*[2])", FitRangeCT[0], FitRangeCT[1]);
    for (int i = 0; i < NCubeCT; i++)
    {
        canvas = new TCanvas();
        hCrosstalkXY[i]->Draw();
        // poissonFix->SetParameters(300, 100, 0.04);
        // hCrosstalkXY->Fit("poissonFix", FitOption, "", FitRangeCT[0], FitRangeCT[1]);
        gStyle->SetOptStat(2210);
        // gPad->SetLogy();
        // gStyle->SetOptFit(111);
        changestatsBoxSize(hCrosstalkXY[i], 0.7, 0.99, 0.65, 0.935);
        figName = TString::Format("%sCrosstalk%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }


    // crosstalk, dark count cut
    for (int i = 0; i < NCubeCT; i++)
    {
        canvas = new TCanvas();
        hCrosstalkXYDarkCut[i]->Draw();
        // poissonFix->SetParameters(300, 100, 0.01);
        // hCrosstalkXYDarkCut->Fit("poissonFix", FitOption, "", FitRangeCTDarkCut[0], FitRangeCTDarkCut[1]);
        gStyle->SetOptStat(2210);
        // gPad->SetLogy();
        // gStyle->SetOptFit(111);
        changestatsBoxSize(hCrosstalkXYDarkCut[i], 0.7, 0.99, 0.65, 0.935);
        figName = TString::Format("%sCrosstalkDarkCut%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SetLogy();
        canvas->SaveAs(figName);
        canvas->Clear();
    }


    // center L.Y. for crosstalk
    TFile histsCenter(TString::Format("%sPECenter_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    canvas = new TCanvas();
    // hPECenterForCTXY->Fit("landau", FitOption, "", FitRangePECenterXY[0], FitRangePECenterXY[1]);
    hPECenterForCTXY->Draw();
    hPECenterForCTXY->Write();
    histsCenter.Close();

    figName = TString::Format("%sPECenter_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    // around L.Y. for crosstalk
    TFile histsAround(TString::Format("%sPEAround_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    for (int i = 0; i < NCubeCT; i++)
    {
        canvas = new TCanvas();
        // hPEAroundForCTXY[i]->Fit();
        hPEAroundForCTXY[i]->Draw();
        hPEAroundForCTXY[i]->Write();

        figName = TString::Format("%sPEAround%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }
    histsAround.Close();


    // // left only
    // TF1* poissonI = new TF1("poissonI", "[0]*TMath::PoissonI(x+0.5, [1])", 0, 10);
    // TFile histsLeft(TString::Format("%sPELeft_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    // canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    // canvas->Divide(nHistHori, nHistVert);
    // canvas->cd(1);
    // poissonI->SetParameters(3000, 0.5);
    // hPELeftForCTXY->Fit("poissonI", FitOption, "", FitRangePELeftXY[0], FitRangePELeftXY[1]);
    // hPELeftForCTXY->Draw();
    // hPELeftForCTXY->Write();
    //
    // canvas->cd(2);
    // poissonI->SetParameters(3000, 0.5);
    // hPELeftForCTXZ->Fit("poissonI", FitOption, "", FitRangePELeftXZ[0], FitRangePELeftXZ[1]);
    // hPELeftForCTXZ->Draw();
    // hPELeftForCTXZ->Write();
    // histsLeft.Close();
    //
    // figName = TString::Format("%sPELeft_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    // canvas->SaveAs(figName);
    // canvas->Clear();

    // Each Cell
    TFile histsCTXYDarkCutEachCell(TString::Format("%sCrosstalkXYDarkCutEachCell_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    canvas = new TCanvas();
    for (int i = 0; i < NCubeCT; i++)
    {
        for (int j = 0; j < NScifiEachHodo; j++)
        {
            for (int k = 0; k < NScifiEachHodo; k++)
            {
                hCrosstalkXYDarkCutEachCell[i][j][k]->Draw();
                hCrosstalkXYDarkCutEachCell[i][j][k]->Write();
                hCrosstalkXYDarkCutMap[i]->SetBinContent(j + 1, k + 1, hCrosstalkXYDarkCutEachCell[i][j][k]->GetMean() * 100);     // *100 means percentile
            }
        }
    }
    histsCTXYDarkCutEachCell.Close();
    canvas->Clear();




    cout << "Straight: " << countStraight << endl;
    cout << "StraightForCT: " << countCTStraight << endl;
    // cout << "UpDown: " << countCTUpDown << endl;
    // cout << "DownOnly: " << countCTDownOnly << endl;
    // cout << "EachCell: " << countCTCell << endl;


    // scatter plot of cross talk (TGraph)
    for (int i = 0; i < NCubeCT; i++)
    {
        canvas = new TCanvas();
        // canvas->Divide(nHistHori, nHistVert);
        // canvas->cd(1);
        gStyle->SetOptFit(111);
        scatterCTXY[i]->Fit("pol1", FitOption, "", MinPEScatterCTCenter, MaxPEScatterCTCenter);
        scatterCTXY[i]->SetTitle(TString::Format("L.Y. %s vs center (using Z readout)", CubeGeometryTitle[i].c_str()));
        scatterCTXY[i]->GetYaxis()->SetTitle(TString::Format("L.Y. %s (ch%d) (p.e.)", CubeGeometryTitle[i].c_str(), CubeChMapXY[i]));
        scatterCTXY[i]->GetXaxis()->SetTitle("L.Y. center (ch8) (p.e.)");
        scatterCTXY[i]->GetXaxis()->SetRangeUser(MinPEScatterCTCenter, MaxPEScatterCTCenter);
        scatterCTXY[i]->GetYaxis()->SetRangeUser(MinPEScatterCTAround, MaxPEScatterCTAround);
        scatterCTXY[i]->Draw("AP");
        // canvas->cd(2);
        // scatterCTXZ->SetTitle("L.Y. left vs center (using Y readout)");
        // scatterCTXZ->GetYaxis()->SetTitle("L.Y. left (ch40) (p.e.)");
        // scatterCTXZ->GetXaxis()->SetTitle("L.Y. center (ch41) (p.e.)");
        // scatterCTXZ->Fit("pol1", FitOption, "", MinPEScatterCTCenter, MaxPEScatterCTCenter);
        // scatterCTXZ->GetXaxis()->SetRangeUser(MinPEScatterCTCenter, MaxPEScatterCTCenter);
        // scatterCTXZ->GetYaxis()->SetRangeUser(MinPEScatterCTAround, MaxPEScatterCTAround);
        // scatterCTXZ->Draw("AP");
        //
        figName = TString::Format("%sCrosstalkScatterPlot%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }

    // scatter plot of crosstalk (TH2D)
    canvas = new TCanvas();
    for (int i = 0; i < NCubeCT; i++)
    {
        hCrosstalkScatterXY[i]->Draw("colz");

        gPad->SetRightMargin(RightMarginForHodoMap);
        changestatsBoxSize(hCrosstalkScatterXY[i], 0.65, 0.88, 0.6, 0.9);
        figName = TString::Format("%sCrosstalkScatterHist%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
    }
    canvas->Clear();

    // ここから重いのでいるときだけコメントアウト外す

    // TFile histsCTScatterXYEachCell(TString::Format("%sCrosstalkScatterHistXYEachCell_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    // canvas = new TCanvas();
    //
    // for (int i = 0; i < NCubeCT; i++)
    // {
    //     for (int j = 0; j < NScifiEachHodo; j++)
    //     {
    //         for (int k = 0; k < NScifiEachHodo; k++)
    //         {
    //             hCrosstalkScatterXYEachCell[i][j][k]->Draw("colz");
    //             gPad->SetRightMargin(RightMarginForHodoMap);
    //             changestatsBoxSize(hCrosstalkScatterXYEachCell[i][j][k], 0.65, 0.88, 0.6, 0.9);
    //             hCrosstalkScatterXYEachCell[i][j][k]->Write();
    //             figName = TString::Format("%sCrosstalkScatterHist%dX%dY%d_%04d_%04d.%s", CrosstalkDir.c_str(), i, j + 1, k + 1, runnum, subrun, outputFileType.c_str());
    //             canvas->SaveAs(figName);
    //         }
    //     }
    // }
    //
    // histsCTScatterXYEachCell.Close();
    // canvas->Clear();

    // ここまで

    // 対角キューブのクロストーク
    const string CubeGeometryFileNameDiagOpp[] = { "LROverUL", "ULOverLR", "LLOverUR", "UROverLL" };
    const string FileNamePEDiagOpp[] = { {}, "ULOverLR", "LLOverUR", "UROverLL" };
    for(int i=0; i<NCubeCTDiagOpp; ++i)
    {
        figName = TString::Format("%sCrosstalkDiagOpp%s_%04d_%04d.%s", ResultDir.c_str(), CubeGeometryFileNameDiagOpp[i].c_str(), runnum, subrun, outputFileType.c_str());
        SaveHist(hCrosstalkXYDarkCutDiagOpp[i], figName, "", true);

        figName = TString::Format("%sCrosstalkScatterHistDiagOpp%s_%04d_%04d.%s", ResultDir.c_str(), CubeGeometryFileNameDiagOpp[i].c_str(), runnum, subrun, outputFileType.c_str());
        SaveHist(hCrosstalkScatterXYDiagOpp[i], figName, "colz");

        figName = TString::Format("%sPEDiagOppHit%d%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        SaveHist(hPEHitForCTXYDiagOpp[i], figName);

        figName = TString::Format("%sPEDiagOppNoHit%d%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        SaveHist(hPENoHitForCTXYDiagOpp[i], figName);
    }





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
    hGapCount1s->SetStats(kFALSE);
    hGapCount2s->SetStats(kFALSE);

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
    histWidth = 1280;
    histHeight = 1200;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);

    canvas->cd(1);
    // hHodoHitMapUp->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapUp->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapUp->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapUp->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapUp->SetStats(kFALSE);
    // hHodoHitMapUp->SetMarkerSize(2.0);
    hHodoHitMapUp->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);

    canvas->cd(2);
    // hHodoHitMapDown->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapDown->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapDown->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapDown->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapDown->SetStats(kFALSE);
    // hHodoHitMapDown->SetMarkerSize(2.0);
    hHodoHitMapDown->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);

    canvas->cd(3);
    hHodoHitMapWithStraightBeam->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->SetStats(kFALSE);
    hHodoHitMapWithStraightBeam->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);

    figName = TString::Format("%sHodoHitMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    // Detection efficiency
    gStyle->SetPaintTextFormat("3.2f");
    for (int hodoX = 0; hodoX < NScifiEachHodo; hodoX++)
    {
        for (int hodoY = 0; hodoY < NScifiEachHodo; hodoY++)
        {
            if (countHodoHitEachCell[hodoX][hodoY] == 0)
            {
                // hDetectionEff[i]->SetBinContent(hodoX+1, hodoY+1, 0.0);
                continue;
            }
            for (int surface = 0; surface < NSurfaceScinti; surface++)
            {
                hDetectionEff[surface]->SetBinContent(hodoX + 1, hodoY + 1, (double) countScintiHitEachSurfaceCell[surface][hodoX][hodoY] / (double) countHodoHitEachCell[hodoX][hodoY]);
            }
        }
    }
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);
    for (int i = 0; i < NSurfaceScinti; i++)
    {
        canvas->cd(i + 1);
        hDetectionEff[i]->GetXaxis()->SetNdivisions(NScifiEachHodo);
        hDetectionEff[i]->GetYaxis()->SetNdivisions(NScifiEachHodo);
        hDetectionEff[i]->SetStats(kFALSE);
        gStyle->SetPadGridX(1);
        gStyle->SetPadGridY(1);
        hDetectionEff[i]->Draw("text colz");
        drawCubeLine("inj", 2);
        gPad->SetRightMargin(RightMarginForHodoMap);
    }
    figName = TString::Format("%sDetectionEfficiency_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    const Int_t NRGBs = 3;
    const Int_t NCont = 255;

    double stops[NRGBs] = { 0.00, .50, 1.00 };
    double red[NRGBs] = { 1.00, 1.0, 1.00 };
    double green[NRGBs] = { 1.00, 0.0, 0.00 };
    double blue[NRGBs] = { 1.00, 0.0, 0.00 };
    TColor::CreateGradientColorTable(2, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Beam hit position dependency of cross talk (Hodomap)
    nHistHori = 1;
    nHistVert = 1;
    histWidth = 1280;
    histHeight = 1200;
    gStyle->SetPaintTextFormat("3.2f");
    for (int i = 0; i < NCubeCT; i++)
    {
        canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
        hCrosstalkXYDarkCutMap[i]->Draw("text colz");
        hCrosstalkXYDarkCutMap[i]->GetXaxis()->SetNdivisions(NScifiEachHodo);
        hCrosstalkXYDarkCutMap[i]->GetYaxis()->SetNdivisions(NScifiEachHodo);
        gPad->SetRightMargin(RightMarginForHodoMap);
        hCrosstalkXYDarkCutMap[i]->SetStats(kFALSE);

        drawCubeLine("inj");

        figName = TString::Format("%sCrosstalkDarkCutMap%d_%04d_%04d.%s", ResultDir.c_str(), i, runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        canvas->Clear();
    }
}

int main(int argc, char** argv)
{
    if (argc <= 2)
    {
        cerr << "run_proto <runnum> <# of files> <shiftHSX1> <shiftHSY1> <shiftHSY2> <shiftHSX2> <min_evt> <max_evt> <fileType>" << endl;

        return -1;
    }

    int runnum = atoi(argv[1]);
    int fileCount = atoi(argv[2]);
    int x1, y1, y2, x2;
    if (argc >= 3 + 4)
    {
        x1 = atoi(argv[3]);
        y1 = atoi(argv[4]);
        y2 = atoi(argv[5]);
        x2 = atoi(argv[6]);
    }
    int min_evt, max_evt;
    if (argc >= 3 + 4 + 2)
    {
        min_evt = atoi(argv[7]);
        max_evt = atoi(argv[8]);
    }

    string fileType;
    if (argc >= 3 + 4 + 2 + 1)
    {
        fileType = argv[9];
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
    else if (argc == 3 + 4 + 2 + 1)
    {
        run_proto(runnum, fileCount, x1, y1, y2, x2, min_evt, max_evt, fileType);
    }


    return 0;
}
