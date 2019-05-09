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
#include <TLine.h>
#include <TEllipse.h>
#include <TGraphErrors.h>

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

bool isGap(int hodoX, int hodoY, bool isHitToCube, E3DOption option, int shiftX=0, int shiftY=0)
{
    // 1cm cube & 2cm cube with parallel placing
    if (isHitToCube && (option == E3DOption::WRef1cmP || option == E3DOption::WRef1cmV || option == E3DOption::WORef1cmP || option == E3DOption::WORef1cmV || option == E3DOption::WRef2cmP) && !(hodoX + shiftX >= 6 && hodoX + shiftX <= 11 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }
    // 2cm cube with vertical placing
    if (isHitToCube && option == E3DOption::WRef2cmV && !(hodoX + shiftX >= 3 && hodoX + shiftX <= 14 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }
    if (!isHitToCube && (option == E3DOption::WRef1cmP || option == E3DOption::WRef1cmV || option == E3DOption::WORef1cmP || option == E3DOption::WORef1cmV || option == E3DOption::WRef2cmP) && (hodoX + shiftX >= 6 && hodoX + shiftX <= 11 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }
    if (!isHitToCube && option == E3DOption::WRef2cmV && (hodoX + shiftX >= 3 && hodoX + shiftX <= 14 && hodoY + shiftY >= 6 && hodoY + shiftY <= 11))
    {
        return true;
    }
    return false;
}

void changeStatsBoxSize(TH1* hist, double x1, double x2, double y1, double y2)
{
    gPad->Update();
    TPaveStats* st = (TPaveStats *) hist->FindObject("stats");
    st->SetX1NDC(x1);
    st->SetX2NDC(x2);
    st->SetY1NDC(y1);
    st->SetY2NDC(y2);
}

void changeStatsBoxSize(TH1* hist, array<double, 4> size)
{
    gPad->Update();
    TPaveStats* st = (TPaveStats *) hist->FindObject("stats");
    st->SetX1NDC(size[0]);
    st->SetX2NDC(size[1]);
    st->SetY1NDC(size[2]);
    st->SetY2NDC(size[3]);
}

// void changeStatsBoxSize(TH2* hist, double x1, double x2, double y1, double y2)
// {
//     gPad->Update();
//     TPaveStats* st = (TPaveStats *) hist->FindObject("stats");
//     st->SetX1NDC(x1);
//     st->SetX2NDC(x2);
//     st->SetY1NDC(y1);
//     st->SetY2NDC(y2);
// }

void drawCubeLine(string config)
{
    // TODO: E3DOptionへの対応
    // 穴は一つしかないから

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

    const int LineColor = 7;
    // 0:white, 1:black, 2:red, 3:green, 4:blue, 5:yellow, 6:magenta, 7:cyan, 8:dark green, 9:purple
    const int LineStyle = 2;
    // 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    const int LineWidth = 8;

    double xShift = 0;
    double yShift = 0;


    if (config == "ex1")
    {
        xShift = -0.2;
        yShift =  0.1;
    }
    else if (config == "ex2")
    {
        xShift = -0.3;
        yShift =  0.3;
    }
    else if (config == "in1")
    {
        xShift =  0.;
        yShift =  0.3;
    }
    else if (config == "ex_w1H")
    {
        xShift = -0.2;
        yShift =  0.2;
    }
    else if (config == "ex_w2H")
    {
        xShift = -0.7;
        yShift =  0.35;
    }
    else if (config == "ex_w1V")
    {
        xShift = -0.2;
        yShift = 0.1;
    }
    else if (config == "ex_w2V")
    {
        xShift = -0.6;
        yShift =  0.45;
    }


    TLine* line1 = new TLine(low_x + xShift, low_y + yShift, low_x + xShift, upp_y + yShift);
    line1->SetLineColor(LineColor);
    line1->SetLineWidth(LineWidth);
    line1->SetLineStyle(LineStyle);
    TLine* line2 = new TLine(upp_x + xShift, low_y + yShift, upp_x + xShift, upp_y + yShift);
    line2->SetLineColor(LineColor);
    line2->SetLineWidth(LineWidth);
    line2->SetLineStyle(LineStyle);
    TLine* line3 = new TLine(low_x + xShift, upp_y + yShift, upp_x + xShift, upp_y + yShift);
    line3->SetLineColor(LineColor);
    line3->SetLineWidth(LineWidth);
    line3->SetLineStyle(LineStyle);
    TLine* line4 = new TLine(low_x + xShift, low_y + yShift, upp_x + xShift, low_y + yShift);
    line4->SetLineColor(LineColor);
    line4->SetLineWidth(LineWidth);
    line4->SetLineStyle(LineStyle);
    TLine* line5 = new TLine(low_fiber_pos + xShift, low_y + yShift, low_fiber_pos + xShift, upp_y + yShift);
    line5->SetLineColor(LineColor);
    line5->SetLineWidth(LineWidth);
    line5->SetLineStyle(LineStyle);
    TLine* line6 = new TLine(upp_fiber_pos + xShift, low_y + yShift, upp_fiber_pos + xShift, upp_y + yShift);
    line6->SetLineColor(LineColor);
    line6->SetLineWidth(LineWidth);
    line6->SetLineStyle(LineStyle);
    TLine* line7 = new TLine(low_x + xShift, low_fiber_pos + yShift, upp_x + xShift, low_fiber_pos + yShift);
    line7->SetLineColor(LineColor);
    line7->SetLineWidth(LineWidth);
    line7->SetLineStyle(LineStyle);
    TLine* line8 = new TLine(low_x + xShift, upp_fiber_pos + yShift, upp_x + xShift, upp_fiber_pos + yShift);
    line8->SetLineColor(LineColor);
    line8->SetLineWidth(LineWidth);
    line8->SetLineStyle(LineStyle);
    TEllipse* circle = new TEllipse(low_x + center_pos + xShift, low_y + center_pos + yShift, Zradius, Zradius);
    circle->SetLineColor(LineColor);
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

// キューブ種類・配置のオプションとrunnumの対応
E3DOption decideOption(int runnum)
{
    if (runnum == 0 || (runnum >= 4 && runnum <= 7) || (runnum >= 9 && runnum <= 13))
        return E3DOption::WRef1cmV; // white refector
    if (runnum >= 14 && runnum <= 23)
        return E3DOption::WRef2cmV;
    if (runnum == 24 || (runnum >= 28 && runnum <= 36))
        return E3DOption::WRef1cmP;
    if ((runnum >= 37 && runnum <= 43) || (runnum >= 49 && runnum <= 51))
        return E3DOption::WRef2cmP;
    if (runnum >= 56 && runnum <= 58 || runnum >= 66 && runnum <= 68)
        return E3DOption::WRef1cmV; // silver ref.
    if (runnum >= 69 && runnum <= 74)
        return E3DOption::WRef1cmP;
    if ((runnum >= 76 && runnum <= 77) || (runnum >= 81 && runnum <= 84))
        return E3DOption::WORef1cmV;
    if (runnum >= 85 && runnum <= 91)
        return E3DOption::WORef1cmP;
    cerr << "Invalid run number!!!" << endl;
    return E3DOption::None;
}

// E3DOptionに応じて、ホドスコープに射影したときに
// 1cmサイズ: true
// 2cmサイズ: false
// 存在しないrunnumはあらかじめ弾かれていることが前提（弾かれていないと存在しないrunnumに対してtrueを返す）
bool is1cmArea(E3DOption option)
{
    if (option == E3DOption::WRef2cmV)
        return false;
    return true;
}

bool isVertical(E3DOption option)
{
    if (option == E3DOption::WRef1cmV || option == E3DOption::WRef2cmV || option == E3DOption::WORef1cmV)
        return true;
    return false;
}

void run_proto(int runnum, int fileCount, int shiftHSX1=0, int shiftHSY1=0, int shiftHSY2=0, int shiftHSX2=0, int min_evt=-1, int max_evt=-1, string outputFileType="png", int gap_point_0=0, int gap_pt2s_0=0, int gap_pt1s_0=0, int gap_hs_0=0, int gap_point_1=0, int gap_pt2s_1=0, int gap_pt1s_1=0, int gap_hs_1=0)
{
    // constants
    gErrorIgnoreLevel = kError;

    const double HodoMapZoomedSize1cm = 10;    // 中央 n ch x n ch
    const double HodoMapZoomedSize2cm = 10;
    const double HodoMapCenter = 8.5;


    const double RightMarginForHodoMap = 0.15;

    int fitRangeLower = 5;
    int fitRangeHigher = 5;
    double minPEEachCellFit = 2.5;
    double maxPEEachCellFit = 18.5;


    const double CutPEForEachCell = 2.0;

    const int hodoLowXFor1cm = 7;
    const int hodoHighXFor1cm = 10;
    const int hodoLowYFor1cm = 7;
    const int hodoHighYFor1cm = 10;

    const int hodoLowXFor2cm = 4;
    const int hodoHighXFor2cm = 13;
    const int hodoLowYFor2cm = 7;
    const int hodoHighYFor2cm = 10;


    int omit_lowx1 = 0;
    int omit_highx1 = 0;
    int omit_lowy1 = 0;
    int omit_highy1 = 0;

    int omit_lowx2 = 0;
    int omit_highx2 = 0;
    int omit_lowy2 = 0;
    int omit_highy2 = 0;

    int subrun = 6;

    const double HodoPEThreshold = 2.5;
    const double ProtoPEThreshold = 2.5;
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
    string resultDirTemp = "3DPrint/";
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

    mkdir(ResultDir.c_str(), 0777);
    mkdir(ScintiPEDir.c_str(), 0777);
    mkdir(EvtDisplayDir.c_str(), 0777);
    mkdir(HodoPEDir.c_str(), 0777);
    mkdir(CellPEDir.c_str(), 0777);


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
    const double MinPEHodo = -0.5;
    const double MaxPEHodo = 69.5;
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

    TH1D* hPEProto;
    TH1D* hPEProtoGood;
    const double MinPEProto = -0.5;
    const double MaxPEProto = 99.5;
    const int NBinPEProto = 100;
    histName = "hPEZY";
    histAxis = "L.Y. with all events (X readout);Light yield (p.e.);Number of events";
    hPEProto = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);
    // hPEProto = new TH1D(histName.c_str(), histAxis.c_str(), 100, -0.5, 9.5);
    histName += "Good";
    histAxis = "L.Y. with good events (X readout);Light yield (p.e.);Number of events";
    hPEProtoGood = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);

    TH1D* hPEProtoGoodCell[NScifiEachHodo][NScifiEachHodo];
    for (int y = 0; y < NScifiEachHodo; y++)
    {
        for (int x = 0; x < NScifiEachHodo; x++)
        {
            histName = "hPECellZY_X" + to_string(x + 1) + "Y" + to_string(y + 1) + "";
            histAxis = "L.Y. ZY Cell (X,Y) = (" + to_string(x + 1) + "," + to_string(y + 1) + ");Light yield (p.e.);Number of events";
            hPEProtoGoodCell[x][y] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);
        }
    }

    TH2D* hPEProtoHodoMap = new TH2D("hPEProtoHodoMap", "L.Y. of scintillator (using landau MPV);cell # along X;cell # along Y; Light yield (p.e.)", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    hPEProtoHodoMap->SetMinimum(0);
    // hPEProtoHodoMap->SetMaximum(MaxPEProto+0.5);
    TH2D* hPEProtoHodoMapMean = new TH2D("hPEProtoHodoMapMean", "L.Y. of scintillator (using mean);cell # along X;cell # along Y; Light yield (p.e.)", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hPEProtoHodoMapMeanZoom = new TH2D("hPEProtoHodoMapMeanZoom", "L.Y. of scintillator (using mean);cell # along X;cell # along Y; Light yield (p.e.)", HodoMapZoomedSize1cm, HodoMapCenter-HodoMapZoomedSize1cm/2, HodoMapCenter+HodoMapZoomedSize1cm/2, HodoMapZoomedSize1cm, HodoMapCenter-HodoMapZoomedSize1cm/2, HodoMapCenter+HodoMapZoomedSize1cm/2);

    TH2D* hDetectionEff = new TH2D("hDetectionEff", "Detection efficiency;cell # along X;cell # along Y; Detection efficiency", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hDetectionEffZoom = new TH2D("hDetectionEffZoom", "Detection efficiency;cell # along X;cell # along Y; Detection efficiency", HodoMapZoomedSize1cm, HodoMapCenter-HodoMapZoomedSize1cm/2, HodoMapCenter+HodoMapZoomedSize1cm/2, HodoMapZoomedSize1cm, HodoMapCenter-HodoMapZoomedSize1cm/2, HodoMapCenter+HodoMapZoomedSize1cm/2);
    hDetectionEff->SetMinimum(0.0);
    hDetectionEff->SetMaximum(1.0);
    hDetectionEffZoom->SetMaximum(1.0);
    hDetectionEffZoom->SetMinimum(0.0);
    int countHodoHitEachCell[NScifiEachHodo][NScifiEachHodo] = { };
    int countScintiHodoHitEachCell[NScifiEachHodo][NScifiEachHodo] = { };

    TH2D* hHodoHitMapWithProtoHitUp = new TH2D("hHodoHitMapWithProtoHitUp", "Upstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapWithProtoHitDown = new TH2D("hHodoHitMapWithProtoHitDown", "Downstream hodoscope hitmap with scinti. hit;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapUp = new TH2D("hHodoHitMapUp", "Upstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);
    TH2D* hHodoHitMapDown = new TH2D("hHodoHitMapDown", "Downstream hodoscope hitmap;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    TH2D* hHodoHitMapWithStraightBeam = new TH2D("hHodoHitMapWithStraightBeam", "Hodoscope hitmap with straight beam event;cell # along X;cell # along Y;Number of events", NScifiEachHodo, MinHodoMap, MaxHodoMap, NScifiEachHodo, MinHodoMap, MaxHodoMap);

    // TH1D* hPEProtoEach[NSurfaceScinti][NChOneSide][NChOneSide];
    // vector<tuple<EEasiroc, int, EScintiSurface, int, int>> usingProtoCh;     // 使用チャンネルを格納 <EASIROC ch, プロトタイプ面, cube # horizontal, cube # vertical>
    // for (int i = 0; i < NEasiroc - 1; i++)
    // {
    //     EEasiroc easiroc;
    //     if (i == 0)
    //         easiroc = EEasiroc::Scinti2;
    //     else
    //         easiroc = EEasiroc::Scinti1;
    //
    //     ToScintiCh* toScintiCh = new ToScintiCh(EScintiType::Printed);
    //
    //
    //     for (int ch = 0; ch < NChEasiroc; ch++)
    //     {
    //         if (toScintiCh->isConnected(easiroc, ch))
    //         {
    //             tuple<EScintiSurface, int, int> scintiCh = toScintiCh->getCh(easiroc, ch);
    //
    //             EScintiSurface surface = get<0> (scintiCh);
    //             int horizontal = get<1> (scintiCh);
    //             int vertical = get<2> (scintiCh);
    //
    //             usingProtoCh.push_back(forward_as_tuple(easiroc, ch, surface, horizontal, vertical));
    //
    //             histName = "hPE";
    //             histAxis = "PE ";
    //             if (surface == EScintiSurface::XY)
    //             {
    //                 histName += "X" + to_string(horizontal) + "Y" + to_string(vertical);
    //                 histAxis += "X=" + to_string(horizontal) + " Y=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
    //             }
    //             else if (surface == EScintiSurface::ZY)
    //             {
    //                 histName += "Z" + to_string(horizontal) + "Y" + to_string(vertical);
    //                 histAxis += "Z=" + to_string(horizontal) + " Y=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
    //             }
    //             else
    //             {
    //                 histName += "X" + to_string(horizontal) + "Z" + to_string(vertical);
    //                 histAxis += "X=" + to_string(horizontal) + " Z=" + to_string(vertical) + ";Light yield (p.e.);Number of events";
    //             }
    //
    //             hPEProtoEach[static_cast<int> (surface)][horizontal][vertical] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);
    //
    //             #ifdef DEBUG
    //                 cout << histName << " created." << endl;
    //             #endif
    //         }
    //     }
    // }


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
    int maxChProtoXOfXZ;
    double maxPEProtoXOfXY;
    double maxPEProtoXOfXZ;

    bool goodEvent;
    bool singleHit;
    bool isStraightBeam;
    bool good1cm;
    bool good2cm;

    bool scintiHit;
    bool hodoScintiHit;

    int countGood = 0;
    int countSingle = 0;
    int countStraight = 0;
    int countGood1cm = 0;
    int countGood2cm = 0;

    int countScintiHit = 0;
    int countHodoScintiHit = 0;

    array<bool, NHodo> hodoHit;

    int gapCount1s = 0;
    int gapCount2s = 0;

    int totalEvt = 0;

    string pt2sfile;
    string pt1sfile;
    string hsfile;


    for (int j = runnum; j < runnum + fileCount; j++)
    {
        // input ADC data
        TString pt2sfile_name = TString::Format("%sproto2s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);
        // TString pt1sfile_name = TString::Format("%sproto1s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);
        TString hsfile_name = TString::Format("%shodo_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, subrun);

        std::string pt2sfile_name_s = std::string(pt2sfile_name);
        // std::string pt1sfile_name_s = std::string(pt1sfile_name);
        std::string hsfile_name_s = std::string(hsfile_name);

        string::size_type pos1 = pt2sfile_name_s.find(".root");
        // string::size_type pos2 = pt1sfile_name_s.find(".root");
        string::size_type pos3 = hsfile_name_s.find(".root");
        string pt2sfile_calib(pt2sfile_name);
        pt2sfile_calib.replace(pos1, 6, "_calib");
        pt2sfile_calib = pt2sfile_calib.substr(13);
        // string pt1sfile_calib(pt1sfile_name);
        // pt1sfile_calib.replace(pos2, 6, "_calib");
        // pt1sfile_calib = pt1sfile_calib.substr(13);
        string hsfile_calib(hsfile_name);
        hsfile_calib.replace(pos3, 6, "_calib");
        hsfile_calib = hsfile_calib.substr(13);

        pt2sfile = pt2sfile_name_s;
        // pt1sfile = pt1sfile_name_s;
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


        // TString pt1s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt1sfile_calib.c_str());
        // ifstream fin2(pt1s_calibname);
        // while (fin2 >> inputCh >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g)
        // {
        //     mean_ped[1][inputCh] = mp;
        //     // sigma_ped[1][inputCh]    = sp;
        //     // chi2_ndf_ped[1][inputCh] = cnp;
        //     // mean_1pe[1][inputCh]     = m1;
        //     // sigma_1pe[1][inputCh]    = s1;
        //     // chi2_ndf_1pe[1][inputCh] = cn1;
        //     gain[1][inputCh] = g;
        // }

        TString hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
        ifstream fin3(hs_calibname);
        if (fin3.fail())
        {
            cerr << "Calibration file does not exist!!!" << endl;
            return;
        }
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
        // TChain* tree2 = new TChain("tree");
        TChain* tree3 = new TChain("tree");

        Int_t adc[NEasiroc][NChEasiroc];
        double pe[NEasiroc][NChEasiroc];

        tree1->Add(pt2sfile.c_str());
        tree1->SetBranchAddress("ADC", &adc[0]);
        // tree2->Add(pt1sfile.c_str());
        // tree2->SetBranchAddress("ADC", &adc[1]);
        tree3->Add(hsfile.c_str());
        tree3->SetBranchAddress("ADC", &adc[2]);

        cout << TString::Format("runnum = %04d", j) << endl;
        cout << "Protype: number of events = " << tree1->GetEntries() << endl;
        cout << "Hodoscope : number of events = " << tree3->GetEntries() << endl;

        if (tree1->GetEntries() == 0 || tree3->GetEntries() == 0)
        {
            cerr << "ROOT files do not exist!!!" << endl;
            return;
        }

        cout << "Select X1:" << "[" << omit_lowx1 + 1 << ", " << 16 - omit_highx1 << "]" << endl;
        cout << "Select Y1:" << "[" << omit_lowy1 + 1 << ", " << 16 - omit_highy1 << "]" << endl;
        cout << "Select Y2:" << "[" << omit_lowy2 + 1 << ", " << 16 - omit_highy2 << "]" << endl;
        cout << "Select X2:" << "[" << omit_lowx2 + 1 << ", " << 16 - omit_highx2 << "]" << endl;

        int less_evt = 0;   // 3つのEASIROCの中で最も少ないイベント数

        less_evt = tree1->GetEntries();

        if (less_evt >= tree3->GetEntries())
        {
            less_evt = tree3->GetEntries();
        }

        E3DOption cubePlaceOption = decideOption(j);
        if (cubePlaceOption == E3DOption::None)
            return;

        else if (cubePlaceOption == E3DOption::WRef2cmP)
        {
            minPEEachCellFit = 15.0;
            maxPEEachCellFit = 35.0;
        }
        else if (cubePlaceOption == E3DOption::WRef1cmV)
        {
            minPEEachCellFit = 3.5;
            maxPEEachCellFit = 17.5;
        }
        else if (cubePlaceOption == E3DOption::WRef2cmV)
        {
            minPEEachCellFit = 9.14 - 4.285;
            maxPEEachCellFit = 9.14 + 9.715;
        }
        else if (cubePlaceOption == E3DOption::WRef1cmP)
        {
            minPEEachCellFit = 8.625 - 4.285;
            maxPEEachCellFit = 8.625 + 9.715;
        }
        if (cubePlaceOption == E3DOption::WORef1cmP || cubePlaceOption == E3DOption::WORef1cmV)
        {
            fitRangeLower = 2;
            fitRangeHigher = 3;
        }



        #ifdef DEBUG
            switch (cubePlaceOption)
            {
                case E3DOption::WRef1cmP:
                    cout << "WRef1cmP" << endl;
                    break;
                case E3DOption::WRef1cmV:
                    cout << "WRef1cmV" << endl;
                    break;
                case E3DOption::WRef2cmP:
                    cout << "WRef2cmP" << endl;
                    break;
                case E3DOption::WRef2cmV:
                    cout << "WRef2cmV" << endl;
                    break;
                case E3DOption::WORef1cmP:
                    cout << "WORef1cmP" << endl;
                    break;
                case E3DOption::WORef1cmV:
                    cout << "WORef1cmV" << endl;
                    break;
                default:
                    break;
            }
        #endif

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
                // tree2->GetEntry(evt);
                tree3->GetEntry(evt);
            }
            else if (evt >= gap_point_0 && evt < gap_point_1)
            {
                tree1->GetEntry(evt - gap_pt2s_0);
                // tree2->GetEntry(evt - gap_pt1s_0);
                tree3->GetEntry(evt - gap_hs_0);
            }
            else if (evt >= gap_point_1)
            {
                tree1->GetEntry(evt - gap_pt2s_1);
                // tree2->GetEntry(evt - gap_pt1s_1);
                tree3->GetEntry(evt - gap_hs_1);
            }


            // Initialize
            hitCountHodo = { };
            trigHodo = { };
            hodoHit = { };
            maxPEHodo = { };
            maxChHodo = { };

            maxChProtoXOfXY = 0;
            maxChProtoXOfXZ = 0;
            maxPEProtoXOfXY = 0.0;
            maxPEProtoXOfXZ = 0.0;

            goodEvent = false;
            singleHit = false;
            isStraightBeam = false;
            good1cm = false;
            good2cm = false;

            scintiHit = false;
            hodoScintiHit = false;

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

            // ホドスコープ全てにヒット→good event
            if (count(trigHodo.begin(), trigHodo.end(), true) == NHodo)
            {
                goodEvent = true;
                countGood++;
            }

            // すべてのホドスコープに1本ずつだけヒット→single hit
            // TODO: single hit イベントを飛ばすかどうか考える
            if (count(hitCountHodo.begin(), hitCountHodo.end(), 1) == NHodo)
            {
                singleHit = true;
                countSingle++;
            }
            else
            {
                // continue;
            }

            // ビームがまっすぐ飛んだ→isStraightBeam
            if (singleHit && maxChHodo[static_cast<int> (EHodoscope::HSX1)] == maxChHodo[static_cast<int> (EHodoscope::HSX2)] && maxChHodo[static_cast<int> (EHodoscope::HSY1)] == maxChHodo[static_cast<int> (EHodoscope::HSY2)])
            {
                isStraightBeam = true;
                countStraight++;
                countHodoHitEachCell[maxChHodo[static_cast<int> (EHodoscope::HSX1)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY1)] - 1]++;
            }


            if (singleHit)
            {
                hHodoHitMapUp->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX1)], maxChHodo[static_cast<int> (EHodoscope::HSY1)]);
                hHodoHitMapDown->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
            }

            if (!goodEvent && singleHit)
            {
                cerr << "!goodEvent && singleHit" << endl;
                return;
            }

            if (isStraightBeam)
            {
                hHodoHitMapWithStraightBeam->Fill(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)]);
            }

            if (singleHit &&
                maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYFor1cm &&
                maxChHodo[static_cast<int> (EHodoscope::HSX1)] >= hodoLowXFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSX1)] <= hodoHighXFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSY1)] >= hodoLowYFor1cm && maxChHodo[static_cast<int> (EHodoscope::HSY1)] <= hodoHighYFor1cm)
            {
                good1cm = true;
                countGood1cm++;
            }
            if (singleHit &&
                maxChHodo[static_cast<int> (EHodoscope::HSX2)] >= hodoLowXFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSX2)] <= hodoHighXFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSY2)] >= hodoLowYFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSY2)] <= hodoHighYFor2cm &&
                maxChHodo[static_cast<int> (EHodoscope::HSX1)] >= hodoLowXFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSX1)] <= hodoHighXFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSY1)] >= hodoLowYFor2cm && maxChHodo[static_cast<int> (EHodoscope::HSY1)] <= hodoHighYFor2cm)
            {
                good2cm = true;
                countGood2cm++;
            }


            // prototype loop
            EEasiroc easiroc = EEasiroc::Scinti2;
            int easirocCh = 9;
            EScintiSurface surface = EScintiSurface::ZY;
            int horizontal = 3;
            int vertical = 3;

            pe[static_cast<int> (easiroc)][easirocCh] = (double(adc[static_cast<int> (easiroc)][easirocCh]) - mean_ped[static_cast<int> (easiroc)][easirocCh]) / gain[static_cast<int> (easiroc)][easirocCh];
            hProto[static_cast<int> (surface)]->SetBinContent(horizontal, vertical, pe[static_cast<int> (easiroc)][easirocCh]);

            if (pe[static_cast<int> (easiroc)][easirocCh] > ProtoPEThreshold)
            {
                hHitRateProto[static_cast<int> (surface)]->Fill(horizontal, vertical);
                hPEProto->Fill(pe[static_cast<int> (easiroc)][easirocCh]);

                scintiHit = true;
                countScintiHit++;
                if (isStraightBeam)
                {
                    countScintiHodoHitEachCell[maxChHodo[static_cast<int> (EHodoscope::HSX1)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY1)] - 1]++;
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

                // 2cm vertical に対するホドスコープ領域ヒットかつシンチヒット
                if (cubePlaceOption == E3DOption::WRef2cmV && good2cm)
                {
                    hodoScintiHit = true;
                    countHodoScintiHit++;
                }
                else if (good1cm)
                {
                    // if(cubePlaceOption == E3DOption::None) cerr << "cube option error!!!" << endl;
                    hodoScintiHit = true;
                    countHodoScintiHit++;
                }
            }

            // hodoscope area good のときのシンチ光量分布
            if (good1cm && is1cmArea(cubePlaceOption))
            {
                hPEProtoGood->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
            }
            if (good2cm && !is1cmArea(cubePlaceOption))
            {
                hPEProtoGood->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
            }

            // セルごとの光量分布
            if (isStraightBeam)
            {
                hPEProtoGoodCell[maxChHodo[static_cast<int> (EHodoscope::HSX2)] - 1][maxChHodo[static_cast<int> (EHodoscope::HSY2)] - 1]->Fill(pe[static_cast<int> (easiroc)][easirocCh]);
            }


            // gap check
            // if (isGap(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)], maxChProtoXOfXY, shiftHSX2, shiftHSY2))
            // {
            //     gapCount1s++;
            // }
            if (isGap(maxChHodo[static_cast<int> (EHodoscope::HSX2)], maxChHodo[static_cast<int> (EHodoscope::HSY2)], scintiHit, cubePlaceOption, shiftHSX2, shiftHSY2))
            {
                gapCount2s++;
            }
            // hGapCount1s->SetBinContent(totalEvt, gapCount1s);
            hGapCount2s->SetBinContent(totalEvt, gapCount2s);


            // Draw event display
            // if (evt >= min_evt && evt <= max_evt && runnum == j && goodEvent && scintiHit)
            if (totalEvt >= min_evt && totalEvt <= max_evt)
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

    cout << "--------------event loop end--------------" << endl;
    cout << "total event: " << totalEvt << endl;
    cout << "good event: " << countGood << endl;
    cout << "single hit event: " << countSingle << endl;
    cout << "straight beam event: " << countStraight << endl;
    cout << "good1cm event: " << countGood1cm << endl;
    cout << "good2cm event: " << countGood2cm << endl;
    cout << "scinti hit: " << countScintiHit << endl;
    cout << "hodoscope & scinti hit: " << countHodoScintiHit << endl;


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
    nHistHori = 3;
    nHistVert = 2;
    hodoHistOrder = { 4, 1, 3, 6 };
    // array<int, NSurfaceScinti> protoHistOrder = { 2, 7 };
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
    canvas->cd(2);
    hPEProto->GetYaxis()->SetTitleSize(TitleSize);
    hPEProto->GetXaxis()->SetTitleSize(TitleSize);
    hPEProto->GetYaxis()->SetLabelSize(LabelSize);
    hPEProto->GetXaxis()->SetLabelSize(LabelSize);
    // hPEProto->GetXaxis()->SetNdivisions(5);
    // hPEProto->GetYaxis()->SetNdivisions(5);
    hPEProto->SetMarkerSize(MarkerSize);
    hPEProto->Draw();

    figName = TString::Format("%sPE_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    // Light yield (each prototype channel)
    canvas = new TCanvas();
    hPEProto->Draw();
    // gPad->SetLogy();
    string chName = "ZY=(3,3)";
    figName = TString::Format("%sPE_%s_%04d_%04d.%s", ScintiPEDir.c_str(), chName.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    // Light yield (scintillator good event)
    TFile histScintiPEGood(TString::Format("%sPEProtoGood_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
    canvas = new TCanvas();
    hPEProtoGood->Draw();
    hPEProtoGood->Write();
    int maximumBin = hPEProtoGood->GetMaximumBin();
    figName = TString::Format("%sPEProtoGood_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    histScintiPEGood.Close();
    canvas->Clear();


    // Light yield (scintillator each cell)
    TFile histsScintiPEEachCell(TString::Format("%sPEEachCell_%04d_%04d.root", CellPEDir.c_str(), runnum, subrun), "RECREATE");
    // TF1* landau = new TF1("");
    gStyle->SetOptFit(111);
    canvas = new TCanvas();
    for (int y = 0; y < NScifiEachHodo; y++)
    {
        for (int x = 0; x < NScifiEachHodo; x++)
        {
            hPEProtoGoodCell[x][y]->Fit("landau", FitOption, "", maximumBin - fitRangeLower, maximumBin + fitRangeHigher);
            TF1* lan = (TF1 *) gROOT->FindObject("landau");
            double mean = hPEProtoGoodCell[x][y]->GetMean();
            double mpv = lan->GetParameter("MPV");
            #ifdef DEBUG
                // cout << "(X,Y)=(" << x+1 << "," << y+1 << "): mean=" << hPEProtoGoodCell[x][y]->GetMean() << ", MPV=" << lan->GetParameter("MPV") << endl;
                cout << TString::Format("(X,Y)") << endl;
            #endif
            if (mean < CutPEForEachCell)
            {
                hPEProtoHodoMap->SetBinContent(x + 1, y + 1, 0);
            }
            else
            {
                hPEProtoHodoMap->SetBinContent(x + 1, y + 1, mpv);
            }
            hPEProtoHodoMapMean->SetBinContent(x + 1, y + 1, mean);
            hPEProtoHodoMapMeanZoom->SetBinContent(x + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, y + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, mean);
            hPEProtoGoodCell[x][y]->Draw();
            hPEProtoGoodCell[x][y]->Write();
            figName = TString::Format("%sPECellXY=(%d,%d)_%04d_%04d.%s", CellPEDir.c_str(), x + 1, y + 1, runnum, subrun, outputFileType.c_str());
            // canvas->SaveAs(figName);
        }
    }
    histsScintiPEEachCell.Close();
    canvas->Clear();


    // Light yield (scintillator eacy Y cell)
    // gStyle->SetOptFit(111);
    if (isVertical(decideOption(runnum)))
    {
        array<int, 2> rangeHSXForEachY;
        array<int, 2> rangeHSYForEachY = { 7, 11 };
        if (is1cmArea(decideOption(runnum)))
        {
            rangeHSXForEachY = { 5, 11 };
        }
        else
        {
            rangeHSXForEachY = { 3, 14 };
        }
        if ((runnum >= 56 && runnum <= 68) || (runnum >= 76 && runnum <= 84))
        {
            rangeHSYForEachY = { 6, 10 };
        }

        TFile histsScintiPEEachY(TString::Format("%sPEEachY_%04d_%04d.root", CellPEDir.c_str(), runnum, subrun), "RECREATE");
        TH1D* hPEProtoGoodY[NScifiEachHodo];
        TGraphErrors* distanceToFiberVsLY = new TGraphErrors();

        canvas = new TCanvas();
        for (int y = rangeHSYForEachY[0] - 1; y < rangeHSYForEachY[1]; y++)
        {
            histName = "hPEEachY_ZY_Y" + to_string(y + 1);
            histAxis = "PE ZY Y = " + to_string(y + 1) + ";Light yield (p.e.);Number of events";
            hPEProtoGoodY[y] = new TH1D(histName.c_str(), histAxis.c_str(), NBinPEProto, MinPEProto, MaxPEProto);

            for (int x = rangeHSXForEachY[0] - 1; x < rangeHSXForEachY[1]; x++)
            {
                hPEProtoGoodY[y]->Add(hPEProtoGoodCell[x][y], 1);
            }
            // hPEProtoGoodY[y]->Fit("landau", FitOption, "", minPEEachCellFit + (y - rangeHSYForEachY[0] + 1) * 0.48, maxPEEachCellFit + (y - rangeHSYForEachY[0] + 1) * 0.48 - 5);
            hPEProtoGoodY[y]->Fit("landau", FitOption, "", maximumBin - fitRangeLower, maximumBin + fitRangeHigher);
            TF1* lan = (TF1 *) gROOT->FindObject("landau");
            double mpv = lan->GetParameter(1);
            distanceToFiberVsLY->SetPoint(y - rangeHSYForEachY[0] + 1, (rangeHSYForEachY[1] - y - 1) * HodoWidth, mpv);
            distanceToFiberVsLY->SetPointError(y - rangeHSYForEachY[0] + 1, HodoWidth / 2, lan->GetParError(1));
            #ifdef DEBUG
                cout << (rangeHSYForEachY[1] - y - 1) << endl;
            #endif
            canvas = new TCanvas();
            hPEProtoGoodY[y]->Draw();
            hPEProtoGoodY[y]->Write();
        }
        histsScintiPEEachY.Close();
        canvas->Clear();

        // Distance of beam hit position to fiber
        TFile graphDisVsLY(TString::Format("%sDisVsLY_%04d_%04d.root", ResultDir.c_str(), runnum, subrun), "RECREATE");
        canvas = new TCanvas();
        distanceToFiberVsLY->SetTitle("Distance of beam hit position to fiber vs. L.Y.");
        distanceToFiberVsLY->GetXaxis()->SetTitle("Distance of beam hit position to fiber (mm)");
        distanceToFiberVsLY->GetYaxis()->SetTitle("Light yield (p.e.)");
        // distanceToFiberVsLY->Fit("pol1", FitOption, "", -0.5, HodoWidth * 5 + 0.5);
        TF1* myPol1 = new TF1("mypol1", "[0]+[1]*x", -0.5, HodoWidth * 5 + 0.5);
        TF1* myExp1 = new TF1("myExp1", "[0]*TMath::Exp(-x/[1])", -0.5, HodoWidth * 5 + 0.5);
        myPol1->SetParameters(8.5, -0.5);
        myExp1->SetParameters(10.0, 12.0);
        // distanceToFiberVsLY->Fit("mypol1", FitOption, "", -0.5, HodoWidth * 5 + 0.5);
        distanceToFiberVsLY->Fit("myExp1", FitOption, "", -0.5, HodoWidth * 5 + 0.5);
        distanceToFiberVsLY->GetXaxis()->SetRangeUser(-0.5, HodoWidth * 5 + 0.5);
        if(runnum <= 74)
        {
            distanceToFiberVsLY->GetYaxis()->SetRangeUser(5, 14);
        }
        else
        {
            distanceToFiberVsLY->GetYaxis()->SetRangeUser(1, 10);
        }
        distanceToFiberVsLY->Draw("AP");
        distanceToFiberVsLY->Write();
        figName = TString::Format("%sDistanceToFiberVsLY_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
        canvas->SaveAs(figName);
        graphDisVsLY.Close();
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


    // gap count
    nHistHori = 1;
    nHistVert = 1;
    histWidth = 800;
    histHeight = 600;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    // canvas->Divide(nHistHori, nHistVert);

    // hGapCount1s->GetXaxis()->SetRangeUser(0, totalEvt);
    hGapCount2s->GetXaxis()->SetRangeUser(0, totalEvt);
    // hGapCount1s->GetXaxis()->SetNdivisions(fileCount, kFALSE);
    hGapCount2s->GetXaxis()->SetNdivisions(fileCount, kFALSE);
    // hGapCount1s->GetXaxis()->SetTickLength(1);
    hGapCount2s->GetXaxis()->SetTickLength(1);
    hGapCount2s->SetStats(kFALSE);

    // canvas->cd(1);
    // hGapCount1s->Draw("P");
    // canvas->cd(2);
    hGapCount2s->Draw("P");
    figName = TString::Format("%sGapCount_%04d_%04d-%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, runnum + fileCount - 1, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();



    // hodoscope hit map
    nHistHori = 3;
    nHistVert = 1;
    histWidth = 640;
    histHeight = 600;
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    canvas->Divide(nHistHori, nHistVert);

    gStyle->SetOptStat(1110);
    const array<double, 4> HodomapStatBoxSize = { 0.65, 0.85, 0.75, 0.93 };

    canvas->cd(1);
    // hHodoHitMapUp->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapUp->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapUp->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapUp->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapUp->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapUp->GetZaxis()->SetTitleOffset(1.3);
    hHodoHitMapUp->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    changeStatsBoxSize(hHodoHitMapUp, HodomapStatBoxSize);

    canvas->cd(2);
    // hHodoHitMapDown->GetYaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetXaxis()->SetTitleSize(0.03);
    // hHodoHitMapDown->GetYaxis()->SetLabelSize(0.04);
    // hHodoHitMapDown->GetXaxis()->SetLabelSize(0.04);
    hHodoHitMapDown->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapDown->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapDown->GetZaxis()->SetTitleOffset(1.3);
    hHodoHitMapDown->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    changeStatsBoxSize(hHodoHitMapDown, HodomapStatBoxSize);


    canvas->cd(3);
    hHodoHitMapWithStraightBeam->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hHodoHitMapWithStraightBeam->GetZaxis()->SetTitleOffset(1.3);
    hHodoHitMapWithStraightBeam->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    gPad->SetRightMargin(RightMarginForHodoMap);
    changeStatsBoxSize(hHodoHitMapWithStraightBeam, HodomapStatBoxSize);

    figName = TString::Format("%sHodoHitMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    // 2D hist のカラーバーを赤色だけにする
    // const Int_t NRGBs = 3;
    // const Int_t NCont = 255;
    //
    // double stops[NRGBs] = { 0.00, .50, 1.00 };
    // double red[NRGBs] = { 1.00, 1.0, 1.00 };
    // double green[NRGBs] = { 1.00, 0.0, 0.00 };
    // double blue[NRGBs] = { 1.00, 0.0, 0.00 };
    // TColor::CreateGradientColorTable(2, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);


// Scintillator PE map
    nHistHori = 1;
    nHistVert = 1;
    gStyle->SetPaintTextFormat("3.2f");
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    hPEProtoHodoMapMean->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hPEProtoHodoMapMean->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hPEProtoHodoMapMean->SetStats(kFALSE);
    hPEProtoHodoMapMean->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    figName = TString::Format("%sScintiPEHodoMapMean_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    const double MarkerSizeHodoMapZoom = 1.8;

    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    hPEProtoHodoMapMeanZoom->GetXaxis()->SetNdivisions(HodoMapZoomedSize1cm);
    hPEProtoHodoMapMeanZoom->GetYaxis()->SetNdivisions(HodoMapZoomedSize1cm);
    hPEProtoHodoMapMeanZoom->SetStats(kFALSE);
    hPEProtoHodoMapMeanZoom->SetMarkerSize(MarkerSizeHodoMapZoom);
    hPEProtoHodoMapMeanZoom->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    figName = TString::Format("%sScintiPEHodoMapMeanZoom_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    hPEProtoHodoMap->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hPEProtoHodoMap->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hPEProtoHodoMap->SetStats(kFALSE);
    hPEProtoHodoMap->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    figName = TString::Format("%sScintiPEHodoMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();


    // Detection efficiency
    // nHistHori = 1;
    // nHistVert = 1;
    // gStyle->SetPaintTextFormat("3.2f");
    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    for (int y = 0; y < NScifiEachHodo; y++)
    {
        for (int x = 0; x < NScifiEachHodo; x++)
        {
            if (countHodoHitEachCell[x][y] == 0)
            {
                hDetectionEff->SetBinContent(x + 1, y + 1, 0.0);
                hDetectionEffZoom->SetBinContent(x + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, y + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, 0.0);
            }
            else
            {
                double de = (double)countScintiHodoHitEachCell[x][y] / (double) countHodoHitEachCell[x][y];
                hDetectionEff->SetBinContent(x + 1, y + 1, de);
                hDetectionEffZoom->SetBinContent(x + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, y + 1 - (NScifiEachHodo - HodoMapZoomedSize1cm)/2.0, de);
            }
            #ifdef DEBUG
                cout << TString::Format("(x,y)=(%2d,%2d): hodo&scinti/hodo = %d/%d", x, y, countScintiHodoHitEachCell[x][y], countHodoHitEachCell[x][y]) << endl;
            #endif
        }
    }
    hDetectionEff->GetXaxis()->SetNdivisions(NScifiEachHodo);
    hDetectionEff->GetYaxis()->SetNdivisions(NScifiEachHodo);
    hDetectionEff->SetStats(kFALSE);
    hDetectionEff->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    figName = TString::Format("%sDetectionEfficiencyHodoMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();

    canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    hDetectionEffZoom->GetXaxis()->SetNdivisions(HodoMapZoomedSize1cm);
    hDetectionEffZoom->GetYaxis()->SetNdivisions(HodoMapZoomedSize1cm);
    hDetectionEffZoom->SetStats(kFALSE);
    hDetectionEffZoom->SetMarkerSize(MarkerSizeHodoMapZoom);
    hDetectionEffZoom->Draw("text colz");
    gPad->SetRightMargin(RightMarginForHodoMap);
    figName = TString::Format("%sDetectionEfficiencyHodoMapZoom_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    canvas->SaveAs(figName);
    canvas->Clear();







    // Beam hit position dependency of cross talk (Hodomap)
    // nHistHori = 2;
    // nHistVert = 1;
    // histWidth = 1200;
    // histHeight = 1200;
    // gStyle->SetPaintTextFormat("3.2f");
    // canvas = new TCanvas("canvas", "", histWidth * nHistHori, histHeight * nHistVert);
    // canvas->Divide(nHistHori, nHistVert);
    // canvas->cd(1);
    // hCrossTalkXYDarkCutMap->Draw("text colz");
    // hCrossTalkXYDarkCutMap->GetXaxis()->SetNdivisions(NScifiEachHodo);
    // hCrossTalkXYDarkCutMap->GetYaxis()->SetNdivisions(NScifiEachHodo);
    // changestatsBoxSize(hCrossTalkXYDarkCutMap, 0.7, 0.9, 0.7, 0.935);
    // drawCubeLine("");
    // canvas->cd(2);
    // hCrossTalkXZDarkCutMap->Draw("text colz");
    // hCrossTalkXZDarkCutMap->GetXaxis()->SetNdivisions(NScifiEachHodo);
    // hCrossTalkXZDarkCutMap->GetYaxis()->SetNdivisions(NScifiEachHodo);
    // changestatsBoxSize(hCrossTalkXZDarkCutMap, 0.7, 0.9, 0.7, 0.935);
    // drawCubeLine("");
    //
    // figName = TString::Format("%sCrossTalkDarkCutMap_%04d_%04d.%s", ResultDir.c_str(), runnum, subrun, outputFileType.c_str());
    // canvas->SaveAs(figName);
    // canvas->Clear();

    cout << endl;
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
    int min_evt, max_evt;
    string fileType;
    if (argc >= 3 + 4)
    {
        x1 = atoi(argv[3]);
        y1 = atoi(argv[4]);
        y2 = atoi(argv[5]);
        x2 = atoi(argv[6]);
        if (argc >= 3 + 4 + 2)
        {
            min_evt = atoi(argv[7]);
            max_evt = atoi(argv[8]);
            if (argc >= 3 + 4 + 2 + 1)
            {
                fileType = argv[9];
            }
        }
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
