#pragma once

// input directories
const string rootfile_dir = "../tree_root/";
const string calibfile_dir = "../calibfile/";
const string evtdisplay_dir = "evtdisplay_dir/";
const string sameHistsDir = "sameHists/";

// time_t rawTime = time(nullptr);
// struct tm timeInfo = *localtime(&rawTime);
// string dateTimeStr = to_string(timeInfo.tm_sec + 100*timeInfo.tm_min + 10000*timeInfo.tm_hour + 1000000*timeInfo.tm_mday + 100000000*(timeInfo.tm_mon+1) + 10000000000*(1900+timeInfo.tm_year));



// string HodoName(int place)
// {
//     switch(place)
//     {
//         case 0: return "HSY1";
//         case 1: return "HSX1";
//         case 2: return "HSX2";
//         case 3: return "HSY2";
//         default: return "";
//     }
// }

enum class EScintiSurface
{
        XY, ZY, XZ, None
};

enum class EHodoscope
{
        HSY1, HSX1, HSX2, HSY2, None
};

enum class EEasiroc
{
        Hodoscope, Scinti1, Scinti2, None
};

enum class EScintiType
{
        WRef3D, WORef3D, OneCube, ThreeCubes, Proto555, Proto446, TwoCubes, NineCubes, Other, None
};

const int ChHodoUpX = 0;
const int ChHodoUpY = 16;
const int ChHodoDownX = 48;
const int ChHodoDownY = 32;

const int NScifi = 64;              //シンチファイバーの本数
const int NHodo = 4;                //ホドスコープの本数
const int NScifiEachHodo = 16;      //ホドスコープ1本あたりのシンチファイバーの本数
const int NScifiEachPlace = 32;     //ホドスコープ1箇所（上流、または下流）あたりのシンチファイバーの本数

const int NTargets = 4;

// const Double_t ThresholdPEHodo = 3.5; //ホドスコープヒストグラムに入れるpeのしきい値
const double ThresholdPESci = -5;    // 各シンチレーターごとの光量分布に使うしきい値
const double ThresholdPESciNHit = 2.5;  // ヒットチャンネル数分布に用いるしきい値
const double ThresholdPEHodoIndi[NHodo] = {3.5, 3.5, 3.5, 3.5};  // HSX2, HSY2, HSY1, HSX1

const double DistanceOfBetHodo = 480;   // unit: mm
const double IntervalOfHodo = 1.7;
const Double_t DistanceOfBetHodoMax = TMath::Sqrt(TMath::Power(DistanceOfBetHodo, 2) + TMath::Power(IntervalOfHodo * NScifiEachHodo * TMath::Sqrt(2), 2)) - 1.0;
const Int_t NBinDistanceHodoHist = 64;
//

// Double_t binmin = ThresholdPEHodo;
const Double_t binmax = 70.;    // イベントディスプレイ　ホドスコープ2次元ヒストグラムの最大値

const Double_t NHitChPerEvtMax = 6.5;
const Int_t NHitChPerEvtNBin = 6;
const Double_t PETotMax = 80.;
const Int_t PETotNBin = 80;


const Double_t PEIndiMin = -5;
const Double_t PEIndiMax = 30.;
const Int_t PEIndiNBin = (PEIndiMax - PEIndiMin)*10;


const Int_t NRGBs = 3;
const Int_t NCont = 255;

const int NDataPEFiberDis = 6;

const Double_t PEFiberDisDrawMinWRef[NTargets] = {3.5, 5.2, 9.3, 7.3};
const Double_t PEFiberDisDrawMaxWRef[NTargets] = {5, 6.8, 11, 8.8};
const Double_t PEFiberDisDrawMinWORef[NTargets] = {1.5, 1.5, 1.8, 1.7};
const Double_t PEFiberDisDrawMaxWORef[NTargets] = {2.2, 2.5, 2.8, 2.5};
const Double_t DistanceErrorAll = 0;




// uniformity関連
// const int UniOmitX1[NTargets] = {5, 5, 5, 5};
// const int UniOmitX2[NTargets] = {5, 5, 0, 0};
// const int UniOmitY1[NTargets] = {5, 5, 5, 5};
// const int UniOmitY2[NTargets] = {5, 5, 5, 5};
const int UniOmitX1[NTargets] = {};
const int UniOmitX2[NTargets] = {};
const int UniOmitY1[NTargets] = {};
const int UniOmitY2[NTargets] = {};
const Double_t PEUniMin = -5;
const Double_t PEUniMax = 30;
const Int_t PEUniNBin = PEUniMax - PEUniMin;
const double ThresholdUniformity = ThresholdPESci;  // ホドスコープ何pe以上でヒストグラムに詰めるか

const Double_t FitMinUniWRef[NTargets] = {2, 3, 5, 5};
const Double_t FitMaxUniWRef[NTargets] = {10, 11, 16, 14};
const Double_t FitMinUniWORef[NTargets] = {1.5, 1.5, 1.5, 1.5};
const Double_t FitMaxUniWORef[NTargets] = {5,5,5,5};

const Double_t PE2DMapMinWRef = 0;
const Double_t PE2DMapMinWORef = -2;
const Double_t PE2DMapMaxWRef = 12;
const Double_t PE2DMapMaxWORef = 7;

const Double_t Chi2NDFMaxWRef = 2;
const Double_t Chi2NDFMaxWORef = 30;

const double ThresholdDE[NTargets] = {2.5, 2.5, 2.5, 2.5};
const Double_t DEHistMax = 1.0;
const Double_t DEHistMin = 0.0;
