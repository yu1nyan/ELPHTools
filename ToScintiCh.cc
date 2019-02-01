#include <iostream>
#include <tuple>
// #include <fstream>
// #include <sstream>
// #include <cstdlib>
// #include <array>
// #include <ctime>

// include ROOT headers
// #include <TFile.h>
// #include <TStyle.h>
#include <TROOT.h>
// #include <TColor.h>
#include <TMath.h>

using namespace std;

#include "const.h"
#include "ToScintiCh.h"




// コンストラクタでシンチレーターのタイプ（プロトタイプか、1キューブか、など）を与える
ToScintiCh::ToScintiCh(EScintiType scintiType)
{
    this->scintiType = scintiType;

    if (scintiType == EScintiType::Proto555)
    {
        // 対応するシンチChがない場合は0を入れる

        // easiroc ch:      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
        proto555Ch.xy_x = { 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5 };
        proto555Ch.xy_y = { 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
        proto555Ch.zy_z = { 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5 };
        proto555Ch.zy_y = { 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
        // proto555Ch.zy_z = {0, 5, 4, 3, 5, 4, 3, 5, 4, 3, 5, 4, 3, 5, 4, 3, 0, 0, 0, 0, 0, 0, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1};
        // proto555Ch.zy_y = {0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};

        // easiroc ch:     32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63
        proto555Ch.xz_x = { 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5 };
        proto555Ch.xz_z = { 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
    }
}

bool ToScintiCh::isConnected(EEasiroc easirocType, int easirocCh)
{
    if (easirocType == EEasiroc::Hodoscope)
    {
        // cerr << "Invalid EASIROC type!!" << endl;

        return true;
    }
    else if (scintiType == EScintiType::Proto555)
    {
        if (easirocCh <= -1 || easirocCh >= 64)
        {
            return false;
        }
        else if (easirocType == EEasiroc::Scinti2)
        {
            return true;
        }
        else if (easirocType == EEasiroc::Scinti1 && easirocCh <= 31)
        {
            return true;
        }
        return false;
    }
    else if (scintiType == EScintiType::TwoCubes)
    {
        if (easirocType == EEasiroc::Scinti2 && (easirocCh == 9 || easirocCh == 40 || easirocCh == 41))
        {
            return true;
        }
        if (easirocType == EEasiroc::Scinti1 && (easirocCh == 8 || easirocCh == 9))
        {
            return true;
        }
        return false;
    }
    else if (scintiType == EScintiType::OneCube)
    {
        if (easirocType == EEasiroc::Scinti2 && (easirocCh == 9 || easirocCh == 41))
        {
            return true;
        }
        if (easirocType == EEasiroc::Scinti1 &&  easirocCh == 9)
        {
            return true;
        }
        return false;
    }
    else if (scintiType == EScintiType::NineCubes)
    {
        if (easirocType == EEasiroc::Scinti2)
        {
            switch (easirocCh)
            {
                case 40:
                case 41:
                case 58:
                case 12:
                case 9:
                case 6:
                    return true;
                default:
                    return false;
            }
        }
        if (easirocType == EEasiroc::Scinti1)
        {
            switch (easirocCh)
            {
                case 5:
                case 6:
                case 24:
                case 8:
                case 9:
                case 26:
                case 11:
                case 12:
                case 28:
                    return true;
                default:
                    return false;
            }
        }
        return false;
    }
    else if (scintiType == EScintiType::Printed)
    {
        if (easirocCh == 9 && easirocType == EEasiroc::Scinti2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// どのEASIROCのどのチャンネルかという情報を与えると、どの平面のどの位置かというのを返す関数
// 1個目のint: 水平方向の位置
// 2個目のint: 垂直方向の位置
tuple<EScintiSurface, int, int> ToScintiCh::getCh(EEasiroc easirocType, int easirocCh)
{
    if (!isConnected(easirocType, easirocCh) || easirocType == EEasiroc::Hodoscope)
    {
        return forward_as_tuple(EScintiSurface::None, 0, 0);
    }
    else if (scintiType == EScintiType::Proto555)
    {
        if (easirocType == EEasiroc::Scinti2)
        {
            if (easirocCh <= 31)
            {
                if (proto555Ch.zy_z[easirocCh] == 0 || proto555Ch.zy_y[easirocCh] == 0)
                {
                    return forward_as_tuple(EScintiSurface::None, 0, 0);
                }

                return forward_as_tuple(EScintiSurface::ZY, proto555Ch.zy_z[easirocCh], proto555Ch.zy_y[easirocCh]);
            }
            else
            {
                if (proto555Ch.xz_x[easirocCh - NScifiEachPlace] == 0 || proto555Ch.xz_z[easirocCh - NScifiEachPlace] == 0)
                {
                    return forward_as_tuple(EScintiSurface::None, 0, 0);
                }

                return forward_as_tuple(EScintiSurface::XZ, proto555Ch.xz_x[easirocCh - NScifiEachPlace], proto555Ch.xz_z[easirocCh - NScifiEachPlace]);
            }
        }
        // else if (easirocType == EEasiroc::Scinti1 && easirocCh <= 31)
        else
        {
            if (proto555Ch.xy_x[easirocCh] == 0 || proto555Ch.xy_y[easirocCh] == 0)
            {
                return forward_as_tuple(EScintiSurface::None, 0, 0);
            }

            return forward_as_tuple(EScintiSurface::XY, proto555Ch.xy_x[easirocCh], proto555Ch.xy_y[easirocCh]);
        }
    }
    else if (scintiType == EScintiType::TwoCubes)
    {
        if (easirocType == EEasiroc::Scinti2)
        {
            switch (easirocCh)
            {
                // case 9: return forward_as_tuple(EScintiSurface::ZY, 1, 1);
                // case 40: return forward_as_tuple(EScintiSurface::XZ, 1, 1);
                // case 41: return forward_as_tuple(EScintiSurface::XZ, 2, 1);
                case 9: return forward_as_tuple(EScintiSurface::ZY, 3, 3);
                case 40: return forward_as_tuple(EScintiSurface::XZ, 2, 3);
                case 41: return forward_as_tuple(EScintiSurface::XZ, 3, 3);
                default: return forward_as_tuple(EScintiSurface::None, 0, 0);
            }
        }
        else
        {
            switch (easirocCh)
            {
                case 8: return forward_as_tuple(EScintiSurface::XY, 2, 3);
                case 9: return forward_as_tuple(EScintiSurface::XY, 3, 3);
                default: return forward_as_tuple(EScintiSurface::None, 0, 0);
            }
        }
    }
    else if (scintiType == EScintiType::Printed)
    {
        if (easirocCh == 9 && easirocType == EEasiroc::Scinti2)
        {
            return forward_as_tuple(EScintiSurface::ZY, 3, 3);
        }
    }

    // cerr << "EASIROC Ch conversion error! (EASIROC -> Scintillator)" << endl;
    return forward_as_tuple(EScintiSurface::None, 0, 0);
}
