#pragma once

#include <array>

class ToScintiCh
{
public:
    ToScintiCh(EScintiType scintiType);
    bool isConnected(EEasiroc easirocType, int easirocCh);
    tuple<EScintiSurface, int, int> getCh(EEasiroc easirocType, int easirocCh);

private:
    EScintiType scintiType;
    static const int NChEasirocOneSide = 32;
    // TODO: 他のScintiTypeのときも使うstructは共通で良い？
    typedef struct ScintiCh
    {
        // int xy_x[NChEasirocOneSide];
        // int xy_y[NChEasirocOneSide];
        // int xz_x[NChEasirocOneSide];
        // int xz_z[NChEasirocOneSide];
        // int zy_z[NChEasirocOneSide];
        // int zy_y[NChEasirocOneSide];
        array<int, NChEasirocOneSide> xy_x;
        array<int, NChEasirocOneSide> xy_y;
        array<int, NChEasirocOneSide> xz_x;
        array<int, NChEasirocOneSide> xz_z;
        array<int, NChEasirocOneSide> zy_z;
        array<int, NChEasirocOneSide> zy_y;

    } ScintiCh_t;
    ScintiCh_t proto555Ch;
    // ScintiCh_t oneCubeCh;
    // ScintiCh_t twoCubeCh;
    // ScintiCh_t nineCubeCh;



};
