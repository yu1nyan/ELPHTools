enum class EScintiSurface
{
    XY, ZY, XZ, None
};

const string SurfaceName[] = {"XY", "ZY", "XZ"};

// EASIROCのCh順に並べる
enum class EHodoscope
{
    HSX1, HSY1, HSY2, HSX2, None
};

const string HodoName[] = {"HSX1", "HSY1", "HSY2", "HSX2"};

enum class EEasiroc
{
    Scinti2, Scinti1, Hodoscope, None
};

enum class EScintiType
{
    WRef3D, WORef3D, OneCube, ThreeCubes, Proto555, Proto446, TwoCubes, NineCubes, Other, None
};
const int NScifi = 64;              //シンチファイバーの本数
const int NHodo = 4;                //ホドスコープの本数
const int NScifiEachHodo = 16;      //ホドスコープ1本あたりのシンチファイバーの本数
const int NScifiEachPlace = 32;     //ホドスコープ1箇所（上流、または下流）あたりのシンチファイバーの本数

const int NEasiroc = 3;
const int NEasirocForScinti = 2;
const int NChEasiroc = 64;
const int NChEasirocHalf = 32;
const int NSurfaceScinti = 3;


// enum class EOutputDir
// {
//     Result, EvtDisplay, ScintiPE, HodoPE
// };
