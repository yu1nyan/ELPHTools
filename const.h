// enum class EScintiSurface
// {
//     XY, ZY, XZ, None
// };

enum class EScintiSurface
{
    ZY, XZ, XY, None
};


// const string SurfaceName[] = {"XY", "ZY", "XZ"};
const string SurfaceName[] = {"ZY", "XZ", "XY"};
const string ReadoutSurfaceName[] = {"X", "Y", "Z"};

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
    Printed, OneCube, ThreeCubes, Proto555, Proto446, TwoCubes, NineCubes, Other, None
};

enum class E3DOption
{
    WRef1cmP, WRef1cmV, WRef2cmP, WRef2cmV, WORef1cmV, WORef1cmP, None
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
const int NScintiOneSide = 5;
const int NCube = 9;
enum class ECubePlace
{
    UpperLeft, Upper, UpperRight, Left, Center, Right, LowerLeft, Lower, LowerRight
};
const string CubeName[] = { "UpperLeft", "Upper", "UpperRight", "Left", "Center", "Right", "LowerLeft", "Lower", "LowerRight" };
const string CubeTitle[] = { "upper left", "upper", "upper right", "left", "center", "right", "lower left", "lower", "lower right" };
const int CubeChMapXY[] = { 11, 12, 28, 8, 9, 26, 5, 6, 24 };

const double HodoWidth = 1.7;   // (mm)


// enum class EOutputDir
// {
//     Result, EvtDisplay, ScintiPE, HodoPE
// };

const double RightMarginForHodoMap = 0.15;
const double TitleOffsetTH2 = 1.6;

const double MarkerSize = 2.2;
