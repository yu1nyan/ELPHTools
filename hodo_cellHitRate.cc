#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <sys/stat.h>

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

using namespace std;

void hodo_cellHitRate(int &runnum, int &nsub, int &file_count, int &u_xgap, int &u_ygap, int &d_xgap, int &d_ygap)
{
    //  gStyle->SetOptStat(1111);
    //  gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(1, 0);

    bool trig_HSX1 = false;
    bool trig_HSY1 = false;
    bool trig_HSY2 = false;
    bool trig_HSX2 = false;

    bool good_event = false;
    bool up_good = false;
    bool down_good = false;
    bool SingleHit = false;

    string rootfile_dir = "../tree_root/";
    string calibfile_dir = "calibfile/";
    string type;

    if (nsub == 3)
        type = "1cube/";
    else
    if (nsub == 4)
        type = "2cubes/";
    else
    if (nsub == 5)
        type = "9cubes/";
    else
    if (nsub == 6)
        type = "3Dprint/";

    string hodoDir = "hodoscope/";
    string resultDir = hodoDir + type;
    mkdir(hodoDir.c_str(), 0777);
    mkdir(resultDir.c_str(), 0777);

    Double_t mean_ped[3][64];
    Double_t mean_1pe[3][64];
    Double_t gain[3][64];
    Double_t sigma_ped[3][64];
    Double_t sigma_1pe[3][64];
    Double_t chi2_ndf_ped[3][64];
    Double_t chi2_ndf_1pe[3][64];

    Double_t pe_h[64];
    Int_t adc_h[64];

    string hsfile;

    int ch;
    double mp, sp, cnp, m1, s1, cn1, g;

    for (int j = runnum; j < runnum + file_count; j++)
    {
        TString hsfile_name = TString::Format("%shodo_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, nsub);

        std::string hsfile_name_s = std::string(hsfile_name);

        string::size_type pos2 = hsfile_name_s.find(".root");
        string hsfile_calib(hsfile_name);
        hsfile_calib.replace(pos2, 6, "_calib");
        hsfile_calib = hsfile_calib.substr(13);

        hsfile = hsfile_name_s;

        TString hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
        ifstream fin2(hs_calibname);
        while (fin2 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g)
        {
            mean_ped[1][ch] = mp;
            sigma_ped[1][ch] = sp;
            chi2_ndf_ped[1][ch] = cnp;
            mean_1pe[1][ch] = m1;
            sigma_1pe[1][ch] = s1;
            chi2_ndf_1pe[1][ch] = cn1;
            gain[1][ch] = g;
        }

        TChain* tree1 = new TChain("tree");

        tree1->Add(hsfile.c_str());
        tree1->SetBranchAddress("ADC", &adc_h);


        int HitCount_HSX1 = 0;
        int HitCount_HSY1 = 0;
        int HitCount_HSY2 = 0;
        int HitCount_HSX2 = 0;

        int max_HSX1 = 0;
        int max_HSY1 = 0;
        int max_HSY2 = 0;
        int max_HSX2 = 0;
        int maxCh_HSX1 = 0;
        int maxCh_HSY1 = 16;
        int maxCh_HSY2 = 32;
        int maxCh_HSX2 = 48;

        int omit_lowx1 = 0;
        int omit_highx1 = 0;
        int omit_lowy1 = 0;
        int omit_highy1 = 0;

        int omit_lowx2 = 0;
        int omit_highx2 = 0;
        int omit_lowy2 = 0;
        int omit_highy2 = 0;

        TString figname1;

        TCanvas* c0 = new TCanvas("c0", "", 2400, 800);

        TH2F* hCellHitRate = new TH2F("", "", 16, 0.5, 16.5, 16, 0.5, 16.5);
        TH2F* huCellHitRate = new TH2F("", "", 16, 0.5, 16.5, 16, 0.5, 16.5);
        TH2F* hdCellHitRate = new TH2F("", "", 16, 0.5, 16.5, 16, 0.5, 16.5);

        for (int y = 0; y < 16; y++)
        {
            for (int x = 0; x < 16; x++)
            {
                omit_lowx1 = 15 - x - u_xgap; // 15-x
                omit_highx1 = x + u_xgap;
                omit_lowy1 = 15 - y - u_ygap;
                omit_highy1 = y + u_ygap;

                omit_lowx2 = x + d_xgap;
                omit_highx2 = 15 - x - d_xgap;
                omit_lowy2 = y + d_ygap;
                omit_highy2 = 15 - y - d_ygap;

                cout << "Hodoscope : number of events = " << tree1->GetEntries() << endl;

                cout << "Select X1:" << "[" << omit_lowx1 + 1 << ", " << 16 - omit_highx1 << "]" << endl;
                cout << "Select Y1:" << "[" << omit_lowy1 + 1 << ", " << 16 - omit_highy1 << "]" << endl;
                cout << "Select Y2:" << "[" << omit_lowy2 + 1 << ", " << 16 - omit_highy2 << "]" << endl;
                cout << "Select X2:" << "[" << omit_lowx2 + 1 << ", " << 16 - omit_highx2 << "]" << endl;

                int good_count = 0;
                int GoodEvent_num = 0;

                for (int evt = 0; evt < tree1->GetEntries(); evt++)
                {
                    tree1->GetEntry(evt);

                    HitCount_HSX1 = 0;
                    HitCount_HSY1 = 0;
                    HitCount_HSY2 = 0;
                    HitCount_HSX2 = 0;

                    trig_HSX1 = false;
                    trig_HSY1 = false;
                    trig_HSY2 = false;
                    trig_HSX2 = false;

                    good_event = false;
                    up_good = false;
                    down_good = false;
                    SingleHit = false;

                    max_HSX1 = 0;
                    max_HSY1 = 0;
                    max_HSY2 = 0;
                    max_HSX2 = 0;
                    maxCh_HSX1 = 0;
                    maxCh_HSY1 = 16;
                    maxCh_HSY2 = 32;
                    maxCh_HSX2 = 48;


                    for (int i = 0; i < 64; i++)
                    {
                        //  pe_h[i] = (double(adc_h[i]));
                        pe_h[i] = (double(adc_h[i]) - mean_ped[1][i]) / gain[1][i];
                        pe_h[i] = double(floor(pe_h[i] * 10) / 10);

                        if (i <= 15)
                        {
                            if (pe_h[i] > 2.5 && i >= (omit_lowx1 + 0) && i <= (15 - omit_highx1))
                            {
                                trig_HSX1 = true;
                            }
                        }
                        else if (i >= 16 && i <= 31)
                        {
                            if (pe_h[i] > 2.5 && i >= (16 + omit_lowy1) && i <= (31 - omit_highy1))
                            {
                                trig_HSY1 = true;
                            }
                        }
                        else if (i >= 32 && i <= 47)
                        {
                            if (pe_h[i] > 2.5 && i >= (32 + omit_lowy2) && i <= (47 - omit_highy2))
                            {
                                trig_HSY2 = true;
                            }
                        }
                        else if (i >= 48 && i <= 63)
                        {
                            if (pe_h[i] > 2.5 && i >= (48 + omit_lowx2) && i <= (63 - omit_highx2))
                            {
                                trig_HSX2 = true;
                            }
                        }
                    }

                    if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2)
                        good_event = true;
                    if (trig_HSX1 && trig_HSY1)
                        up_good = true;
                    if (trig_HSX2 && trig_HSY2)
                        down_good = true;

                    if (good_event)
                    {
                        good_count++;
                        GoodEvent_num = evt;
                        hCellHitRate->Fill(x + 1, y + 1, 1);
                    }
                    if (up_good)
                    {
                        huCellHitRate->Fill(x + 1, y + 1, 1);
                    }
                    if (down_good)
                    {
                        hdCellHitRate->Fill(x + 1, y + 1, 1);
                    }
                }
            } // end of x direction loop
        } // end of y direction loop

        c0->Divide(3, 1);
        c0->cd(3);
        hCellHitRate->GetYaxis()->SetTitleOffset(1.0);
        hCellHitRate->GetYaxis()->SetTitle("Cell# along Y");
        hCellHitRate->GetXaxis()->SetTitle("Cell# along X");
        hCellHitRate->GetYaxis()->SetTitleSize(0.04);
        hCellHitRate->GetXaxis()->SetTitleSize(0.04);
        hCellHitRate->GetYaxis()->SetLabelSize(0.04);
        hCellHitRate->GetXaxis()->SetLabelSize(0.04);
        hCellHitRate->GetXaxis()->SetNdivisions(16);
        hCellHitRate->GetYaxis()->SetNdivisions(16);
        hCellHitRate->SetTitle("hitrate : good events");
        hCellHitRate->SetMarkerSize(1.5);
        hCellHitRate->Draw("text colz");
        hCellHitRate->SetLineColor(2);
        hCellHitRate->SetFillColor(2);
        gPad->RedrawAxis();
        gPad->SetTicks();
        gPad->SetGrid();


        c0->cd(1);
        huCellHitRate->GetYaxis()->SetTitleOffset(1.0);
        huCellHitRate->GetYaxis()->SetTitle("Cell# along Y");
        huCellHitRate->GetXaxis()->SetTitle("Cell# along X");
        huCellHitRate->GetYaxis()->SetTitleSize(0.04);
        huCellHitRate->GetXaxis()->SetTitleSize(0.04);
        huCellHitRate->GetYaxis()->SetLabelSize(0.04);
        huCellHitRate->GetXaxis()->SetLabelSize(0.04);
        huCellHitRate->GetXaxis()->SetNdivisions(16);
        huCellHitRate->GetYaxis()->SetNdivisions(16);
        huCellHitRate->SetTitle("Hitrate : upstream cells");
        huCellHitRate->SetMarkerSize(1.5);
        huCellHitRate->Draw("text colz");
        huCellHitRate->SetLineColor(2);
        huCellHitRate->SetFillColor(2);
        gPad->RedrawAxis();
        gPad->SetTicks();
        gPad->SetGrid();



        c0->cd(2);
        hdCellHitRate->GetYaxis()->SetTitleOffset(1.0);
        hdCellHitRate->GetYaxis()->SetTitle("Cell# along Y");
        hdCellHitRate->GetXaxis()->SetTitle("Cell# along X");
        hdCellHitRate->GetYaxis()->SetTitleSize(0.04);
        hdCellHitRate->GetXaxis()->SetTitleSize(0.04);
        hdCellHitRate->GetYaxis()->SetLabelSize(0.04);
        hdCellHitRate->GetXaxis()->SetLabelSize(0.04);
        hdCellHitRate->GetXaxis()->SetNdivisions(16);
        hdCellHitRate->GetYaxis()->SetNdivisions(16);
        hdCellHitRate->SetTitle("Hitrate : downstream cells");
        hdCellHitRate->SetMarkerSize(1.5);
        hdCellHitRate->Draw("text colz");
        hdCellHitRate->SetLineColor(2);
        hdCellHitRate->SetFillColor(2);
        gPad->RedrawAxis();
        gPad->SetTicks();
        gPad->SetGrid();


        TString hratename = TString::Format("%scell_HitRate_%04d_%04d.pdf", resultDir.c_str(), j, nsub);
        c0->SaveAs(hratename);

        delete hCellHitRate;
        delete huCellHitRate;
        delete hdCellHitRate;
        delete c0;
    } // end of each file loop
}

int main(int argc, char** argv)
{
    if (argc != 8)
    {
        cerr << "hodo_cellHitRate <1st run#> <nsub> <# of files> <u_xgap> <u_ygap> <d_xgap> <d_ygap>" <<
            endl;

        return -1;
    }

    int runnum = atoi(argv[1]);
    int nsub = atoi(argv[2]);
    int file_count = atoi(argv[3]);
    int u_xgap = atoi(argv[4]);
    int u_ygap = atoi(argv[5]);
    int d_xgap = atoi(argv[6]);
    int d_ygap = atoi(argv[7]);

    hodo_cellHitRate(runnum, nsub, file_count,
                     u_xgap, u_ygap, d_xgap, d_ygap);

    return 0;
}
