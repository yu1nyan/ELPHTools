#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>

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

void gap_check_2E(int &runnum, int &nsub, int &file_count, int &min_evt, int &max_evt, int &gap_point_0, int &gap_pt_0, int &gap_hs_0, int &gap_point_1, int &gap_pt_1, int &gap_hs_1)
{
    gStyle->SetOptStat(1111);
    //  gStyle->SetOptFit(1111);

    TF1 *f1 = new TF1("landau_orig", "[0]*TMath::Landau(x,[1],[2])");

    int map[25] = { 31, 30, 15, 14, 13, // SLOT1 A1,B1,C1,D1,E1
                    29, 28, 12, 11, 10, //      A2,B2,C2,D2,E2
                    27, 26, 9,  8, 7, //      A3,B3,C3,D3,E3
                    25, 24, 6,  5, 4, //      A4,B4,C4,D4,E4
                    23, 22, 3,  2, 1 }; //      A5,B5,C5,D5,E5

    int Bij_X[9] = { 1, 1, 1,
                     1, 2, 3,
                     1, 2, 3 };

    int Bij_Y[9] = { 1, 2, 3,
                     1, 1, 1,
                     1, 2, 3 };

    Double_t binmin = -3.; // before calib : 700. / after calib : 0.
    Double_t binmax = 30.; // before calib : 2000. / after calib : 30

    TH2F *hHSX1 = new TH2F("hHSX1", "HSX1", 16, 0.5, 16.5, 1, 0.5, 1.5);
    hHSX1->SetMinimum(binmin);
    hHSX1->SetMaximum(binmax);
    TH2F *hHSY1 = new TH2F("hHSY1", "HSY1", 1, 0.5, 1.5, 16, 0.5, 16.5);
    hHSY1->SetMinimum(binmin);
    hHSY1->SetMaximum(binmax);
    TH2F *hHSY2 = new TH2F("hHSY2", "HSY2", 1, 0.5, 1.5, 16, 0.5, 16.5);
    hHSY2->SetMinimum(binmin);
    hHSY2->SetMaximum(binmax);
    TH2F *hHSX2 = new TH2F("hHSX2", "HSX2", 16, 0.5, 16.5, 1, 0.5, 1.5);
    hHSX2->SetMinimum(binmin);
    hHSX2->SetMaximum(binmax);
    TH2F *hCubeX = new TH2F("hCubeX", "Several cubes setup : X readout (event-by-event)",
                            1, 0.5, 1.5, 3, 0.5, 3.5);
    hCubeX->SetMinimum(binmin);
    hCubeX->SetMaximum(binmax + 20);
    TH2F *hCubeY = new TH2F("hCubeY", "Several cubes setup : Y readout (event-by-event)",
                            3, 0.5, 3.5, 1, 0.5, 1.5);
    hCubeY->SetMinimum(binmin);
    hCubeY->SetMaximum(binmax + 20);
    TH2F *hCubeZ = new TH2F("hCubeZ", "Several cubes setup : Z readout (event-by-event)",
                            3, 0.5, 3.5, 3, 0.5, 3.5);
    hCubeZ->SetMinimum(binmin);
    hCubeZ->SetMaximum(binmax + 20);

    TH1F *hHitRate_HSX1 = new TH1F("hHitRate_HSX1", "HitRate_HSX1", 16, 0.5, 16.5);
    hHitRate_HSX1->SetMinimum(0);
    TH1F *hHitRate_HSY1 = new TH1F("hHitRate_HSY1", "HitRate_HSY1", 16, 0.5, 16.5);
    hHitRate_HSY1->SetMinimum(0);
    TH1F *hHitRate_HSY2 = new TH1F("hHitRate_HSY2", "HitRate_HSY2", 16, 0.5, 16.5);
    hHitRate_HSY2->SetMinimum(0);
    TH1F *hHitRate_HSX2 = new TH1F("hHitRate_HSX2", "HitRate_HSX2", 16, 0.5, 16.5);
    hHitRate_HSX2->SetMinimum(0);

    TH2F *hHitRate_X = new TH2F("hHitRate_X", "HitRate_X", 1, 0.5, 1.5, 3, 0.5, 3.5);
    TH2F *hHitRate_Y = new TH2F("hHitRate_Y", "HitRate_Y", 3, 0.5, 3.5, 1, 0.5, 1.5);
    TH2F *hHitRate_Z = new TH2F("hHitRate_Z", "HitRate_Z", 3, 0.5, 3.5, 3, 0.5, 3.5);

    TH1F *hNHit_HSX1 = new TH1F("hNHit_HSX1", "Number of Hits per Event (HSX1)", 16, 0.5, 16.5);
    TH1F *hNHit_HSY1 = new TH1F("hNHit_HSY1", "Number of Hits per Event (HSY1)", 16, 0.5, 16.5);
    TH1F *hNHit_HSY2 = new TH1F("hNHit_HSY2", "Number of Hits per Event (HSY2)", 16, 0.5, 16.5);
    TH1F *hNHit_HSX2 = new TH1F("hNHit_HSX2", "Number of Hits per Event (HSX2)", 16, 0.5, 16.5);


    TCanvas *c0 = new TCanvas("c0", "", 1800, 600);

    TH1F    *hPE[25]      = { };
    TH1F    *hPE_h[64]    = { };
    Double_t array_PE[25] = { };
    //  Double_t array_PE_h[64] = {};

    for (int x = 0; x < 25; x++) {
        TString histname1 = TString::Format("Position No.%d ", x);
        // hPE[x] = new TH1F(histname1, histname1, 1800,700,2500);
        hPE[x] = new TH1F(histname1, histname1, 200, -5, 95);
    }
    for (int x = 0; x < 64; x++) {
        TString histname2 = TString::Format("Hodoscope Ch%d ", x);
        // hPE_h[x] = new TH1F(histname2, histname2, 1800,700,2500);
        hPE_h[x] = new TH1F(histname2, histname2, 200, -5, 95);
    }

    TH1F *hDist = new TH1F("hDist", "", 40, 0., 20.);

    bool trig_HSX1        = false;
    bool trig_HSY1        = false;
    bool trig_HSY2        = false;
    bool trig_HSX2        = false;
    bool trig_X           = false;
    bool trig_Y           = false;
    bool trig_Z           = false;
    int  all_count_x1     = 0;
    int  all_count_y1     = 0;
    int  all_count_y2     = 0;
    int  all_count_x2     = 0;
    int  trig_count       = 0;
    int  andtrig_count    = 0;
    int  nontrig_count    = 0;
    int  andnontrig_count = 0;

    bool good_event = false;
    bool DpassEvent = false;

    string rootfile_dir   = "../tree_root/";
    string calibfile_dir  = "calibfile/";
    string evtdisplay_dir = "evtdisplay_dir/";

    Double_t mean_ped[3][64];
    Double_t mean_1pe[3][64];
    Double_t gain[3][64];
    Double_t sigma_ped[3][64];
    Double_t sigma_1pe[3][64];
    Double_t chi2_ndf_ped[3][64];
    Double_t chi2_ndf_1pe[3][64];

    Double_t pe_p[64];
    Double_t pe_h[64];
    Int_t    adc_p[64];
    Int_t    adc_h[64];

    int less_evt_sum    = 0;
    int good_count_sum  = 0;
    int Dpass_count_sum = 0;

    string ptfile;
    string hsfile;

    int    ch;
    double mp, sp, cnp, m1, s1, cn1, g;

    string type_dir;
    if (nsub == 3) { type_dir = "1cube/"; }
    else if (nsub == 4) { type_dir = "2cubes/"; }
    else if (nsub == 5) { type_dir = "9cubes/"; }
    else if (nsub == 6) { type_dir = "3Dprint/"; }

    for (int j = runnum; j < runnum + file_count; j++) {
        TString ptfile_name = TString::Format("%sproto_%04d_%04d_tree2.root",
                                              rootfile_dir.c_str(), j, nsub);
        TString hsfile_name = TString::Format("%shodo_%04d_%04d_tree2.root",
                                              rootfile_dir.c_str(), j, nsub);
        /*
        if (j == ) {
          gap_point_0 = ;
          gap_pt2s_0 = ;
          gap_pt1s_0 = ;
          gap_hs_0 = ;
          gap_point_1 = ;
          gap_pt2s_1 = ;
          gap_pt1s_1 = ;
          gap_hs_1 = ;
        } else if {
          gap_point_0 = ;
          gap_pt2s_0 = ;
          gap_pt1s_0 = ;
          gap_hs_0 = ;
          gap_point_1 = ;
          gap_pt2s_1 = ;
          gap_pt1s_1 = ;
          gap_hs_1 = ;
        }
        */

        std::string ptfile_name_s = std::string(ptfile_name);
        std::string hsfile_name_s = std::string(hsfile_name);

        string::size_type pos1 = ptfile_name_s.find(".root");
        string::size_type pos2 = hsfile_name_s.find(".root");
        string            ptfile_calib(ptfile_name);
        ptfile_calib.replace(pos1, 6, "_calib");
        ptfile_calib = ptfile_calib.substr(13);
        string hsfile_calib(hsfile_name);
        hsfile_calib.replace(pos2, 6, "_calib");
        hsfile_calib = hsfile_calib.substr(13);

        ptfile = ptfile_name_s;
        hsfile = hsfile_name_s;

        TString  pt_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), ptfile_calib.c_str());
        ifstream fin1(pt_calibname);
        while (fin1 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
            mean_ped[0][ch]     = mp;
            sigma_ped[0][ch]    = sp;
            chi2_ndf_ped[0][ch] = cnp;
            mean_1pe[0][ch]     = m1;
            sigma_1pe[0][ch]    = s1;
            chi2_ndf_1pe[0][ch] = cn1;
            gain[0][ch]         = g;
        }

        TString  hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
        ifstream fin2(hs_calibname);
        while (fin2 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
            mean_ped[1][ch]     = mp;
            sigma_ped[1][ch]    = sp;
            chi2_ndf_ped[1][ch] = cnp;
            mean_1pe[1][ch]     = m1;
            sigma_1pe[1][ch]    = s1;
            chi2_ndf_1pe[1][ch] = cn1;
            gain[1][ch]         = g;
        }

        TChain *tree1 = new TChain("tree");
        TChain *tree2 = new TChain("tree");

        tree1->Add(ptfile.c_str());
        tree1->SetBranchAddress("ADC", &adc_p);
        tree2->Add(hsfile.c_str());
        tree2->SetBranchAddress("ADC", &adc_h);


        int HitCount_HSX1 = 0;
        int HitCount_HSY1 = 0;
        int HitCount_HSY2 = 0;
        int HitCount_HSX2 = 0;

        int max_HSX1   = 0;
        int max_HSY1   = 0;
        int max_HSY2   = 0;
        int max_HSX2   = 0;
        int maxCh_HSX1 = 0;
        int maxCh_HSY1 = 16;
        int maxCh_HSY2 = 32;
        int maxCh_HSX2 = 48;

        double dist_X = 0;
        double dist_Y = 0;
        double dist   = 0;

        TString figname1;

        int omit_lowx1  = 0;
        int omit_highx1 = 0;
        int omit_lowy1  = 0;
        int omit_highy1 = 0;

        int omit_lowx2  = 0;
        int omit_highx2 = 0;
        int omit_lowy2  = 0;
        int omit_highy2 = 0;

        cout << "Protype : number of events = " << tree1->GetEntries() << endl;
        cout << "Hodoscope : number of events = " << tree2->GetEntries() << endl;

        cout << "Select X1:" << "[" << omit_lowx1 + 1 << ", " << 16 - omit_highx1 << "]" << endl;
        cout << "Select Y1:" << "[" << omit_lowy1 + 1 << ", " << 16 - omit_highy1 << "]" << endl;
        cout << "Select Y2:" << "[" << omit_lowy2 + 1 << ", " << 16 - omit_highy2 << "]" << endl;
        cout << "Select X2:" << "[" << omit_lowx2 + 1 << ", " << 16 - omit_highx2 << "]" << endl;

        int less_evt          = 0;
        int good_count        = 0;
        int Dpass_count       = 0;
        int DpassEvent_num[2] = { };
        int GoodEvent_num     = 0;

        if (tree1->GetEntries() >= tree2->GetEntries()) { less_evt = tree2->GetEntries(); }
        else { less_evt = tree1->GetEntries(); }

        for (int evt = 0; evt < less_evt; evt++) {
            if (evt < gap_point_0) {
                tree1->GetEntry(evt);
                tree2->GetEntry(evt);
            } else if (evt >= gap_point_0 && evt < gap_point_1) {
                tree1->GetEntry(evt - gap_pt_0);
                tree2->GetEntry(evt - gap_hs_0);
            } else if (evt >= gap_point_1) {
                tree1->GetEntry(evt - gap_pt_1);
                tree2->GetEntry(evt - gap_hs_1);
            }

            HitCount_HSX1 = 0;
            HitCount_HSY1 = 0;
            HitCount_HSY2 = 0;
            HitCount_HSX2 = 0;

            trig_HSX1 = false;
            trig_HSY1 = false;
            trig_HSY2 = false;
            trig_HSX2 = false;

            trig_X     = false;
            trig_Y     = false;
            trig_Z     = false;
            DpassEvent = false;
            good_event = false;

            max_HSX1   = 0;
            max_HSY1   = 0;
            max_HSY2   = 0;
            max_HSX2   = 0;
            maxCh_HSX1 = 0;
            maxCh_HSY1 = 16;
            maxCh_HSY2 = 32;
            maxCh_HSX2 = 48;

            dist_X = 0;
            dist_Y = 0;
            dist   = 0;

            figname1 = TString::Format("%s/evtdisplay/hit_%d_%d.pdf", type_dir.c_str(), runnum, evt);

            for (int i = 0; i < 64; i++) {
                //  pe_h[i] = (double(adc_h[i]));
                pe_h[i] = (double(adc_h[i]) - mean_ped[1][i]) / gain[1][i];
                pe_h[i] = double(floor(pe_h[i] * 10) / 10);

                if (i <= 15) {
                    hHSX1->SetBinContent(16 - i, 1, pe_h[i]);
                    if (pe_h[i] > 7.5 &&
                        i >= (omit_lowx1 + 0) && i <= (15 - omit_highx1)) {
                        hHitRate_HSX1->Fill(16 - i);
                        HitCount_HSX1++;
                        all_count_x1++;
                        trig_HSX1 = true;
                        max_HSX1  = pe_h[0];
                        if (pe_h[i] > max_HSX1) { maxCh_HSX1 = 15 - i; }
                    }
                }
                else if (i >= 16 && i <= 31) {
                    hHSY1->SetBinContent(1, 32 - i, pe_h[i]);
                    if (pe_h[i] > 7.5 &&
                        i >= (16 + omit_lowy1) && i <= (31 - omit_highy1)) {
                        hHitRate_HSY1->Fill(32 - i);
                        HitCount_HSY1++;
                        all_count_y1++;
                        trig_HSY1 = true;
                        max_HSY1  = pe_h[16];
                        if (pe_h[i] > max_HSY1) { maxCh_HSY1 = 31 - i; }
                    }
                }
                else if (i >= 32 && i <= 47) {
                    hHSY2->SetBinContent(1, i - 31, pe_h[i]);
                    if (pe_h[i] > 7.5 &&
                        i >= (32 + omit_lowy2) && i <= (47 - omit_highy2)) {
                        hHitRate_HSY2->Fill(i - 31);
                        HitCount_HSY2++;
                        all_count_y2++;
                        trig_HSY2 = true;
                        max_HSY2  = pe_h[32];
                        if (pe_h[i] > max_HSY2) { maxCh_HSY2 = i - 32; }
                    }
                }
                else if (i >= 48 && i <= 63) {
                    hHSX2->SetBinContent(i - 47, 1, pe_h[i]);
                    if (pe_h[i] > 7.5 &&
                        i >= (48 + omit_lowx2) && i <= (63 - omit_highx2)) {
                        hHitRate_HSX2->Fill(i - 47);
                        HitCount_HSX2++;
                        all_count_x2++;
                        trig_HSX2 = true;
                        max_HSX2  = pe_h[48];
                        if (pe_h[i] > max_HSX2) { maxCh_HSX2 = i - 48; }
                    }
                }
            }

            hNHit_HSX1->Fill(HitCount_HSX1);
            hNHit_HSY1->Fill(HitCount_HSY1);
            hNHit_HSY2->Fill(HitCount_HSY2);
            hNHit_HSX2->Fill(HitCount_HSX2);

            if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2) { good_event = true; }


            for (int i = 0; i < 32; i++) {
                // pe_t[map[i]] = (double(adc_t[map[i]]));
                pe_p[map[i]] = (double(adc_p[map[i]]) - mean_ped[0][map[i]]) / gain[0][map[i]];
                pe_p[map[i]] = double(floor(pe_p[map[i]] * 10) / 10);

                if (i <= 2) {
                    hCubeX->SetBinContent(1, i + 1, pe_p[map[i]]);
                    if (pe_p[map[i]] > 9.5) {
                        hHitRate_X->Fill(1, i + 1);
                        trig_X = true;
                    }
                }
                if (i >= 3 && i <= 5) {
                    hCubeY->SetBinContent(i - 2, 1, pe_p[map[i]]);
                    if (pe_p[map[i]] > 9.5) {
                        hHitRate_Y->Fill(i - 2, 1);
                        trig_Y = true;
                    }
                }
                if (i >= 6 && i <= 14) {
                    hCubeZ->SetBinContent(Bij_X[i], Bij_Y[i], pe_p[map[i]]);
                    if (pe_p[map[i]] > 9.5) {
                        hHitRate_Z->Fill(Bij_X[i], Bij_Y[i]);
                        trig_Z = true;
                    }
                }
            }

            if (trig_X && trig_Y && trig_Z) {
                Dpass_count++;
                DpassEvent = true;
                // cout << "D-passing event# : " << evt << endl;
                if (DpassEvent_num[0] == 0) { DpassEvent_num[0] = evt; }
                else {
                    DpassEvent_num[1] = DpassEvent_num[0];
                    DpassEvent_num[0] = evt;
                }
                if (DpassEvent_num[1] - GoodEvent_num < 0) {
                    cout << "Diff between PT & HS (bigger): " << DpassEvent_num[0] - GoodEvent_num << endl;
                }
            }

            if (trig_HSX1 || trig_HSY1) { trig_count++; }
            if (trig_HSX1 && trig_HSY1) { andtrig_count++; }
            if (trig_HSX2 || trig_HSY2) { nontrig_count++; }
            if (trig_HSX2 && trig_HSY2) { andnontrig_count++; }

            if (good_event) {
                good_count++;
                GoodEvent_num = evt;
                // cout << "GOOD event# : " << evt << endl;

                cout << "Diff between PT & HS (smaller): " << DpassEvent_num[0] - GoodEvent_num << endl;

                dist_X = abs(maxCh_HSX2 - maxCh_HSX1);
                dist_Y = abs(maxCh_HSY2 - maxCh_HSY1);
                dist   = sqrt(pow(dist_X, 2) + pow(dist_Y, 2));

                hDist->Fill(dist);

                if (evt >= min_evt && evt <= max_evt) {
                    gStyle->SetOptStat(0);
                    c0->Clear();
                    c0->Divide(4, 2);

                    c0->cd(5);
                    hHSX1->GetYaxis()->SetTitle();
                    hHSX1->GetXaxis()->SetTitle();
                    hHSX1->GetYaxis()->SetTitleSize(0.03);
                    hHSX1->GetXaxis()->SetTitleSize(0.03);
                    hHSX1->GetYaxis()->SetLabelSize(0.04);
                    hHSX1->GetXaxis()->SetLabelSize(0.04);
                    hHSX1->GetXaxis()->SetNdivisions(16);
                    hHSX1->GetYaxis()->SetNdivisions(0);
                    hHSX1->SetMarkerSize(2.0);
                    hHSX1->Draw("text60 colz");

                    c0->cd(1);
                    hHSY1->GetYaxis()->SetTitle();
                    hHSY1->GetXaxis()->SetTitle();
                    hHSY1->GetYaxis()->SetTitleSize(0.03);
                    hHSY1->GetXaxis()->SetTitleSize(0.03);
                    hHSY1->GetYaxis()->SetLabelSize(0.04);
                    hHSY1->GetXaxis()->SetLabelSize(0.04);
                    hHSY1->GetXaxis()->SetNdivisions(0);
                    hHSY1->GetYaxis()->SetNdivisions(16);
                    hHSY1->SetMarkerSize(2.0);
                    hHSY1->Draw("text colz");

                    c0->cd(4);
                    hHSY2->GetYaxis()->SetTitle();
                    hHSY2->GetXaxis()->SetTitle();
                    hHSY2->GetYaxis()->SetTitleSize(0.03);
                    hHSY2->GetXaxis()->SetTitleSize(0.03);
                    hHSY2->GetYaxis()->SetLabelSize(0.04);
                    hHSY2->GetXaxis()->SetLabelSize(0.04);
                    hHSY2->GetXaxis()->SetNdivisions(0);
                    hHSY2->GetYaxis()->SetNdivisions(16);
                    hHSY2->SetMarkerSize(2.0);
                    hHSY2->Draw("text colz");

                    c0->cd(8);
                    hHSX2->GetYaxis()->SetTitle();
                    hHSX2->GetXaxis()->SetTitle();
                    hHSX2->GetYaxis()->SetTitleSize(0.03);
                    hHSX2->GetXaxis()->SetTitleSize(0.03);
                    hHSX2->GetYaxis()->SetLabelSize(0.04);
                    hHSX2->GetXaxis()->SetLabelSize(0.04);
                    hHSX2->GetXaxis()->SetNdivisions(16);
                    hHSX2->GetYaxis()->SetNdivisions(0);
                    hHSX2->SetMarkerSize(2.0);
                    hHSX2->Draw("text60 colz");

                    c0->cd(2);
                    hCubeX->GetYaxis()->SetTitle();
                    hCubeX->GetXaxis()->SetTitle();
                    hCubeX->GetYaxis()->SetTitleSize(0.03);
                    hCubeX->GetXaxis()->SetTitleSize(0.03);
                    hCubeX->GetYaxis()->SetLabelSize(0.04);
                    hCubeX->GetXaxis()->SetLabelSize(0.04);
                    hCubeX->GetXaxis()->SetNdivisions(0);
                    hCubeX->GetYaxis()->SetNdivisions(3);
                    hCubeX->SetMarkerSize(2.0);
                    hCubeX->Draw("text colz");

                    c0->cd(7);
                    hCubeY->GetYaxis()->SetTitle();
                    hCubeY->GetXaxis()->SetTitle();
                    hCubeY->GetYaxis()->SetTitleSize(0.03);
                    hCubeY->GetXaxis()->SetTitleSize(0.03);
                    hCubeY->GetYaxis()->SetLabelSize(0.04);
                    hCubeY->GetXaxis()->SetLabelSize(0.04);
                    hCubeY->GetXaxis()->SetNdivisions(3);
                    hCubeY->GetYaxis()->SetNdivisions(0);
                    hCubeY->SetMarkerSize(2.0);
                    hCubeY->Draw("text colz");

                    c0->cd(3);
                    hCubeZ->GetYaxis()->SetTitle();
                    hCubeZ->GetXaxis()->SetTitle();
                    hCubeZ->GetYaxis()->SetTitleSize(0.03);
                    hCubeZ->GetXaxis()->SetTitleSize(0.03);
                    hCubeZ->GetYaxis()->SetLabelSize(0.04);
                    hCubeZ->GetXaxis()->SetLabelSize(0.04);
                    hCubeZ->GetXaxis()->SetNdivisions(3);
                    hCubeZ->GetYaxis()->SetNdivisions(3);
                    hCubeZ->SetMarkerSize(2.0);
                    hCubeZ->Draw("text colz");

                    c0->SaveAs(figname1);
                }
            }
        }


        TString  figname3 = TString::Format("%sdistance_%d.pdf", type_dir.c_str(), runnum);
        TCanvas *c1       = new TCanvas("c1", "", 800, 600);
        gStyle->SetOptStat(1111);
        hDist->GetYaxis()->SetTitleOffset(1.0);
        hDist->GetYaxis()->SetTitle("# of events");
        hDist->GetXaxis()->SetTitle("distance [cells]");
        hDist->GetYaxis()->SetTitleSize(0.05);
        hDist->GetXaxis()->SetTitleSize(0.05);
        hDist->GetYaxis()->SetLabelSize(0.03);
        hDist->GetXaxis()->SetLabelSize(0.04);
        hDist->GetXaxis()->SetNdivisions(511);
        hDist->GetYaxis()->SetNdivisions(511);
        hDist->SetTitle("Distance between up-hit cell and down-hit cell");
        hDist->Draw();

        c1->SaveAs(figname3);

        cout << all_count_x1  << " events passed X1" << endl;
        cout << all_count_y1  << " events passed Y1" << endl;
        cout << all_count_y2  << " events passed Y2" << endl;
        cout << all_count_x2  << " events passed X2" << endl;

        cout << "DAQ trigger counts : " << trig_count << endl;
        cout << "1st layer counts : " << HitCount_HSX1 << endl;
        cout << "non-trigger hodoscope OR counts : " << nontrig_count << endl;
        cout << "trigger hodoscope AND counts : " << andtrig_count << endl;
        cout << "non-trigger hodoscope AND counts : " << andnontrig_count << endl;

        less_evt_sum    += less_evt;
        good_count_sum  += good_count;
        Dpass_count_sum += Dpass_count;

        cout << good_count << " good events out of " << less_evt <<
            " successfully processed" << endl;
        cout << "(in total : " << good_count_sum << " good events out of " << less_evt_sum <<
            " )" << endl;
        cout << Dpass_count << " D-passing events out of " << less_evt <<
            " successfully processed" << endl;
        cout << "(in total : " << Dpass_count_sum << " D-passing events out of " << less_evt_sum <<
            " )" << endl;
    } // end of each file loop

    TString  text5 = TString::Format("%s75ch/array_PE.tsv", type_dir.c_str());
    ofstream outputtext5(text5);
    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    TString  pdfname1;
    gStyle->SetOptFit(1111);

    for (int x = 0; x < 75; x++) {
        f1->SetParameters(10, 60, 1.);
        hPE[x]->Fit("landau_orig", "", "", 30, 75);
        array_PE[x] = f1->GetParameter(1);
        pdfname1    = TString::Format("%s75ch/run%d_PE_position%d.pdf", type_dir.c_str(), runnum, x);

        hPE[x]->GetYaxis()->SetTitle("Number of Events");
        hPE[x]->GetXaxis()->SetTitle("p.e.");
        hPE[x]->GetYaxis()->SetTitleSize(0.04);
        hPE[x]->GetXaxis()->SetTitleSize(0.04);
        hPE[x]->GetYaxis()->SetLabelSize(0.05);
        hPE[x]->GetXaxis()->SetLabelSize(0.05);
        hPE[x]->Draw();
        c2->SetLogy();
        c2->SaveAs(pdfname1);
        c2->Clear();

        outputtext5 << x << " " << array_PE[x] << endl;
    }

    TCanvas *c3 = new TCanvas("c3", "", 1800, 600);
    c3->Clear();
    c3->Divide(4, 2);

    c3->cd(3);
    hHitRate_X->GetXaxis()->SetTitle("Cube# along X");
    hHitRate_X->GetYaxis()->SetTitle("Cube# along Y");
    hHitRate_X->GetYaxis()->SetTitleSize(0.03);
    hHitRate_X->GetXaxis()->SetTitleSize(0.03);
    hHitRate_X->GetYaxis()->SetLabelSize(0.04);
    hHitRate_X->GetXaxis()->SetLabelSize(0.04);
    hHitRate_X->GetXaxis()->SetNdivisions(0);
    hHitRate_X->GetYaxis()->SetNdivisions(3);
    hHitRate_X->SetMarkerSize(2.0);
    hHitRate_X->Draw("text colz");

    c3->cd(7);
    hHitRate_Y->GetXaxis()->SetTitle("Cube# along X");
    hHitRate_Y->GetYaxis()->SetTitle("Cube# along Z");
    hHitRate_Y->GetYaxis()->SetTitleSize(0.03);
    hHitRate_Y->GetXaxis()->SetTitleSize(0.03);
    hHitRate_Y->GetYaxis()->SetLabelSize(0.04);
    hHitRate_Y->GetXaxis()->SetLabelSize(0.04);
    hHitRate_Y->GetXaxis()->SetNdivisions(3);
    hHitRate_Y->GetYaxis()->SetNdivisions(0);
    hHitRate_Y->SetMarkerSize(2.0);
    hHitRate_Y->Draw("text colz");

    c3->cd(2);
    hHitRate_Z->GetYaxis()->SetTitle("Cube# along Y");
    hHitRate_Z->GetXaxis()->SetTitle("Cube# along Z");
    hHitRate_Z->GetYaxis()->SetTitleSize(0.03);
    hHitRate_Z->GetXaxis()->SetTitleSize(0.03);
    hHitRate_Z->GetYaxis()->SetLabelSize(0.04);
    hHitRate_Z->GetXaxis()->SetLabelSize(0.04);
    hHitRate_Z->GetXaxis()->SetNdivisions(3);
    hHitRate_Z->GetYaxis()->SetNdivisions(3);
    hHitRate_Z->SetMarkerSize(2.0);
    hHitRate_Z->Draw("text colz");


    c3->cd(5);
    hHitRate_HSX1->GetYaxis()->SetTitle();
    hHitRate_HSX1->GetXaxis()->SetTitle();
    hHitRate_HSX1->GetYaxis()->SetTitleSize(0.03);
    hHitRate_HSX1->GetXaxis()->SetTitleSize(0.03);
    hHitRate_HSX1->GetYaxis()->SetLabelSize(0.04);
    hHitRate_HSX1->GetXaxis()->SetLabelSize(0.04);
    hHitRate_HSX1->GetXaxis()->SetNdivisions(16);
    hHitRate_HSX1->SetMarkerSize(2.0);
    hHitRate_HSX1->Draw();

    c3->cd(1);
    hHitRate_HSY1->GetYaxis()->SetTitle();
    hHitRate_HSY1->GetXaxis()->SetTitle();
    hHitRate_HSY1->GetYaxis()->SetTitleSize(0.03);
    hHitRate_HSY1->GetXaxis()->SetTitleSize(0.03);
    hHitRate_HSY1->GetYaxis()->SetLabelSize(0.04);
    hHitRate_HSY1->GetXaxis()->SetLabelSize(0.04);
    hHitRate_HSY1->GetXaxis()->SetNdivisions(16);
    hHitRate_HSY1->SetMarkerSize(2.0);
    hHitRate_HSY1->Draw("hbar");

    c3->cd(4);
    hHitRate_HSY2->GetYaxis()->SetTitle();
    hHitRate_HSY2->GetXaxis()->SetTitle();
    hHitRate_HSY2->GetYaxis()->SetTitleSize(0.03);
    hHitRate_HSY2->GetXaxis()->SetTitleSize(0.03);
    hHitRate_HSY2->GetYaxis()->SetLabelSize(0.04);
    hHitRate_HSY2->GetXaxis()->SetLabelSize(0.04);
    hHitRate_HSY2->GetXaxis()->SetNdivisions(16);
    hHitRate_HSY2->SetMarkerSize(2.0);
    hHitRate_HSY2->Draw("hbar");

    c3->cd(8);
    hHitRate_HSX2->GetYaxis()->SetTitle();
    hHitRate_HSX2->GetXaxis()->SetTitle();
    hHitRate_HSX2->GetYaxis()->SetTitleSize(0.03);
    hHitRate_HSX2->GetXaxis()->SetTitleSize(0.03);
    hHitRate_HSX2->GetYaxis()->SetLabelSize(0.04);
    hHitRate_HSX2->GetXaxis()->SetLabelSize(0.04);
    hHitRate_HSX2->GetXaxis()->SetNdivisions(16);
    hHitRate_HSX2->SetMarkerSize(2.0);
    hHitRate_HSX2->Draw();

    TString hitratename = TString::Format("%shitRate_allRun.pdf", type_dir.c_str());
    c3->SaveAs(hitratename);


    TCanvas *c4 = new TCanvas("c4", "", 800, 600);
    c4->Clear();
    c4->Divide(2, 2);

    c4->cd(3);
    hNHit_HSX1->GetYaxis()->SetTitle("# of Events");
    hNHit_HSX1->GetXaxis()->SetTitle("# of Hits");
    hNHit_HSX1->GetYaxis()->SetTitleSize(0.03);
    hNHit_HSX1->GetXaxis()->SetTitleSize(0.03);
    hNHit_HSX1->GetYaxis()->SetLabelSize(0.04);
    hNHit_HSX1->GetXaxis()->SetLabelSize(0.04);
    hNHit_HSX1->Draw();

    c4->cd(1);
    hNHit_HSY1->GetYaxis()->SetTitle("# of Events");
    hNHit_HSY1->GetXaxis()->SetTitle("# of Hits");
    hNHit_HSY1->GetYaxis()->SetTitleSize(0.03);
    hNHit_HSY1->GetXaxis()->SetTitleSize(0.03);
    hNHit_HSY1->GetYaxis()->SetLabelSize(0.04);
    hNHit_HSY1->GetXaxis()->SetLabelSize(0.04);
    hNHit_HSY1->Draw();

    c4->cd(2);
    hNHit_HSY2->GetYaxis()->SetTitle("# of Events");
    hNHit_HSY2->GetXaxis()->SetTitle("# of Hits");
    hNHit_HSY2->GetYaxis()->SetTitleSize(0.03);
    hNHit_HSY2->GetXaxis()->SetTitleSize(0.03);
    hNHit_HSY2->GetYaxis()->SetLabelSize(0.04);
    hNHit_HSY2->GetXaxis()->SetLabelSize(0.04);
    hNHit_HSY2->Draw();

    c4->cd(4);
    hNHit_HSX2->GetYaxis()->SetTitle("# of Events");
    hNHit_HSX2->GetXaxis()->SetTitle("# of Hits");
    hNHit_HSX2->GetYaxis()->SetTitleSize(0.03);
    hNHit_HSX2->GetXaxis()->SetTitleSize(0.03);
    hNHit_HSX2->GetYaxis()->SetLabelSize(0.04);
    hNHit_HSX2->GetXaxis()->SetLabelSize(0.04);
    hNHit_HSX2->Draw();

    TString multiname = TString::Format("%shodoscope_HitMultiplicity_allRun.pdf", type_dir.c_str());
    c4->SaveAs(multiname);
}

int main(int argc, char **argv)
{
    if (argc != 12) {
        cerr << "gap_check_2E <1st run#> <sub#> <# of files> <min_evt> <max_evt> <1st gap point & gap> <2nd>" <<
            endl;

        return -1;
    }

    int runnum      = atoi(argv[1]);
    int nsub        = atoi(argv[2]);
    int file_count  = atoi(argv[3]);
    int min_evt     = atoi(argv[4]);
    int max_evt     = atoi(argv[5]);
    int gap_point_0 = atoi(argv[6]);
    int gap_pt_0    = atoi(argv[7]);
    int gap_hs_0    = atoi(argv[8]);
    int gap_point_1 = atoi(argv[9]);
    int gap_pt_1    = atoi(argv[10]);
    int gap_hs_1    = atoi(argv[11]);

    gap_check_2E(runnum, nsub, file_count, min_evt, max_evt,
                 gap_point_0, gap_pt_0, gap_hs_0,
                 gap_point_1, gap_pt_1, gap_hs_1);

    return 0;
}
