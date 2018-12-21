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

void gap_check_3E(int &runnum, int &file_count, int &nsub, int &gap_point_0, int &gap_pt2s_0, int &gap_pt1s_0, int &gap_hs_0, int &gap_point_1, int &gap_pt2s_1, int &gap_pt1s_1, int &gap_hs_1)
{
    gStyle->SetOptStat(1111);
    //  gStyle->SetOptFit(1111);

    int map[75] = { 31, 30, 15, 14, 13, // SLOT1 A1,B1,C1,D1,E1 (For readout X)
                    29, 28, 12, 11, 10, //      A2,B2,C2,D2,E2
                    27, 26, 9,  8, 7, //      A3,B3,C3,D3,E3
                    25, 24, 6,  5, 4, //      A4,B4,C4,D4,E4
                    23, 22, 3,  2, 1, //      A5,B5,C5,D5,E5
                    63, 62, 47, 46, 45, // SLOT2 A1,B1,C1,D1,E1 (For readout Y)
                    61, 60, 44, 43, 42, //      A2,B2,C2,D2,E2
                    59, 58, 41, 40, 39, //      A3,B3,C3,D3,E3
                    57, 56, 38, 37, 36, //      A4,B4,C4,D4,E4
                    55, 54, 35, 34, 33, //      A5,B5,C5,D5,E5
                    31, 30, 15, 14, 13, // SLOT3 A1,B1,C1,D1,E1 (For readout Z)
                    29, 28, 12, 11, 10, //      A2,B2,C2,D2,E2
                    27, 26, 9,  8, 7, //      A3,B3,C3,D3,E3
                    25, 24, 6,  5, 4, //      A4,B4,C4,D4,E4
                    23, 22, 3,  2, 1 }; //      A5,B5,C5,D5,E5

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);

    bool trig_HSX1 = false;
    bool trig_HSY1 = false;
    bool trig_HSY2 = false;
    bool trig_HSX2 = false;

    bool good_event_hs = false;
    bool trig_up_2s    = false;
    bool trig_down_2s  = false;
    bool good_event_2s = false;
    bool trig_1s       = false;
    bool good_event_1s = false;

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

    Double_t pe_pt2s[64];
    Double_t pe_pt1s[64];
    Double_t pe_h[64];
    Int_t    adc_pt2s[64];
    Int_t    adc_pt1s[64];
    Int_t    adc_h[64];

    string pt2sfile;
    string pt1sfile;
    string hsfile;

    int    ch;
    double mp, sp, cnp, m1, s1, cn1, g;


    for (int j = runnum; j < runnum + file_count; j++) {
        TString pt2sfile_name = TString::Format("%sproto2s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, nsub);
        TString pt1sfile_name = TString::Format("%sproto1s_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, nsub);
        TString hsfile_name = TString::Format("%shodo_%04d_%04d_tree2.root", rootfile_dir.c_str(), j, nsub);

        TString  gaptextname = TString::Format("gap_check/3E/gap_check_%d.tsv", j);
        ofstream outputfile(gaptextname);

        std::string pt2sfile_name_s = std::string(pt2sfile_name);
        std::string pt1sfile_name_s = std::string(pt1sfile_name);
        std::string hsfile_name_s   = std::string(hsfile_name);

        string::size_type pos1 = pt2sfile_name_s.find(".root");
        string::size_type pos2 = pt1sfile_name_s.find(".root");
        string::size_type pos3 = hsfile_name_s.find(".root");
        string            pt2sfile_calib(pt2sfile_name);
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
        hsfile   = hsfile_name_s;

        TString pt2s_calibname = TString::Format("%s%s.tsv",
                                                 calibfile_dir.c_str(), pt2sfile_calib.c_str());
        ifstream fin1(pt2s_calibname);
        while (fin1 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
            mean_ped[0][ch]     = mp;
            sigma_ped[0][ch]    = sp;
            chi2_ndf_ped[0][ch] = cnp;
            mean_1pe[0][ch]     = m1;
            sigma_1pe[0][ch]    = s1;
            chi2_ndf_1pe[0][ch] = cn1;
            gain[0][ch]         = g;
        }

        TString pt1s_calibname = TString::Format("%s%s.tsv",
                                                 calibfile_dir.c_str(), pt1sfile_calib.c_str());
        ifstream fin2(pt1s_calibname);
        while (fin2 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
            mean_ped[1][ch]     = mp;
            sigma_ped[1][ch]    = sp;
            chi2_ndf_ped[1][ch] = cnp;
            mean_1pe[1][ch]     = m1;
            sigma_1pe[1][ch]    = s1;
            chi2_ndf_1pe[1][ch] = cn1;
            gain[1][ch]         = g;
        }

        TString  hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
        ifstream fin3(hs_calibname);
        while (fin3 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
            mean_ped[2][ch]     = mp;
            sigma_ped[2][ch]    = sp;
            chi2_ndf_ped[2][ch] = cnp;
            mean_1pe[2][ch]     = m1;
            sigma_1pe[2][ch]    = s1;
            chi2_ndf_1pe[2][ch] = cn1;
            gain[2][ch]         = g;
        }

        TChain *tree1 = new TChain("tree");
        TChain *tree2 = new TChain("tree");
        TChain *tree3 = new TChain("tree");

        tree1->Add(pt2sfile.c_str());
        tree1->SetBranchAddress("ADC", &adc_pt2s);
        tree2->Add(pt1sfile.c_str());
        tree2->SetBranchAddress("ADC", &adc_pt1s);
        tree3->Add(hsfile.c_str());
        tree3->SetBranchAddress("ADC", &adc_h);

        int omit_lowx1  = 0;
        int omit_highx1 = 0;
        int omit_lowy1  = 0;
        int omit_highy1 = 0;

        int omit_lowx2  = 0;
        int omit_highx2 = 0;
        int omit_lowy2  = 0;
        int omit_highy2 = 0;

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

        cout << "Protype 2s: number of events = " << tree1->GetEntries() << endl;
        cout << "Protype 1s: number of events = " << tree2->GetEntries() << endl;
        cout << "Hodoscope : number of events = " << tree3->GetEntries() << endl;

        cout << "Select X1:" << "[" << omit_lowx1 + 1 << ", " << 16 - omit_highx1 << "]" << endl;
        cout << "Select Y1:" << "[" << omit_lowy1 + 1 << ", " << 16 - omit_highy1 << "]" << endl;
        cout << "Select Y2:" << "[" << omit_lowy2 + 1 << ", " << 16 - omit_highy2 << "]" << endl;
        cout << "Select X2:" << "[" << omit_lowx2 + 1 << ", " << 16 - omit_highx2 << "]" << endl;

        int less_evt          = 0;
        int Dpass_count       = 0;
        int good_count        = 0;
        int all_trig_num      = 0;
        int DpassEvent_num[4] = { };
        int GoodEvent_num     = 0;

        int trig_count       = 0;
        int andtrig_count    = 0;
        int nontrig_count    = 0;
        int andnontrig_count = 0;

        TH1F *hDist = new TH1F("hDist", "", 40, 0., 20.);

        if (tree1->GetEntries() >= tree2->GetEntries()) { less_evt = tree2->GetEntries(); }
        else { less_evt = tree1->GetEntries(); }

        if (less_evt >= tree3->GetEntries()) { less_evt = tree3->GetEntries(); }

        for (int evt = 0; evt < less_evt; evt++) {
            if (evt < gap_point_0) {
                tree1->GetEntry(evt);
                tree2->GetEntry(evt);
                tree3->GetEntry(evt);
            } else if (evt >= gap_point_0 && evt < gap_point_1) {
                tree1->GetEntry(evt - gap_pt2s_0);
                tree2->GetEntry(evt - gap_pt1s_0);
                tree3->GetEntry(evt - gap_hs_0);
            } else if (evt >= gap_point_1) {
                tree1->GetEntry(evt - gap_pt2s_1);
                tree2->GetEntry(evt - gap_pt1s_1);
                tree3->GetEntry(evt - gap_hs_1);
            }

            trig_HSX1 = false;
            trig_HSY1 = false;
            trig_HSY2 = false;
            trig_HSX2 = false;

            trig_up_2s    = false;
            trig_down_2s  = false;
            trig_1s       = false;
            good_event_hs = false;
            good_event_2s = false;
            good_event_1s = false;
            max_HSX1      = 0;
            max_HSY1      = 0;
            max_HSY2      = 0;
            max_HSX2      = 0;
            maxCh_HSX1    = 0;
            maxCh_HSY1    = 16;
            maxCh_HSY2    = 32;
            maxCh_HSX2    = 48;

            dist_X = 0;
            dist_Y = 0;
            dist   = 0;

            for (int i = 0; i < 64; i++) {
                //  pe_h[i] = (double(adc_h[i]));
                pe_h[i] = (double(adc_h[i]) - mean_ped[2][i]) / gain[2][i];
                pe_h[i] = double(floor(pe_h[i] * 10) / 10);

                if (i <= 15) {
                    if (pe_h[i] > 7.5 && i >= (omit_lowx1 + 0) && i <= (15 - omit_highx1)) {
                        trig_HSX1 = true;
                        max_HSX1  = pe_h[0];
                        if (pe_h[i] > max_HSX1) { maxCh_HSX1 = 15 - i; }
                    }
                }
                else if (i >= 16 && i <= 31) {
                    if (pe_h[i] > 7.5 && i >= (16 + omit_lowy1) && i <= (31 - omit_highy1)) {
                        trig_HSY1 = true;
                        max_HSY1  = pe_h[16];
                        if (pe_h[i] > max_HSY1) { maxCh_HSY1 = 31 - i; }
                    }
                }
                else if (i >= 32 && i <= 47) {
                    if (pe_h[i] > 7.5 && i >= (32 + omit_lowy2) && i <= (47 - omit_highy2)) {
                        trig_HSY2 = true;
                        max_HSY2  = pe_h[32];
                        if (pe_h[i] > max_HSY2) { maxCh_HSY2 = i - 32; }
                    }
                }
                else if (i >= 48 && i <= 63) {
                    if (pe_h[i] > 7.5 && i >= (48 + omit_lowx2) && i <= (63 - omit_highx2)) {
                        trig_HSX2 = true;
                        max_HSX2  = pe_h[48];
                        if (pe_h[i] > max_HSX2) { maxCh_HSX2 = i - 48; }
                    }
                }
            }

            if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2) { good_event_hs = true; }

            for (int i = 0; i < 75; i++) {
                // pe_pt2s[map[i]] = (double(adc_pt2s[map[i]]));
                if (i < 50) {
                    pe_pt2s[map[i]] = (double(adc_pt2s[map[i]]) - mean_ped[0][map[i]]) / gain[0][map[i]];
                    pe_pt2s[map[i]] = double(floor(pe_pt2s[map[i]] * 10) / 10);

                    if (i < 25) {
                        if (i % 5 == 0 && pe_pt2s[map[i]] > 9.5) {
                            trig_up_2s = true;
                        }
                        if (i % 5 == 4 && pe_pt2s[map[i]] > 9.5) {
                            trig_down_2s = true;
                        }
                    }
                } else {
                    // pe_pt1s[map[i]] = (double(adc_pt1s[map[i]]));
                    pe_pt1s[map[i]] = (double(adc_pt1s[map[i]]) - mean_ped[1][map[i]]) / gain[1][map[i]];
                    pe_pt1s[map[i]] = double(floor(pe_pt1s[map[i]] * 10) / 10);

                    if ( (i == 56 || i == 57 || i == 58 || i == 61 || i == 62 || i == 63 || i == 66 || i == 67 || i == 68) && pe_pt1s[map[i]] > 39.5) {
                        trig_1s = true;
                    }
                }
            }

            if (trig_up_2s && trig_down_2s) {
                good_event_2s = true;
                Dpass_count++;
                // outputfile << "2s D-passing event# : " << evt << endl;
                if (DpassEvent_num[0] == 0) { DpassEvent_num[0] = evt; }
                else {
                    DpassEvent_num[1] = DpassEvent_num[0];
                    DpassEvent_num[0] = evt;
                }
                if (DpassEvent_num[1] - GoodEvent_num < 0) {
                    outputfile << "Diff between PT & HS (bigger:2s): " <<
                        DpassEvent_num[0] - GoodEvent_num << endl;
                }
            }


            if (trig_1s) {
                good_event_1s = true;
                // outputfile << "1s D-passing event# : " << evt << endl;
                if (DpassEvent_num[2] == 0) { DpassEvent_num[2] = evt; }
                else {
                    DpassEvent_num[3] = DpassEvent_num[2];
                    DpassEvent_num[2] = evt;
                }
                if (DpassEvent_num[3] - GoodEvent_num < 0) {
                    outputfile << "Diff between PT & HS (bigger:1s): " <<
                        DpassEvent_num[2] - GoodEvent_num << endl;
                }
            }

            if (trig_HSX1 || trig_HSY1) { trig_count++; }
            if (trig_HSX1 && trig_HSY1) { andtrig_count++; }
            if (trig_HSX2 || trig_HSY2) { nontrig_count++; }
            if (trig_HSX2 && trig_HSY2) { andnontrig_count++; }

            if (good_event_hs) {
                outputfile << "GOOD event# : " << evt << endl;
                good_count++;
                GoodEvent_num = evt;
                outputfile << "Diff between PT & HS (smaller): " << DpassEvent_num[0] - GoodEvent_num <<
                    " / " << DpassEvent_num[2] - GoodEvent_num << endl;
                all_trig_num++;

                dist_X = abs(maxCh_HSX2 - maxCh_HSX1);
                dist_Y = abs(maxCh_HSY2 - maxCh_HSY1);
                dist   = sqrt(pow(dist_X, 2) + pow(dist_Y, 2));

                hDist->Fill(dist);
            }
        }


        TString figname3 = TString::Format("gap_check/3E/distance_%d.pdf", j);
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
        delete hDist;

        cout << "run# = " << j << endl;
        cout << "DAQ trigger counts : " << trig_count << endl;
        cout << "non-trigger hodoscope OR counts : " << nontrig_count << endl;
        cout << "trigger hodoscope AND counts : " << andtrig_count << endl;
        cout << "non-trigger hodoscope AND counts : " << andnontrig_count << endl;

        cout << good_count << " good events out of " << tree1->GetEntries() <<
            " successfully processed" << endl;
        cout << Dpass_count << " D-passing events out of " << tree1->GetEntries() <<
            " successfully processed" << endl;
    } // end of run loop
}

int main(int argc, char **argv)
{
    if (argc != 12) {
        cerr << "gap_check_3E <runnum> <nsub> <# of files> <gap_point_0 & gap> <gap_point_1 & gap>" << endl;

        return -1;
    }

    int runnum      = atoi(argv[1]);
	int nsub		= atoi(argv[2]);
    int file_count  = atoi(argv[3]);
    int gap_point_0 = atoi(argv[4]);
    int gap_pt2s_0  = atoi(argv[5]);
    int gap_pt1s_0  = atoi(argv[6]);
    int gap_hs_0    = atoi(argv[7]);
    int gap_point_1 = atoi(argv[8]);
    int gap_pt2s_1  = atoi(argv[9]);
    int gap_pt1s_1  = atoi(argv[10]);
    int gap_hs_1    = atoi(argv[11]);

    gap_check_3E(runnum, file_count, nsub, gap_point_0, gap_pt2s_0, gap_pt1s_0, gap_hs_0, gap_point_1, gap_pt2s_1, gap_pt1s_1, gap_hs_1);

    return 0;
}
