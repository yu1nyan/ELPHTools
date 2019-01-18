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

void run_proto(int &runnum, int &min_evt, int &max_evt,
	       int &gap_point_0, int &gap_pt2s_0, int &gap_pt1s_0, int &gap_hs_0,
	       int &gap_point_1, int &gap_pt2s_1, int &gap_pt1s_1, int &gap_hs_1) {

  gStyle->SetOptStat(1111);
  //  gStyle->SetOptFit(1111);

  TF1 *f1 = new TF1("landau_orig","[0]*TMath::Landau(x,[1],[2])");

  int map[75] = {31,30,15,14,13,  //SLOT1 A1,B1,C1,D1,E1 (For readout X)
		 29,28,12,11,10,  //      A2,B2,C2,D2,E2
		 27,26, 9, 8, 7,  //      A3,B3,C3,D3,E3
		 25,24, 6, 5, 4,  //      A4,B4,C4,D4,E4
		 23,22, 3, 2, 1,  //      A5,B5,C5,D5,E5
		 63,62,47,46,45,  //SLOT2 A1,B1,C1,D1,E1 (For readout Y)
		 61,60,44,43,42,  //      A2,B2,C2,D2,E2
		 59,58,41,40,39,  //      A3,B3,C3,D3,E3
		 57,56,38,37,36,  //      A4,B4,C4,D4,E4
		 55,54,35,34,33,  //      A5,B5,C5,D5,E5
		 31,30,15,14,13,  //SLOT3 A1,B1,C1,D1,E1 (For readout Z)
		 29,28,12,11,10,  //      A2,B2,C2,D2,E2
		 27,26, 9, 8, 7,  //      A3,B3,C3,D3,E3
		 25,24, 6, 5, 4,  //      A4,B4,C4,D4,E4
		 23,22, 3, 2, 1}; //      A5,B5,C5,D5,E5

  int Bij_YZ_z[25] = {5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1};

  int Bij_YZ_y[25] = {5, 5, 5, 5, 5,
                      4, 4, 4, 4, 4,
                      3, 3, 3, 3, 3,
                      2, 2, 2, 2, 2,
                      1, 1, 1, 1, 1};
  
  int Bij_ZX_x[25] = {5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1};

  int Bij_ZX_z[25] = {5, 5, 5, 5, 5,
                      4, 4, 4, 4, 4,
                      3, 3, 3, 3, 3,
                      2, 2, 2, 2, 2,
                      1, 1, 1, 1, 1};

  int Bij_XY_x[25] = {5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1,
                      5, 4, 3, 2, 1};

  int Bij_XY_y[25] = {5, 5, 5, 5, 5,
                      4, 4, 4, 4, 4,
                      3, 3, 3, 3, 3,
                      2, 2, 2, 2, 2,
                      1, 1, 1, 1, 1};
  
  Double_t binmin = -3.; // before calib : 700. / after calib : 0.
  Double_t binmax = 30.; // before calib : 2000. / after calib : 30
 
  TH2F* hHSX1 = new TH2F("hHSX1", "HSX1", 16, 0.5, 16.5, 1, 0.5, 1.5);
  hHSX1->SetMinimum(binmin);
  hHSX1->SetMaximum(binmax);
  TH2F* hHSY1 = new TH2F("hHSY1", "HSY1", 1, 0.5, 1.5, 16, 0.5, 16.5);
  hHSY1->SetMinimum(binmin);
  hHSY1->SetMaximum(binmax);
  TH2F* hHSY2 = new TH2F("hHSY2", "HSY2", 1, 0.5, 1.5, 16, 0.5, 16.5);
  hHSY2->SetMinimum(binmin);
  hHSY2->SetMaximum(binmax);
  TH2F* hHSX2 = new TH2F("hHSX2", "HSX2", 16, 0.5, 16.5, 1, 0.5, 1.5);
  hHSX2->SetMinimum(binmin);
  hHSX2->SetMaximum(binmax);
  TH2F* hProto_x = new TH2F("hProto_x", "Prototype : X axis", 5, 0.5, 5.5, 5, 0.5, 5.5);
  hProto_x->SetMinimum(binmin);
  hProto_x->SetMaximum(binmax+20);
  TH2F* hProto_y = new TH2F("hProto_y", "Prototype : Y axis", 5, 0.5, 5.5, 5, 0.5, 5.5);
  hProto_y->SetMinimum(binmin);
  hProto_y->SetMaximum(binmax+20);
  TH2F* hProto_z = new TH2F("hProto_z", "Prototype : Z axis", 5, 0.5, 5.5, 5, 0.5, 5.5);
  hProto_z->SetMinimum(binmin);
  hProto_z->SetMaximum(binmax+20);

  TH1F* hHitRate_HSX1 = new TH1F("hHitRate_HSX1", "HitRate_HSX1", 16, 0.5, 16.5);
  hHitRate_HSX1->SetMinimum(0);
  TH1F* hHitRate_HSY1 = new TH1F("hHitRate_HSY1", "HitRate_HSY1", 16, 0.5, 16.5);
  hHitRate_HSY1->SetMinimum(0);
  TH1F* hHitRate_HSY2 = new TH1F("hHitRate_HSY2", "HitRate_HSY2", 16, 0.5, 16.5);
  hHitRate_HSY2->SetMinimum(0);
  TH1F* hHitRate_HSX2 = new TH1F("hHitRate_HSX2", "HitRate_HSX2", 16, 0.5, 16.5);
  hHitRate_HSX2->SetMinimum(0);

  TH2F* hHitRate_XY = new TH2F("hHitRate_XY", "HitRate_XY", 5, 0.5, 5.5, 5, 0.5, 5.5);
  TH2F* hHitRate_ZX = new TH2F("hHitRate_ZX", "HitRate_ZX", 5, 0.5, 5.5, 5, 0.5, 5.5);
  TH2F* hHitRate_YZ = new TH2F("hHitRate_YZ", "HitRate_YZ", 5, 0.5, 5.5, 5, 0.5, 5.5);

  TH1F* hPE_X = new TH1F("hPE_X","PE_X", 130,-10,120);
  TH1F* hPE_Y = new TH1F("hPE_Y","PE_Y", 130,-10,120);
  TH1F* hPE_Z = new TH1F("hPE_Z","PE_Z", 130,-10,120);

  TH1F *hNHit_HSX1 = new TH1F("hNHit_HSX1","Number of Hits per Event (HSX1)",16,0.5,16.5);
  TH1F *hNHit_HSY1 = new TH1F("hNHit_HSY1","Number of Hits per Event (HSY1)",16,0.5,16.5);
  TH1F *hNHit_HSY2 = new TH1F("hNHit_HSY2","Number of Hits per Event (HSY2)",16,0.5,16.5);
  TH1F *hNHit_HSX2 = new TH1F("hNHit_HSX2","Number of Hits per Event (HSX2)",16,0.5,16.5);

  TCanvas*c0 = new TCanvas("c0","", 1800,600);
  
  TH1F* hPE[75] = {};
  TH1F* hPE_h[64] = {};
  Double_t array_PE[75] = {};
  //  Double_t array_PE_h[64] = {};

  for (int x=0; x<75; x++) {
      TString histname1 = TString::Format("Position No.%d ", x);
      // hPE[x] = new TH1F(histname1, histname1, 1800,700,2500);
      hPE[x] = new TH1F(histname1, histname1, 200,-5,95);
    }
  for (int x=0; x<64; x++) {
      TString histname2 = TString::Format("Hodoscope Ch%d ", x);
      // hPE_h[x] = new TH1F(histname2, histname2, 1800,700,2500);
      hPE_h[x] = new TH1F(histname2, histname2, 200,-5,95);
    }

  TH1F* hDist = new TH1F("hDist","",40,0.,20.);
  
  bool trig_HSX1 = false;
  bool trig_HSY1 = false;
  bool trig_HSY2 = false;
  bool trig_HSX2 = false;
  int all_count_x1 = 0;
  int all_count_y1 = 0;
  int all_count_y2 = 0;
  int all_count_x2 = 0;
  int trig_count = 0;
  int andtrig_count = 0;
  int nontrig_count = 0;
  int andnontrig_count = 0;
  
  bool trig_up_2s = false; 
  bool trig_down_2s = false;
  bool trig_1s = false;
  bool good_event_hs = false;
  bool DpassEvent_2s = false;
  bool DpassEvent_1s = false;
  
  string rootfile_dir = "../tree_root/";
  string calibfile_dir = "calibfile/";
  string evtdisplay_dir = "evtdisplay_dir/";

  Double_t mean_ped[3][64];
  Double_t mean_1pe[3][64];
  Double_t gain[3][64];
  Double_t sigma_ped[3][64];
  Double_t sigma_1pe[3][64];
  Double_t chi2_ndf_ped[3][64];
  Double_t chi2_ndf_1pe[3][64];
  
  Double_t pe_2s[64];
  Double_t pe_1s[64];
  Double_t pe_h[64];
  Int_t adc_2s[64];
  Int_t adc_1s[64];
  Int_t adc_h[64];

  int count_sum = 0;

  string pt2sfile;
  string pt1sfile;
  string hsfile;
  
  int ch;
  double mp, sp, cnp, m1, s1, cn1, g;

  for (int j=runnum-1; j<runnum; j++) {

    TString pt2sfile_name = TString::Format("%sproto2s_%04d_0002_tree2.root",
					    rootfile_dir.c_str(), j+1);
    TString pt1sfile_name = TString::Format("%sproto1s_%04d_0002_tree2.root",
					    rootfile_dir.c_str(), j+1);
    TString hsfile_name = TString::Format("%shodo_%04d_0002_tree2.root",
                                           rootfile_dir.c_str(), j+1);
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
    
    std::string pt2sfile_name_s = std::string(pt2sfile_name);
    std::string pt1sfile_name_s = std::string(pt1sfile_name);
    std::string hsfile_name_s = std::string(hsfile_name);
    
    string::size_type pos1 = pt2sfile_name_s.find(".root");
    string::size_type pos2 = pt1sfile_name_s.find(".root");
    string::size_type pos3 = hsfile_name_s.find(".root");
    string pt2sfile_calib(pt2sfile_name);
    pt2sfile_calib.replace(pos1,6,"_calib");
    pt2sfile_calib = pt2sfile_calib.substr(13);
    string pt1sfile_calib(pt1sfile_name);
    pt1sfile_calib.replace(pos2,6,"_calib");
    pt1sfile_calib = pt1sfile_calib.substr(13);
    string hsfile_calib(hsfile_name);
    hsfile_calib.replace(pos3,6,"_calib");
    hsfile_calib = hsfile_calib.substr(13);
    
    pt2sfile = pt2sfile_name_s;
    pt1sfile = pt1sfile_name_s;
    hsfile = hsfile_name_s;

  TString pt2s_calibname = TString::Format("%s%s.tsv",
					   calibfile_dir.c_str(), pt2sfile_calib.c_str());
  ifstream fin1(pt2s_calibname);
  while (fin1 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
    mean_ped[0][ch] = mp;
    sigma_ped[0][ch] = sp;
    chi2_ndf_ped[0][ch] = cnp;
    mean_1pe[0][ch] = m1;
    sigma_1pe[0][ch] = s1;
    chi2_ndf_1pe[0][ch] = cn1;
    gain[0][ch] = g;
  }
  
  TString pt1s_calibname = TString::Format("%s%s.tsv",
					   calibfile_dir.c_str(), pt1sfile_calib.c_str());
  ifstream fin2(pt1s_calibname);
  while (fin2 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
    mean_ped[1][ch] = mp;
    sigma_ped[1][ch] = sp;
    chi2_ndf_ped[1][ch] = cnp;
    mean_1pe[1][ch] = m1;
    sigma_1pe[1][ch] = s1;
    chi2_ndf_1pe[1][ch] = cn1;
    gain[1][ch] = g;
  }

  TString hs_calibname = TString::Format("%s%s.tsv",calibfile_dir.c_str(), hsfile_calib.c_str());
  ifstream fin3(hs_calibname);
  while (fin3 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
    mean_ped[2][ch] = mp;
    sigma_ped[2][ch] = sp;
    chi2_ndf_ped[2][ch] = cnp;
    mean_1pe[2][ch] = m1;
    sigma_1pe[2][ch] = s1;
    chi2_ndf_1pe[2][ch] = cn1;
    gain[2][ch] = g;
  }
  
    TChain *tree1 = new TChain("tree");
    TChain *tree2 = new TChain("tree");
    TChain *tree3 = new TChain("tree");
 
    tree1->Add(pt2sfile.c_str());
    tree1->SetBranchAddress("ADC", &adc_2s);
    tree2->Add(pt1sfile.c_str());
    tree2->SetBranchAddress("ADC", &adc_1s);
    tree3->Add(hsfile.c_str());
    tree3->SetBranchAddress("ADC", &adc_h);
   
    int omit_lowx1 = 0;
    int omit_highx1 = 0;
    int omit_lowy1 = 0;
    int omit_highy1 = 0;
    
    int omit_lowx2 = 0;
    int omit_highx2 = 0;
    int omit_lowy2 = 0;
    int omit_highy2 = 0;

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
    
    double dist_X = 0;
    double dist_Y = 0;
    double dist = 0;
    
    TString figname1;
     
 cout << "Protype+Hodoscope (Down): number of events = " << tree1->GetEntries() << endl;
 cout << "Hodoscope (Up): number of events = " << tree2->GetEntries() << endl;

 cout << "Select X1:" << "[" << omit_lowx1+1 << ", " <<16-omit_highx1 <<"]" << endl;
 cout << "Select Y1:" << "[" << omit_lowy1+1 << ", " <<16-omit_highy1 <<"]" << endl;
 cout << "Select Y2:" << "[" << omit_lowy2+1 << ", " <<16-omit_highy2 <<"]" << endl;
 cout << "Select X2:" << "[" << omit_lowx2+1 << ", " <<16-omit_highx2 <<"]" << endl;

  int less_evt = 0;
  int count = 0;
  int all_trig_num = 0;
  int DpassEvent_num[4] = {};
  int GoodEvent_num = 0;
  
  if (tree1->GetEntries() >= tree2->GetEntries()) less_evt = tree2->GetEntries();
  else less_evt = tree1->GetEntries();

  if (less_evt >= tree3->GetEntries()) less_evt = tree3->GetEntries();
  
  for(int evt = 0; evt < less_evt; evt++) {
    if (evt < gap_point_0) {
      tree1->GetEntry(evt);
      tree2->GetEntry(evt);
      tree3->GetEntry(evt);
    } else if (evt >= gap_point_0 && evt < gap_point_1) {
      tree1->GetEntry(evt-gap_pt2s_0);
      tree2->GetEntry(evt-gap_pt1s_0);
      tree3->GetEntry(evt-gap_hs_0);
    } else if (evt >= gap_point_1) {
      tree1->GetEntry(evt-gap_pt2s_1);
      tree2->GetEntry(evt-gap_pt1s_1);
      tree3->GetEntry(evt-gap_hs_1);
    }

     HitCount_HSX1 = 0;
     HitCount_HSY1 = 0;
     HitCount_HSY2 = 0;
     HitCount_HSX2 = 0;
 
    trig_HSX1 = false;
    trig_HSY1 = false;
    trig_HSY2 = false;
    trig_HSX2 = false;

    trig_up_2s = false;
    trig_down_2s = false;
    trig_1s = false;
    DpassEvent_2s = false;
    DpassEvent_1s = false;
    good_event_hs = false;
    
    max_HSX1 = 0;
    max_HSY1 = 0;
    max_HSY2 = 0;
    max_HSX2 = 0;
    maxCh_HSX1 = 0;
    maxCh_HSY1 = 16;
    maxCh_HSY2 = 32;
    maxCh_HSX2 = 48;

    dist_X = 0;
    dist_Y = 0;
    dist = 0;

    figname1 = TString::Format("cosmic_test/evtdisplay/hit_%d_%d.pdf",runnum,evt);
    
    for ( int i=0; i<64; i++ ) {
      //  pe_h[i] = (double(adc_h[i]));
      //      cout << gain[1][i] << endl;
      pe_h[i] = (double(adc_h[i]) - mean_ped[2][i]) / gain[2][i];
      pe_h[i] = double(floor(pe_h[i]*10)/10);
      
      if ( i <= 15 ) {
	hHSX1->SetBinContent(16-i,1,pe_h[i]);
	if (pe_h[i] > 7.5
	    && i >=(omit_lowx1+0) && i<=(15-omit_highx1)) {
	  hHitRate_HSX1->Fill(16-i);
	  HitCount_HSX1++;
	  all_count_x1++;
	  trig_HSX1 = true;
	  max_HSX1 = pe_h[0];
	  if (pe_h[i] > max_HSX1) maxCh_HSX1 = 15-i;
	}
      }
      else if ( i>=16 && i <= 31 ) {
	hHSY1->SetBinContent(1,32-i,pe_h[i]);
	if (pe_h[i] > 7.5
	    && i >=(16+omit_lowy1) && i <= (31-omit_highy1)) {
	  hHitRate_HSY1->Fill(32-i);
	  HitCount_HSY1++;
	  all_count_y1++;
	  trig_HSY1 = true;
	  max_HSY1 = pe_h[16];
	  if (pe_h[i] > max_HSY1) maxCh_HSY1 = 31-i;
	}
      }
      else if ( i >=32 && i <= 47 ) {
	hHSY2->SetBinContent(1,i-31,pe_h[i]);
	if (pe_h[i] > 7.5
	    && i >=(32+omit_lowy2) && i <= (47-omit_highy2)) {
	   hHitRate_HSY2->Fill(i-31);
	   HitCount_HSY2++;
	   all_count_y2++;
	  trig_HSY2 = true;
	  max_HSY2 = pe_h[32];
	  if (pe_h[i] > max_HSY2) maxCh_HSY2 = i-32;
	}
      }
      else if ( i >= 48 && i <= 63 ) {
	hHSX2->SetBinContent(i-47,1,pe_h[i]);
	if (pe_h[i] > 7.5
	    && i >= (48+omit_lowx2) && i <= (63-omit_highx2)) {
	  hHitRate_HSX2->Fill(i-47);
	  HitCount_HSX2++;
	  all_count_x2++;
	  trig_HSX2 = true;
	  max_HSX2 = pe_h[48];
	  if (pe_h[i] > max_HSX2) maxCh_HSX2 = i-48;
	}
      }
      
    }

    hNHit_HSX1->Fill(HitCount_HSX1);
    hNHit_HSY1->Fill(HitCount_HSY1);
    hNHit_HSY2->Fill(HitCount_HSY2);
    hNHit_HSX2->Fill(HitCount_HSX2);
    
    if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2) good_event_hs = true;

    
    for ( int i = 0; i < 75; i++ ) {

      // pe_t[map[i]] = (double(adc_t[map[i]]));
      if (i < 50) {
	pe_2s[map[i]] = (double(adc_2s[map[i]]) - mean_ped[0][map[i]]) / gain[0][map[i]];
	pe_2s[map[i]] = double(floor(pe_2s[map[i]]*10)/10);
	
        if (i < 25) {
	  hProto_x->SetBinContent(i%5+1, i/5+1, pe_2s[map[i]]);
	  if ( i%5 == 0 && pe_2s[map[i]] > 9.5 ) {
	    trig_up_2s = true;
	  }
	  if ( i%5 == 4 && pe_2s[map[i]] > 9.5 ) {
	    trig_down_2s = true;
	  }	  
	  if (pe_2s[map[i]] > 9.5) {
	    hHitRate_YZ->Fill(Bij_YZ_z[i], Bij_YZ_y[i]);
	    if (good_event_hs) {
	      hPE_X->Fill(pe_2s[map[i]]);
	    }
	  } 
	} else if (i >= 25 && i < 50) {
	  hProto_y->SetBinContent(i%5+1, i/5-4, pe_2s[map[i]]);
	  if (pe_2s[map[i]] > 9.5) {
	    hHitRate_ZX->Fill(Bij_ZX_x[i-25], Bij_ZX_z[i-25]);
	    if (good_event_hs) {
	      hPE_Y->Fill(pe_2s[map[i]]);
	    }
	  }
	}

      } else if (i >= 50) {
	pe_1s[map[i]] = (double(adc_1s[map[i]]) - mean_ped[1][map[i]]) / gain[1][map[i]];
	pe_1s[map[i]] = double(floor(pe_1s[map[i]]*10)/10);
	hProto_z->SetBinContent(i%5+1, i/5-9, pe_1s[map[i]]);
	
	if (pe_1s[map[i]] > 9.5) {
	  hHitRate_XY->Fill(Bij_XY_x[i-50], Bij_XY_y[i-50]);
	  if (good_event_hs) {
	    hPE_Z->Fill(pe_1s[map[i]]);
	  }
	}

	if ( (i ==  6 || i ==  7 || i ==  8 ||
	      i == 11 || i == 12 || i == 13 ||
	      i == 16 || i == 17 || 1 == 18 ) && pe_1s[map[i]] > 29.5 ) {
	  trig_1s = true;
	}
      }
    }

    if (trig_up_2s && trig_down_2s) {
      DpassEvent_2s = true;
      // cout << "D-passing event# : " << evt << endl;
      if (DpassEvent_num[0] == 0) DpassEvent_num[0] = evt;
      else {
	DpassEvent_num[1] = DpassEvent_num[0];
	DpassEvent_num[0] = evt;
      }
      count++;
      if (DpassEvent_num[1]-GoodEvent_num < 0) {
	cout << "Diff between PT & HS (bigger:2s): " << DpassEvent_num[0] - GoodEvent_num << endl;
      }
    }

    
    if (trig_1s) {
      DpassEvent_1s = true;
      // cout << "D-passing event# : " << evt << endl;
      if (DpassEvent_num[2] == 0) DpassEvent_num[0] = evt;
      else {
	DpassEvent_num[3] = DpassEvent_num[2];
	DpassEvent_num[2] = evt;
      }
      count++;
      if (DpassEvent_num[3]-GoodEvent_num < 0) {
	cout << "Diff between PT & HS (bigger:1s): " << DpassEvent_num[2] - GoodEvent_num << endl;
      }
    }
 
    for ( int i=0; i<75; i++ ) {
      if (DpassEvent_2s && DpassEvent_1s) {
	if (i < 50) hPE[i]->Fill(pe_2s[map[i]]);
	else hPE[i]->Fill(pe_1s[map[i]]);
      }
    }

    if (trig_HSX1 || trig_HSY1) trig_count++;
    if (trig_HSX1 && trig_HSY1) andtrig_count++;
    if (trig_HSX2 || trig_HSY2) nontrig_count++;
    if (trig_HSX2 && trig_HSY2) andnontrig_count++;
    
    if (good_event_hs) {
      // cout << "GOOD event# : " << evt << endl;
      GoodEvent_num = evt;
      cout << "Diff between PT & HS (smaller): " << DpassEvent_num[0] - GoodEvent_num
	   << " / " << DpassEvent_num[2] - GoodEvent_num << endl;
      all_trig_num++;

      dist_X = abs(maxCh_HSX2-maxCh_HSX1);
      dist_Y = abs(maxCh_HSY2-maxCh_HSY1);
      dist = sqrt(pow(dist_X,2) + pow(dist_Y,2));

      hDist->Fill(dist);

      if (evt >= min_evt && evt <= max_evt) {

	gStyle->SetOptStat(0);
	c0->Clear();
	c0->Divide(4,2);

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
	hProto_x->GetYaxis()->SetTitle();
	hProto_x->GetXaxis()->SetTitle();
	hProto_x->GetYaxis()->SetTitleSize(0.03);
	hProto_x->GetXaxis()->SetTitleSize(0.03);
	hProto_x->GetYaxis()->SetLabelSize(0.04);
	hProto_x->GetXaxis()->SetLabelSize(0.04);
	hProto_x->GetXaxis()->SetNdivisions(5);
	hProto_x->GetYaxis()->SetNdivisions(5);
	hProto_x->SetMarkerSize(2.0);
	hProto_x->Draw("text colz");

	c0->cd(7);
	hProto_y->GetYaxis()->SetTitle();
	hProto_y->GetXaxis()->SetTitle();
	hProto_y->GetYaxis()->SetTitleSize(0.03);
	hProto_y->GetXaxis()->SetTitleSize(0.03);
	hProto_y->GetYaxis()->SetLabelSize(0.04);
	hProto_y->GetXaxis()->SetLabelSize(0.04);
	hProto_y->GetXaxis()->SetNdivisions(5);
	hProto_y->GetYaxis()->SetNdivisions(5);
	hProto_y->SetMarkerSize(2.0);
	hProto_y->Draw("text colz");
	
	c0->cd(3);
	hProto_z->GetYaxis()->SetTitle();
	hProto_z->GetXaxis()->SetTitle();
	hProto_z->GetYaxis()->SetTitleSize(0.03);
	hProto_z->GetXaxis()->SetTitleSize(0.03);
	hProto_z->GetYaxis()->SetLabelSize(0.04);
	hProto_z->GetXaxis()->SetLabelSize(0.04);
	hProto_z->GetXaxis()->SetNdivisions(5);
	hProto_z->GetYaxis()->SetNdivisions(5);
	hProto_z->SetMarkerSize(2.0);
	hProto_z->Draw("text colz");
	
	c0->SaveAs(figname1);
      }
    }
  }


  TString figname3 = TString::Format("prototype/distance_%d.pdf",runnum);
  TCanvas*c1 = new TCanvas("c1","",800,600);
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

  cout << all_count_x1  << " events passed X1" <<endl;
  cout << all_count_y1  << " events passed Y1" <<endl;
  cout << all_count_y2  << " events passed Y2" <<endl;
  cout << all_count_x2  << " events passed X2" <<endl;
  
  cout << "DAQ trigger counts : " << trig_count << endl;
  cout << "1st layer counts : " << HitCount_HSX1 << endl;
  cout << "non-trigger hodoscope OR counts : " << nontrig_count << endl;
  cout << "trigger hodoscope AND counts : " << andtrig_count << endl;
  cout << "non-trigger hodoscope AND counts : " << andnontrig_count << endl;
  count_sum += count;

  cout << count << " good events out of " << tree1->GetEntries() << " successfully processed" << endl;
  cout << "events which passed through all hodoscpe layers : " << all_trig_num << endl;
  //  cout << "(in total : " << count_sum << " good events out of " << max_evt_sum << " )" << endl;

  }

ofstream outputtext5("cosmic_test/75ch/array_PE.tsv");
TCanvas*c2 = new TCanvas("c2","",800,600);
TString pdfname1;

  for (int x=0; x<75; x++) {
    f1->SetParameters(10, 60, 1.);
    hPE[x]->Fit("landau_orig", "", "", 30, 75);
    array_PE[x] = f1->GetParameter(1);
    pdfname1 = TString::Format("cosmic_test/75ch/run%d_PE_position%d.pdf", runnum, x);    

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

    TCanvas*c3 = new TCanvas("c3","",1800,600);
  c3->Clear();
  c3->Divide(4,2);

  c3->cd(3);
  hHitRate_XY->GetXaxis()->SetTitle("Cube# along X");
  hHitRate_XY->GetYaxis()->SetTitle("Cube# along Y");
  hHitRate_XY->GetYaxis()->SetTitleSize(0.03);
  hHitRate_XY->GetXaxis()->SetTitleSize(0.03);
  hHitRate_XY->GetYaxis()->SetLabelSize(0.04);
  hHitRate_XY->GetXaxis()->SetLabelSize(0.04);
  hHitRate_XY->GetXaxis()->SetNdivisions(5);
  hHitRate_XY->GetYaxis()->SetNdivisions(5);
  hHitRate_XY->SetMarkerSize(2.0);
  hHitRate_XY->Draw("text colz");

  c3->cd(7);
  hHitRate_ZX->GetXaxis()->SetTitle("Cube# along X");
  hHitRate_ZX->GetYaxis()->SetTitle("Cube# along Z");
  hHitRate_ZX->GetYaxis()->SetTitleSize(0.03);
  hHitRate_ZX->GetXaxis()->SetTitleSize(0.03);
  hHitRate_ZX->GetYaxis()->SetLabelSize(0.04);
  hHitRate_ZX->GetXaxis()->SetLabelSize(0.04);
  hHitRate_ZX->GetXaxis()->SetNdivisions(5);
  hHitRate_ZX->GetYaxis()->SetNdivisions(5);
  hHitRate_ZX->SetMarkerSize(2.0);
  hHitRate_ZX->Draw("text colz");

  c3->cd(2);
  hHitRate_YZ->GetYaxis()->SetTitle("Cube# along Y");
  hHitRate_YZ->GetXaxis()->SetTitle("Cube# along Z");
  hHitRate_YZ->GetYaxis()->SetTitleSize(0.03);
  hHitRate_YZ->GetXaxis()->SetTitleSize(0.03);
  hHitRate_YZ->GetYaxis()->SetLabelSize(0.04);
  hHitRate_YZ->GetXaxis()->SetLabelSize(0.04);
  hHitRate_YZ->GetXaxis()->SetNdivisions(5);
  hHitRate_YZ->GetYaxis()->SetNdivisions(5);
  hHitRate_YZ->SetMarkerSize(2.0);
  hHitRate_YZ->Draw("text colz");


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

  c3->SaveAs("prototype/hitRate_allRun.pdf");


  TCanvas*c4 = new TCanvas("c4","",800,600);
  c4->Clear();
  c4->Divide(2,2);
  
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

  c4->SaveAs("prototype/hodoscope_HitMultiplicity_allRun.pdf");

  TCanvas*c5 = new TCanvas("c5","",800,600);
  c5->Clear();
  c5->Divide(2,2);

  c5->cd(1);
  hPE_X->Fit("landau","","",25,50);
  hPE_X->GetYaxis()->SetTitle("Number of Events");
  hPE_X->GetXaxis()->SetTitle("p.e.");
  hPE_X->GetYaxis()->SetTitleSize(0.04);
  hPE_X->GetXaxis()->SetTitleSize(0.04);
  hPE_X->GetYaxis()->SetLabelSize(0.05);
  hPE_X->GetXaxis()->SetLabelSize(0.05);
  hPE_X->Draw();

  c5->cd(4);
  hPE_Y->Fit("landau","","",25,50);
  hPE_Y->GetYaxis()->SetTitle("Number of Events");
  hPE_Y->GetXaxis()->SetTitle("p.e.");
  hPE_Y->GetYaxis()->SetTitleSize(0.04);
  hPE_Y->GetXaxis()->SetTitleSize(0.04);
  hPE_Y->GetYaxis()->SetLabelSize(0.05);
  hPE_Y->GetXaxis()->SetLabelSize(0.05);
  hPE_Y->Draw();

  c5->cd(2);
  hPE_Z->Fit("landau","","",25,50);
  hPE_Z->GetYaxis()->SetTitle("Number of Events");
  hPE_Z->GetXaxis()->SetTitle("p.e.");
  hPE_Z->GetYaxis()->SetTitleSize(0.04);
  hPE_Z->GetXaxis()->SetTitleSize(0.04);
  hPE_Z->GetYaxis()->SetLabelSize(0.05);
  hPE_Z->GetXaxis()->SetLabelSize(0.05);
  hPE_Z->Draw();

  c5->SaveAs("prototype/proto_PE_allRun.pdf");
  
}

int main(int argc, char** argv)
{
  if (argc != 12) {
    cerr << "run_proto <runnum> <min_evt> <max_evt> <gap_point_0 & gap> <gap_point_1 & gap>"
	 << endl;
    return -1;
  }

  int runnum = atoi(argv[1]);
  int min_evt = atoi(argv[2]);
  int max_evt = atoi(argv[3]);
  int gap_point_0 = atoi(argv[4]);
  int gap_pt2s_0 = atoi(argv[5]);
  int gap_pt1s_0 = atoi(argv[6]);
  int gap_hs_0 = atoi(argv[7]);
  int gap_point_1 = atoi(argv[8]);
  int gap_pt2s_1 = atoi(argv[9]);
  int gap_pt1s_1 = atoi(argv[10]);
  int gap_hs_1 = atoi(argv[11]);

run_proto(runnum, min_evt, max_evt,
	  gap_point_0, gap_pt2s_0, gap_pt1s_0, gap_hs_0,
	  gap_point_1, gap_pt2s_1, gap_pt1s_1, gap_hs_1);
return 0;
}
