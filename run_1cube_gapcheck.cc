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

void run_1cube(int &runnum, int &subnum ,/*string &dir*/ int &file_count, int &min_evt, int &max_evt, int &loop) {

  gStyle->SetOptStat(1111);
  //  gStyle->SetOptFit(1111);

  TF1 *f1 = new TF1("landau_orig","[0]*TMath::Landau(x,[1],[2])");
  TF1 *f2 = new TF1("gaus_orig","[0]*TMath::Gaus(x,[1],[2])");
  
  int map[3] = {9,41,9 }; //

  Double_t binmin = -3.; // before calib : 700. / after calib : 0.
  Double_t binmax = 30.; // before calib : 2000. / after calib : 30

  
  int gapcheck_hsx1_ch=0;
  int gapcheck_hsy1_ch=0;
  int gapcheck_hsx2_ch=0;
  int gapcheck_hsy2_ch=0;

  

  int gap_count = 0;
  int gap_cons_count = 0;
  int gap_cons_count_l = 0;
  int gap_count_x = 0;
  int gap_cons_count_x =0;
  int gap_cons_count_l_x = 0;
	 
  TH1F* hPE_cell_all[16][16] = {};
  TH1F* hPE_cell_x[16][16] = {};
  TH1F* hPE_cell_y[16][16] = {};
  TH1F* hPE_cell_z[16][16] = {};
  Double_t cell_PE_all[16][16] = {};
  Double_t cell_PE_x[16][16] = {};
  Double_t cell_PE_y[16][16] = {};
  Double_t cell_PE_z[16][16] = {};
  Double_t eff_all[16][16] = {};
  Double_t eff_x[16][16] = {};
  Double_t eff_y[16][16] = {};
  Double_t eff_z[16][16] = {};
  Double_t allpeak[16][16];
  Double_t xpeak[16][16];
  Double_t ypeak[16][16];
  Double_t zpeak[16][16];
  
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
  TH2F* hCubeX = new TH2F("hCubeX", "1cube : X readout (event-by-event)", 1, 0.5, 1.5, 1, 0.5, 1.5);
  hCubeX->SetMinimum(binmin);
  hCubeX->SetMaximum(binmax+20);
  TH2F* hCubeY = new TH2F("hCubeY", "1cube : Y readout (event-by-event)", 1, 0.5, 1.5, 1, 0.5, 1.5);
  hCubeY->SetMinimum(binmin);
  hCubeY->SetMaximum(binmax+20);
  TH2F* hCubeZ = new TH2F("hCubeZ", "1cube : Z readout (event-by-event)", 1, 0.5, 1.5, 1, 0.5, 1.5);
  hCubeZ->SetMinimum(binmin);
  hCubeZ->SetMaximum(binmax+20);

  TH1F* hHitRate_HSX1 = new TH1F("hHitRate_HSX1", "HitRate_HSX1", 16, 0.5, 16.5);
  hHitRate_HSX1->SetMinimum(0);
  TH1F* hHitRate_HSY1 = new TH1F("hHitRate_HSY1", "HitRate_HSY1", 16, 0.5, 16.5);
  hHitRate_HSY1->SetMinimum(0);
  TH1F* hHitRate_HSY2 = new TH1F("hHitRate_HSY2", "HitRate_HSY2", 16, 0.5, 16.5);
  hHitRate_HSY2->SetMinimum(0);
  TH1F* hHitRate_HSX2 = new TH1F("hHitRate_HSX2", "HitRate_HSX2", 16, 0.5, 16.5);
  hHitRate_HSX2->SetMinimum(0);

  TH2F* hHitRate_X = new TH2F("hHitRate_X", "HitRate_X", 1, 0.5, 1.5, 1, 0.5, 1.5);
  TH2F* hHitRate_Y = new TH2F("hHitRate_Y", "HitRate_Y", 1, 0.5, 1.5, 1, 0.5, 1.5);
  TH2F* hHitRate_Z = new TH2F("hHitRate_Z", "HitRate_Z", 1, 0.5, 1.5, 1, 0.5, 1.5);

  const Int_t GapGraphMax = 50000;
  //  TString gap_cons_count_h_name = TString::Format("gap_cons_count_z_%D+%D_%D",runnum, file_count -1 ,nasub);
  
  //TH1D *gap_cons_count_h = new TH1D(gap_cons_count_h_name.Data(),gap_cons_count_h_name,Data(),GapGraphMax,0,GapGraphMax);
  
  //add gapcheck 18/11/15 11:47

  TH1D* gap_count_h_1s2s[100];
  TH1D* gap_count_h_hodo1s[100];
  TH2F* hHSC_gapcheck_2D[100];
  TString gap_count_h_1s2s_name;
  TString gap_count_h_hodo1s_name;
  TString hHSC_gapcheck_2D_name;
  for(int i=0;i<file_count;i++){
    gap_count_h_1s2s_name="";
    gap_count_h_1s2s_name = TString::Format("gap_count_1s2s_%d_%d",runnum+i ,subnum);
    gap_count_h_1s2s[i] = new TH1D(gap_count_h_1s2s_name,gap_count_h_1s2s_name,GapGraphMax,0,GapGraphMax);
    gap_count_h_hodo1s_name="";
    gap_count_h_hodo1s_name = TString::Format("gap_count_hodo1s_%d_%d",runnum+i ,subnum);
    gap_count_h_hodo1s[i] = new TH1D(gap_count_h_hodo1s_name,gap_count_h_hodo1s_name,GapGraphMax,0,GapGraphMax);
    hHSC_gapcheck_2D_name="";
    hHSC_gapcheck_2D_name = TString::Format("HSC_gapcheck_%d_%d",runnum+i,subnum);
    hHSC_gapcheck_2D[i] = new TH2F(hHSC_gapcheck_2D_name,hHSC_gapcheck_2D_name,16,0.5,16.5,16,0.5,16.5);
  }
  //
  
  TH2F *hHSCProf2D[16][16];
  TString name;

  for (int nx=0;nx<16; nx++){
    for (int ny=0;ny<16; ny++){
      name = TString::Format("HSCProf2D_%d_%d",nx,ny);
      hHSCProf2D[nx][ny] = new TH2F(name,name,16,0.5,16.5,16,0.5,16.5);
    }
  }

  TH2F *hPEdist2D[16][16];
  TString pname;

  for(int nx =0;nx<16; nx++){
    for (int ny =0;ny<16;ny++){
      pname = TString::Format("PEdist2D_%d_%d",nx,ny);
      hPEdist2D[nx][ny] = new TH2F(pname,pname,16,0.5,16.5,16,0.5,16.5);
    }
  }
  
  TH1F* hPE_X = new TH1F("hPE_X","PE_X", 130,-10,120);
  TH1F* hPE_Y = new TH1F("hPE_Y","PE_Y", 130,-10,120);
  TH1F* hPE_Z = new TH1F("hPE_Z","PE_Z", 130,-10,120);

  TH1F *hNHit_HSX1 = new TH1F("hNHit_HSX1","Number of Hits per Event (HSX1)",16,0.5,16.5);
  TH1F *hNHit_HSY1 = new TH1F("hNHit_HSY1","Number of Hits per Event (HSY1)",16,0.5,16.5);
  TH1F *hNHit_HSY2 = new TH1F("hNHit_HSY2","Number of Hits per Event (HSY2)",16,0.5,16.5);
  TH1F *hNHit_HSX2 = new TH1F("hNHit_HSX2","Number of Hits per Event (HSX2)",16,0.5,16.5);

  for (int x=0; x<16; x++) {
    for (int y=0; y<16; y++) {
      TString cellhistname1 = TString::Format("cell(all) x:%d y:%d", x, y);
      TString cellhistname2 = TString::Format("cell(x) x:%d y:%d", x, y);
      TString cellhistname3 = TString::Format("cell(y) x:%d y:%d", x, y);
      TString cellhistname4 = TString::Format("cell(z) x:%d y:%d", x, y);
      hPE_cell_all[x][y] = new TH1F(cellhistname1, cellhistname1, 60,-10,80);
      hPE_cell_x[x][y] = new TH1F(cellhistname2, cellhistname2, 45,-10,80);
      hPE_cell_y[x][y] = new TH1F(cellhistname3, cellhistname3, 45,-10,80);
      hPE_cell_z[x][y] = new TH1F(cellhistname4, cellhistname4, 45,-10,80);
    }
  }

  TCanvas*c0 = new TCanvas("c0","", 1800,600);
  
  TH1F* hPE[25] = {};
  TH1F* hPE_h[64] = {};
  //  Double_t array_PE_h[64] = {};

  for (int x=0; x<25; x++) {
      TString histname1 = TString::Format("Position No.%d ", x);
      // hPE[x] = new TH1F(histname1, histname1, 1800,700,2500);
      hPE[x] = new TH1F(histname1, histname1, 200,-5,95);
    }
  for (int x=0; x<64; x++) {
      TString histname2 = TString::Format("Hodoscope Ch%d ", x);
      // hPE_h[x] = new TH1F(histname2, histname2, 1800,700,2500);
      hPE_h[x] = new TH1F(histname2, histname2, 200,-5,95);
    }

  //  TH1F* hDist = new TH1F("hDist","",40,0.,20.);
  
  bool trig_HSX1 = false;
  bool trig_HSY1 = false;
  bool trig_HSY2 = false;
  bool trig_HSX2 = false;
  bool trig_X = false;
  bool trig_Y = false;
  bool trig_Z = false;
  int all_count_x1 = 0;
  int all_count_y1 = 0;
  int all_count_y2 = 0;
  int all_count_x2 = 0;
  int trig_count = 0;
  int andtrig_count = 0;
  int nontrig_count = 0;
  int andnontrig_count = 0;
  
  bool good_event = false;
  bool DpassEvent = false;
  bool SingleHit = false;
  
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
  
  Double_t pe_p1[64];
  Double_t pe_p2[64];
  Double_t pe_h[64];
  Int_t adc_p1[64];
  Int_t adc_p2[64];
  Int_t adc_h[64];

  Int_t cell_x[10][10] = {};
  Int_t cell_y[10][10] = {};
  Int_t cell_z[10][10] = {};

  int less_evt_sum = 0;
  int good_count_sum = 0;
  int Dpass_count_sum = 0;

  string pt1sfile;
  string pt2sfile;
  string hsfile;

  int loop_num = 1;
  
  int ch;
  double mp, sp, cnp, m1, s1, cn1, g;

  for (int j=runnum; j<runnum+file_count; j++) {

    //add gap_count 18/11/15 11:47
    int gap_count_1s2s = 0;
    int gap_count_hodo1s=0;
    //add gapcheck 18/11/15 15:10
    int Gcheck_hsx1_ch=0;
    int Gcheck_hsy1_ch=0;
    int Gcheck_hsx2_ch=0;
    int Gcheck_hsy2_ch=0;
    //

    TString pt1sfile_name = TString::Format("%sproto1s_%04d_000%d_tree2.root",rootfile_dir.c_str(), j,subnum);

    //for test
    //TString pt1sfile_name = TString::Format("%sproto1s_0005_000%d_tree2.root",rootfile_dir.c_str(), subnum);

    TString pt2sfile_name = TString::Format("%sproto2s_%04d_000%d_tree2.root", rootfile_dir.c_str(),j,subnum);  

    //for test
    //TString pt2sfile_name = TString::Format("%sproto2s_0005_000%d_tree2.root", rootfile_dir.c_str(),subnum);  

    TString hsfile_name = TString::Format("%shodo_%04d_000%d_tree2.root", rootfile_dir.c_str(), j,subnum);

    //for test
    //TString hsfile_name = TString::Format("%shodo_0005_000%d_tree2.root", rootfile_dir.c_str(), subnum);

    /*
      int gap_point_0 = 0;
      int gap_pt_0 = 0;
      int gap_hs_0 = 0;
      int gap_point_1 = 0;
      int gap_pt_1 = 0;
      int gap_hs_1 = 0;

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
    
    std::string pt1sfile_name_s = std::string(pt1sfile_name);
    std::string pt2sfile_name_s = std::string(pt2sfile_name);
    std::string hsfile_name_s = std::string(hsfile_name);
    
    //string::size_type pos1 = pt1sfile_name_s.find(".root");
    //string::size_type pos2 = pt2sfile_name_s.find(".root");
    //string::size_type pos3 = hsfile_name_s.find(".root");
    /*string pt1sfile_calib(pt1sfile_name);
      pt1sfile_calib.replace(pos1,6,"_calib");
      pt1sfile_calib = pt1sfile_calib.substr(13);
      string pt2sfile_calib(pt2sfile_name);
      pt2sfile_calib.replace(pos2,6,"_calib");
      pt2sfile_calib = pt2sfile_calib.substr(13);
      string hsfile_calib(hsfile_name);
      hsfile_calib.replace(pos3,6,"_calib");
      hsfile_calib = hsfile_calib.substr(13);*/

    string pt1sfile_calib = "proto1s_0100_0002_tree2_calib";
    string pt2sfile_calib ="proto2s_0100_0002_tree2_calib";
    string hsfile_calib ="hodo_0100_0002_tree2_calib";
    
    pt1sfile = pt1sfile_name_s;
    pt2sfile = pt2sfile_name_s;
    hsfile = hsfile_name_s;
    TString pt1s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt1sfile_calib.c_str());
    ifstream fin1(pt1s_calibname);
    while (fin1 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
      mean_ped[0][ch] = mp;
      sigma_ped[0][ch] = sp;
      chi2_ndf_ped[0][ch] = cnp;
      mean_1pe[0][ch] = m1;
      sigma_1pe[0][ch] = s1;
      chi2_ndf_1pe[0][ch] = cn1;
      gain[0][ch] = g;
    }
    TString pt2s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt2sfile_calib.c_str());
    ifstream fin2(pt2s_calibname);
    while (fin2 >> ch >> mp >> sp >> cnp >> m1 >> s1 >> cn1 >> g) {
      mean_ped[1][ch] = mp;
      sigma_ped[1][ch] = sp;
      chi2_ndf_ped[1][ch] = cnp;
      mean_1pe[1][ch] = m1;
      sigma_1pe[1][ch] = s1;
      chi2_ndf_1pe[1][ch] = cn1;
      gain[1][ch] = g;
    }
  
    TString hs_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), hsfile_calib.c_str());
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
  
    tree1->Add(pt1sfile.c_str());
    tree1->SetBranchAddress("ADC", &adc_p1);
    tree2->Add(pt2sfile.c_str());
    tree2->SetBranchAddress("ADC", &adc_p2);
    tree3->Add(hsfile.c_str());
    tree3->SetBranchAddress("ADC", &adc_h);
  
  
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
  
    if (loop == 1) loop_num = 16;
  
    for (int y=0; y<loop_num; y++) {
      for (int x=0; x<loop_num; x++) {
	int omit_lowx1 = 0;
	int omit_highx1 = 0;
	int omit_lowy1 = 0;
	int omit_highy1 = 0;
      
	int omit_lowx2 = 0;
	int omit_highx2 = 0;
	int omit_lowy2 = 0;
	int omit_highy2 = 0;
      
	if (loop_num == 16) {
	  omit_lowx1 = 15-x;
	  omit_highx1 = x;
	  omit_lowy1 = 15-y;
	  omit_highy1 = y;
	
	  omit_lowx2 = x;
	  omit_highx2 = 15-x;
	  omit_lowy2 = y;
	  omit_highy2 = 15-y;
	}
      
	cout << "Protype   : number of events = " << tree1->GetEntries() << endl;
	cout << "Hodoscope : number of events = " << tree3->GetEntries() << endl;
      
	cout << "Select X1:" << "[" << omit_lowx1+1 << ", " <<16-omit_highx1 <<"]" << endl;
	cout << "Select Y1:" << "[" << omit_lowy1+1 << ", " <<16-omit_highy1 <<"]" << endl;
	cout << "Select Y2:" << "[" << omit_lowy2+1 << ", " <<16-omit_highy2 <<"]" << endl;
	cout << "Select X2:" << "[" << omit_lowx2+1 << ", " <<16-omit_highx2 <<"]" << endl;
      
	int less_evt = 0;
	int good_count = 0;
	int Dpass_count = 0;
	int DpassEvent_num[2] = {};
	int GoodEvent_num = 0;
      
	if (tree1->GetEntries() >= tree3->GetEntries()) less_evt = tree3->GetEntries();
	else less_evt = tree1->GetEntries();
      
	for(int evt = 0; evt < less_evt; evt++) {
	  tree1->GetEntry(evt);
	  tree2->GetEntry(evt);
	  tree3->GetEntry(evt);
	  /*
	    if (evt < gap_point_0) {
	    tree1->GetEntry(evt);
	    tree2->GetEntry(evt);
	    } else if (evt >= gap_point_0 && evt < gap_point_1) {
	    tree1->GetEntry(evt-gap_pt_0);
	    tree2->GetEntry(evt-gap_hs_0);
	    } else if (evt >= gap_point_1) {
	    tree1->GetEntry(evt-gap_pt_1);
	    tree2->GetEntry(evt-gap_hs_1);
	    }
	  */
	  HitCount_HSX1 = 0;
	  HitCount_HSY1 = 0;
	  HitCount_HSY2 = 0;
	  HitCount_HSX2 = 0;
	
	  trig_HSX1 = false;
	  trig_HSY1 = false;
	  trig_HSY2 = false;
	  trig_HSX2 = false;
	
	  trig_X = false;
	  trig_Y = false;
	  trig_Z = false;
	  DpassEvent = false;
	  good_event = false;
	  SingleHit = false;
	
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
	
	
	  figname1 = TString::Format("1cube/evtdisplay/hit_%d_%d_%d.pdf",runnum,subnum,evt);
	
	  //HODOSCOPE
	  const double hodoThre = 2.5;
	
	  double max_hsx1=0;
	  double max_hsy1=0;
	  double max_hsx2=0;
	  double max_hsy2=0;
	
	  for ( int i=0; i<64; i++ ) {
	    //  pe_h[i] = (double(adc_h[i]));
	    pe_h[i] = (double(adc_h[i]) - mean_ped[2][i]) / gain[2][i];
	    pe_h[i] = double(floor(pe_h[i]*10)/10);
	  
	    //edit gapcheck_hs 18/11/15 12:02
	    if (pe_h[i] > hodoThre && i >=(omit_lowx1+0) && i<=(15-omit_highx1)) {
	      if ( i <= 15 ) {
		if(max_hsx1<pe_h[i]){
		  max_hsx1=pe_h[i];
		  Gcheck_hsx1_ch=i;
		}
	      }
	    }
	    if (pe_h[i] > hodoThre && i >=(16+omit_lowy1) && i <= (31-omit_highy1)) {
	      if ( i>=16 && i <= 31 ) {
		if(max_hsy1<pe_h[i]){
		  max_hsy1=pe_h[i];
		  Gcheck_hsy1_ch=i;
		}
	      }
	    }
	    if (pe_h[i] > hodoThre && i >=(32+omit_lowy2) && i <= (47-omit_highy2)) {
	      if ( i >=32 && i <= 47 ) {
		if(max_hsy2<pe_h[i]){
		  max_hsy2=pe_h[i];
		  Gcheck_hsy2_ch =i;
		}
	      }
	    }
	    if (pe_h[i] > hodoThre && i >= (48+omit_lowx2) && i <= (63-omit_highx2)) {
	      if(i>=48){
		if(max_hsx2<pe_h[i]){
		  max_hsx2=pe_h[i];
		  Gcheck_hsx2_ch=i;
		}
	      }
	    }
	    //
	  
	    if ( i <= 15 ) {
	      if (loop == 0) hHSX1->SetBinContent(16-i,1,pe_h[i]);
	      if (pe_h[i] > hodoThre
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
	      if (loop == 0) hHSY1->SetBinContent(1,32-i,pe_h[i]);
	      if(max_hsy1<pe_h[i]){
		max_hsy1=pe_h[i];
		gapcheck_hsy1_ch=i;
	      }
	      if (pe_h[i] > hodoThre
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
	      if (loop == 0) hHSY2->SetBinContent(1,i-31,pe_h[i]);
	      if(max_hsy2<pe_h[i]){
		max_hsy2=pe_h[i];
		gapcheck_hsy2_ch =i;
	      }
	      if (pe_h[i] > hodoThre
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
	      if(max_hsx2<pe_h[i]){
		max_hsx2=pe_h[i];
		gapcheck_hsx2_ch=i;
	      }
	      if (loop == 0) hHSX2->SetBinContent(i-47,1,pe_h[i]);
	      if (pe_h[i] > hodoThre
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
	  //hHSCProf2D fill                                                                                    
	  if (gapcheck_hsx1_ch==(15-x) && gapcheck_hsy1_ch==(31-y) )
	    {
	      hHSCProf2D[x][y]->Fill(gapcheck_hsx2_ch-47,gapcheck_hsy2_ch-31);
	    }
	  

	  if (loop == 0) {
	    hNHit_HSX1->Fill(HitCount_HSX1);
	    hNHit_HSY1->Fill(HitCount_HSY1);
	    hNHit_HSY2->Fill(HitCount_HSY2);
	    hNHit_HSX2->Fill(HitCount_HSX2);
	  }
	  
	  if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2) good_event = true;
	  if (HitCount_HSX1 == 1 && HitCount_HSY1 == 1
	      && HitCount_HSY2 == 1 && HitCount_HSX2 == 1) SingleHit = true;

	  //PROTO 	  
	  double pThre = 9.5;//Thresholod of prototype PE
	  for ( int i = 0; i < 3; i++ ) {
	    
	    // pe_t[map[i]] = (double(adc_t[map[i]]));

	    pe_p1[map[i]] = double(floor(pe_p1[map[i]]*10)/10);
	    pe_p2[map[i]] = double(floor(pe_p2[map[i]]*10)/10);
	    
	    if ( i == 0 ) {
	      pe_p2[map[i]] = (double(adc_p2[map[i]]) - mean_ped[1][map[i]]) / gain[1][map[i]];
	      if (loop == 0) hCubeX->SetBinContent(1, 1, pe_p2[map[i]]);
	      if (good_event) {
		hPE_cell_all[x][y]->Fill(pe_p2[map[i]]);
		hPE_cell_x[x][y]->Fill(pe_p2[map[i]]);
		if (pe_p2[map[i]] > pThre) {
		  trig_X = true;
		  hHitRate_X->Fill(1, 1);
		  hPE_X->Fill(pe_p2[map[i]]);
		  cell_x[x][y]++;
		}
	      }
	    }
	    if ( i == 1 ) {
	      pe_p2[map[i]] = (double(adc_p2[map[i]]) - mean_ped[1][map[i]]) / gain[1][map[i]];
	      if (loop == 0) hCubeY->SetBinContent(1, 1, pe_p2[map[i]]);
	      if (good_event) {
		hPE_cell_all[x][y]->Fill(pe_p2[map[i]]);
		hPE_cell_y[x][y]->Fill(pe_p2[map[i]]);
		if (pe_p2[map[i]] > pThre) {
		  trig_Y = true;
		  hHitRate_Y->Fill(1, 1);
		  hPE_Y->Fill(pe_p2[map[i]]);
		  cell_y[x][y]++;
		}
	      }
	    }
	    if ( i == 2 ) {
	      pe_p1[map[i]] = (double(adc_p1[map[i]]) - mean_ped[0][map[i]]) / gain[0][map[i]];
	      if (loop == 0) hCubeZ->SetBinContent(1, 1, pe_p1[map[i]]);
	      if (good_event) {
		hPE_cell_all[x][y]->Fill(pe_p1[map[i]]);
		hPE_cell_z[x][y]->Fill(pe_p1[map[i]]);
		if (pe_p1[map[i]] > pThre) {
		  trig_Z = true;
		  hHitRate_Z->Fill(1, 1);
		  hPE_Z->Fill(pe_p1[map[i]]);
		  cell_z[x][y]++;
		}
	      }
	    }
	    
	  }


	  //add count gap_count & Fill 18/11/15 12:02
	  if(trig_X || trig_Y || trig_Z){
	    if((trig_X && trig_Y) || trig_Z){
	      
	      hHSC_gapcheck_2D[j-runnum]->Fill(15-Gcheck_hsx1_ch+1,31-Gcheck_hsy1_ch+1);
	      
	      if((15-Gcheck_hsx1_ch+1)>=6 && (15-Gcheck_hsx1_ch+1)<=11 && (31-Gcheck_hsy1_ch+1)>=6 && (31-Gcheck_hsy1_ch+1)<=11){
		//do nothing
	      }
	      else{
		gap_count_hodo1s++;
	      }
	    }
	    if(trig_X && trig_Y && trig_Z){
	      //do nothing
	    }
	    else{
	      gap_count_1s2s++;
	    }
	  }

	  //
	       
	  if (trig_X && trig_Y && trig_Z) {
	    Dpass_count++;
	    DpassEvent = true;
	    // cout << "D-passing event# : " << evt << endl;
	    if (DpassEvent_num[0] == 0) DpassEvent_num[0] = evt;
	    else {
	      DpassEvent_num[1] = DpassEvent_num[0];
	      DpassEvent_num[0] = evt;
	    }
	    if (DpassEvent_num[1]-GoodEvent_num < 0) {
	      // cout << "Diff between PT & HS (bigger): " << DpassEvent_num[0] - GoodEvent_num << endl;
	    }
	  }
	  
	  if (trig_HSX1 || trig_HSY1) trig_count++;
	  if (trig_HSX1 && trig_HSY1) andtrig_count++;
	  if (trig_HSX2 || trig_HSY2) nontrig_count++;
	  if (trig_HSX2 && trig_HSY2) andnontrig_count++;
	  
	  if (good_event && SingleHit) {
	    good_count++;
	    GoodEvent_num = evt;
	    // cout << "GOOD event# : " << evt << endl;
	    
	    // cout << "Diff between PT & HS (smaller): " << DpassEvent_num[0] - GoodEvent_num << endl;
	    /*
	      dist_X = abs(maxCh_HSX2-maxCh_HSX1);
	      dist_Y = abs(maxCh_HSY2-maxCh_HSY1);
	      dist = sqrt(pow(dist_X,2) + pow(dist_Y,2));
	      
	      hDist->Fill(dist);
	    */
	    
	    if (loop == 0 && evt >= min_evt && evt <= max_evt) {
	      
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
	      hCubeX->GetYaxis()->SetTitle();
	      hCubeX->GetXaxis()->SetTitle();
	      hCubeX->GetYaxis()->SetTitleSize(0.03);
	      hCubeX->GetXaxis()->SetTitleSize(0.03);
	      hCubeX->GetYaxis()->SetLabelSize(0.04);
	      hCubeX->GetXaxis()->SetLabelSize(0.04);
	      hCubeX->GetXaxis()->SetNdivisions(0);
	      hCubeX->GetYaxis()->SetNdivisions(0);
	      hCubeX->SetMarkerSize(2.0);
	      hCubeX->Draw("text colz");
	      
	      c0->cd(7);
	      hCubeY->GetYaxis()->SetTitle();
	      hCubeY->GetXaxis()->SetTitle();
	      hCubeY->GetYaxis()->SetTitleSize(0.03);
	      hCubeY->GetXaxis()->SetTitleSize(0.03);
	      hCubeY->GetYaxis()->SetLabelSize(0.04);
	      hCubeY->GetXaxis()->SetLabelSize(0.04);
	      hCubeY->GetXaxis()->SetNdivisions(0);
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
	      hCubeZ->GetXaxis()->SetNdivisions(0);
	      hCubeZ->GetYaxis()->SetNdivisions(0);
	      hCubeZ->SetMarkerSize(2.0);
	      hCubeZ->Draw("text colz");
	      
	      c0->SaveAs(figname1);
	    }
	  }
	  
	  //add gap_count hist 18/11/15 12:02
	  gap_count_h_1s2s[j-runnum]->SetBinContent(evt,gap_count_1s2s);
	  gap_count_h_hodo1s[j-runnum]->SetBinContent(evt,gap_count_hodo1s);
	  
	}
	
	cout << all_count_x1  << " events passed X1" <<endl;
	cout << all_count_y1  << " events passed Y1" <<endl;
	cout << all_count_y2  << " events passed Y2" <<endl;
	cout << all_count_x2  << " events passed X2" <<endl;
  
	cout << "DAQ trigger counts : " << trig_count << endl;
	cout << "non-trigger hodoscope OR counts : " << nontrig_count << endl;
	cout << "trigger hodoscope AND counts : " << andtrig_count << endl;
	cout << "non-trigger hodoscope AND counts : " << andnontrig_count << endl;

	less_evt_sum += less_evt;
	good_count_sum += good_count;
	Dpass_count_sum += Dpass_count;
	
	cout << good_count << " good events out of " << less_evt
	     << " successfully processed" << endl;
	cout << "(in total : " << good_count_sum << " good events out of " << less_evt_sum
	     << " )" << endl;
	cout << Dpass_count << " D-passing events out of " << less_evt
	     << " successfully processed" << endl;
	cout << "(in total : " << Dpass_count_sum << " D-passing events out of " << less_evt_sum
	     << " )" << endl;
	
      } // end of x direction loop
    } // end of y direction loop
  } // end of each file loop
  
  ofstream outputtext5("1cube/75ch/array_PE.tsv");
  TCanvas*c1 = new TCanvas("c1","",800,600);
  TCanvas*c2 = new TCanvas("c2","",800,600);
  TString pdfname1;
  gStyle->SetOptFit(1111);

  if (loop == 1) {
    ofstream outputtext1("1cube/cellhit_dist.tsv");
    ofstream outputtext2("1cube/cell_PE_all.tsv");
    ofstream outputtext3("1cube/cell_PE_x.tsv");
    ofstream outputtext4("1cube/cell_PE_y.tsv");
    ofstream outputtext5("1cube/cell_PE_z.tsv");

    TCanvas*c6 = new TCanvas("c6","");
    c6->Clear();
    c6->Divide(1,2);
    for (int x=0; x<16; x++) {
      for (int y=0; y<16; y++) {

	eff_all[x][y] = (cell_x[x][y]+cell_y[x][y]+cell_z[x][y])/hPE_cell_all[x][y]->GetEntries();
	eff_x[x][y] = cell_x[x][y]/hPE_cell_x[x][y]->GetEntries();
	eff_y[x][y] = cell_y[x][y]/hPE_cell_y[x][y]->GetEntries();
	eff_z[x][y] = cell_z[x][y]/hPE_cell_z[x][y]->GetEntries();
       
	TString cellfilename1 = TString::Format("1cube/cell_PE_all/cell_all_%d_%d_PE.pdf", x, y);
	TString cellfilename2 = TString::Format("1cube/cell_PE_XYZ/cell_XYZ_%d_%d_PE.pdf", x, y);
       
	if ( x <= 5 || x >= 11 || y <= 5 || y >= 11 ) {
	  f2->SetParameters(3000, 0, 1.);
	  hPE_cell_all[x][y]->Fit("gaus_orig", "", "", -3, 3);
	  hPE_cell_x[x][y]->Fit("gaus_orig", "", "", -3, 3);
	  hPE_cell_y[x][y]->Fit("gaus_orig", "", "", -3, 3);
	  hPE_cell_z[x][y]->Fit("gaus_orig", "", "", -3, 3);
	} else {
	  f1->SetParameters(1000, 30, 1.);
	  //select fit range
	  // allpeak[x][y]=  hPE_cell_all[x][y]->GetMaximumBin();
	  //xpeak[x][y]=  hPE_cell_x[x][y]->GetMaximumBin();
	  //ypeak[x][y]=  hPE_cell_y[x][y]->GetMaximumBin();
	  //zpeak[x][y]=  hPE_cell_z[x][y]->GetMaximumBin();
	  hPE_cell_all[x][y]->Fit("landau_orig", "", "", 10, 60);
	  hPE_cell_x[x][y]->Fit("landau_orig", "", "",10 , 60);
	  hPE_cell_y[x][y]->Fit("landau_orig", "", "", 10, 60);
	  hPE_cell_z[x][y]->Fit("landau_orig", "", "", 10, 60);
	 
	}
	cell_PE_all[x][y] = f1->GetParameter(1);
	cell_PE_x[x][y] = f1->GetParameter(1);
	cell_PE_y[x][y] = f1->GetParameter(1);
	cell_PE_z[x][y] = f1->GetParameter(1);
       
	c1->cd();
	hPE_cell_all[x][y]->GetYaxis()->SetTitle("# of Events");
	hPE_cell_all[x][y]->GetXaxis()->SetTitle("p.e.");
	hPE_cell_all[x][y]->GetYaxis()->SetTitleSize(0.04);
	hPE_cell_all[x][y]->GetXaxis()->SetTitleSize(0.04);
	hPE_cell_all[x][y]->GetYaxis()->SetLabelSize(0.05);
	hPE_cell_all[x][y]->GetXaxis()->SetLabelSize(0.05);
	hPE_cell_all[x][y]->Draw();
       
	c1->SaveAs(cellfilename1);
	c1->Clear();

	c2->cd();
	c2->Clear();
	c2->Divide(2,2);

	c2->cd(1);
	hPE_cell_x[x][y]->GetYaxis()->SetTitle("# of Events");
	hPE_cell_x[x][y]->GetXaxis()->SetTitle("p.e.");
	hPE_cell_x[x][y]->GetYaxis()->SetTitleSize(0.04);
	hPE_cell_x[x][y]->GetXaxis()->SetTitleSize(0.04);
	hPE_cell_x[x][y]->GetYaxis()->SetLabelSize(0.05);
	hPE_cell_x[x][y]->GetXaxis()->SetLabelSize(0.05);
	hPE_cell_x[x][y]->Draw();

	c2->cd(4);
	hPE_cell_y[x][y]->GetYaxis()->SetTitle("# of Events");
	hPE_cell_y[x][y]->GetXaxis()->SetTitle("p.e.");
	hPE_cell_y[x][y]->GetYaxis()->SetTitleSize(0.04);
	hPE_cell_y[x][y]->GetXaxis()->SetTitleSize(0.04);
	hPE_cell_y[x][y]->GetYaxis()->SetLabelSize(0.05);
	hPE_cell_y[x][y]->GetXaxis()->SetLabelSize(0.05);
	hPE_cell_y[x][y]->Draw();

	c2->cd(2);
	hPE_cell_z[x][y]->GetYaxis()->SetTitle("# of Events");
	hPE_cell_z[x][y]->GetXaxis()->SetTitle("p.e.");
	hPE_cell_z[x][y]->GetYaxis()->SetTitleSize(0.04);
	hPE_cell_z[x][y]->GetXaxis()->SetTitleSize(0.04);
	hPE_cell_z[x][y]->GetYaxis()->SetLabelSize(0.05);
	hPE_cell_z[x][y]->GetXaxis()->SetLabelSize(0.05);
	hPE_cell_z[x][y]->Draw();     

	c2->SaveAs(cellfilename2);

       
	outputtext1 << x << " " << y << " " << cell_x[x][y] << " "
		    << cell_y[x][y] << " " << cell_z[x][y] << endl;
	outputtext2 << x << " " << y << " " << cell_PE_all[x][y] << " " << eff_all[x][y] << endl;
	outputtext3 << x << " " << y << " " << cell_PE_x[x][y] << " " << eff_x[x][y] << endl;
	outputtext4 << x << " " << y << " " << cell_PE_y[x][y] << " " << eff_y[x][y] << endl;
	outputtext5 << x << " " << y << " " << cell_PE_z[x][y] << " " << eff_z[x][y] << endl;
	c6 -> cd();
	hHSCProf2D[x][y]->GetYaxis()->SetTitle("Fiber along Y axis");
	hHSCProf2D[x][y]->GetXaxis()->SetTitle("Fiber along X axis");
	hHSCProf2D[x][y]->GetYaxis()->SetTitleSize(0.03);
	hHSCProf2D[x][y]->GetXaxis()->SetTitleSize(0.03);
	hHSCProf2D[x][y]->GetYaxis()->SetLabelSize(0.04);
	hHSCProf2D[x][y]->GetXaxis()->SetLabelSize(0.04);
	hHSCProf2D[x][y]->SetMarkerSize(2.0);
	hHSCProf2D[x][y]->Draw("TEXT COLZ");
	name = TString::Format("gap_check/Prof2D_%d_%d_run%d+%d.pdf",x,y,runnum,file_count-1);
	c6->SaveAs(name);
   
      } // end of y loop
    } // end of x loop
   
  } else {

    TCanvas*c3 = new TCanvas("c3","",1800,600);
    c3->Clear();
    c3->Divide(4,2);

    c3->cd(2);
    hHitRate_X->GetXaxis()->SetTitle("Cube# along Z");
    hHitRate_X->GetYaxis()->SetTitle("Cube# along Y");
    hHitRate_X->GetYaxis()->SetTitleSize(0.03);
    hHitRate_X->GetXaxis()->SetTitleSize(0.03);
    hHitRate_X->GetYaxis()->SetLabelSize(0.04);
    hHitRate_X->GetXaxis()->SetLabelSize(0.04);
    hHitRate_X->GetXaxis()->SetNdivisions(0);
    hHitRate_X->GetYaxis()->SetNdivisions(0);
    hHitRate_X->SetMarkerSize(2.0);
    hHitRate_X->Draw("text colz");

    c3->cd(7);
    hHitRate_Y->GetXaxis()->SetTitle("Cube# along X");
    hHitRate_Y->GetYaxis()->SetTitle("Cube# along Z");
    hHitRate_Y->GetYaxis()->SetTitleSize(0.03);
    hHitRate_Y->GetXaxis()->SetTitleSize(0.03);
    hHitRate_Y->GetYaxis()->SetLabelSize(0.04);
    hHitRate_Y->GetXaxis()->SetLabelSize(0.04);
    hHitRate_Y->GetXaxis()->SetNdivisions(0);
    hHitRate_Y->GetYaxis()->SetNdivisions(0);
    hHitRate_Y->SetMarkerSize(2.0);
    hHitRate_Y->Draw("text colz");

    c3->cd(3);
    hHitRate_Z->GetYaxis()->SetTitle("Cube# along X");
    hHitRate_Z->GetXaxis()->SetTitle("Cube# along Y");
    hHitRate_Z->GetYaxis()->SetTitleSize(0.03);
    hHitRate_Z->GetXaxis()->SetTitleSize(0.03);
    hHitRate_Z->GetYaxis()->SetLabelSize(0.04);
    hHitRate_Z->GetXaxis()->SetLabelSize(0.04);
    hHitRate_Z->GetXaxis()->SetNdivisions(0);
    hHitRate_Z->GetYaxis()->SetNdivisions(0);
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

    TString hratename = TString::Format("1cube/hitRate_%d+%d_%d.pdf", runnum, file_count-1,subnum);
    c3->SaveAs(hratename);


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

    TString multiname = TString::Format("1cube/hodoscope_HitMultiplicity_%d+%d_%d.pdf",
					runnum, file_count-1,subnum);
    c4->SaveAs(multiname);

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

    TString PEname = TString::Format("1cube/proto_PE_%d+%d_%d.pdf",
				     runnum,file_count-1,subnum);
    c5->SaveAs(PEname);
  }  
  //add c12 for gapcheck//add c12 for gapcheck
  for(int j=0;j<file_count;j++){
   
    TCanvas* c12 = new TCanvas("c12","",1800,1500);
    c12->Clear();
    c12->Divide(2,2);
    c12->cd(1);
  
    gap_count_h_hodo1s[j]->GetXaxis()->SetTitle("Number of Events");
    gap_count_h_hodo1s[j]->GetYaxis()->SetTitle("gap_cons_count");
    gap_count_h_hodo1s[j]->GetYaxis()->SetTitleSize(0.04);
    gap_count_h_hodo1s[j]->GetXaxis()->SetTitleSize(0.04);
    gap_count_h_hodo1s[j]->GetYaxis()->SetLabelSize(0.05);
    gap_count_h_hodo1s[j]->GetXaxis()->SetLabelSize(0.05);
    gap_count_h_hodo1s[j]->GetYaxis()->SetRangeUser(0,20000);
    gap_count_h_hodo1s[j]->Draw("P");
   
    c12->cd(2);
    gap_count_h_1s2s[j]->GetXaxis()->SetTitle("Number of Events");
    gap_count_h_1s2s[j]->GetYaxis()->SetTitle("gap_cons_count");
    gap_count_h_1s2s[j]->GetYaxis()->SetTitleSize(0.04);
    gap_count_h_1s2s[j]->GetXaxis()->SetTitleSize(0.04);
    gap_count_h_1s2s[j]->GetYaxis()->SetLabelSize(0.05);
    gap_count_h_1s2s[j]->GetXaxis()->SetLabelSize(0.05);
    gap_count_h_1s2s[j]->GetYaxis()->SetRangeUser(0,20000);
    gap_count_h_1s2s[j]->Draw("P");
   
    c12->cd(3);
    hHSC_gapcheck_2D[j]->GetYaxis()->SetTitle();
    hHSC_gapcheck_2D[j]->GetXaxis()->SetTitle();
    hHSC_gapcheck_2D[j]->GetYaxis()->SetTitleSize(0.03);
    hHSC_gapcheck_2D[j]->GetXaxis()->SetTitleSize(0.03);
    hHSC_gapcheck_2D[j]->GetYaxis()->SetLabelSize(0.04);
    hHSC_gapcheck_2D[j]->GetXaxis()->SetLabelSize(0.04);
    hHSC_gapcheck_2D[j]->GetXaxis()->SetNdivisions(5);
    hHSC_gapcheck_2D[j]->GetYaxis()->SetNdivisions(5);
    hHSC_gapcheck_2D[j]->SetMarkerSize(2.0);
    hHSC_gapcheck_2D[j]->Draw("text colz");
    
    TString c12_name = TString::Format("./1cube/gap_count_%d_%d.pdf",runnum+j,subnum);
    c12->SaveAs(c12_name.Data());
    
    //
  }
}

int main(int argc, char** argv)
{
  if (argc != 7) {
    cerr << "run_1cube <1st run#> <subnum>  <# of files> <min_evt> <max_evt> <loop>"
	 << endl;
    return -1;
  }

  int runnum = atoi(argv[1]);
  int subnum = atoi(argv[2]);
  //  string dir = argv[3];
  int file_count = atoi(argv[3]);
  int min_evt = atoi(argv[4]);
  int max_evt = atoi(argv[5]);
  int loop = atoi(argv[6]);

  run_1cube(runnum, subnum, file_count, min_evt, max_evt, loop);

  return 0;
}

