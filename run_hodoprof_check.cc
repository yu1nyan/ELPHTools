#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <TMath.h>
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

void run_hodoprof_check(int &runnum, int subnum, int &file_count) {

  gStyle->SetOptStat(1110);
  TH2F *hNHit_Cell_U = new TH2F("hNHit_Cell_U","Up",16,0.5,16.5,16,0.5,16.5);
  TH2F *hNHit_Cell_D = new TH2F("hNHit_Cell_D","Down",16,0.5,16.5,16,0.5,16.5);
  TH2F *hNHit_Cell_G = new TH2F("hNHit_Cell_G","GoodEvent",16,0.5,16.5,16,0.5,16.5);



 
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
  
  Double_t pe_h[64];
  Int_t adc_2s[64];
  Int_t adc_1s[64];
  Int_t adc_h[64];

  string pt1sfile;
  string pt2sfile;
  string hsfile;

  int loop_num = 1;
  
  int ch;
  double mp, sp, cnp, m1, s1, cn1, g;

  for (int j=runnum; j<runnum+file_count; j++) {

    TString pt1sfile_name = TString::Format("%sproto1s_%04d_000%d_tree2.root",
					  rootfile_dir.c_str(), j, subnum);
    TString pt2sfile_name = TString::Format("%sproto2s_%04d_000%d_tree2.root",
					  rootfile_dir.c_str(), j, subnum);
    TString hsfile_name = TString::Format("%shodo_%04d_000%d_tree2.root",
					  rootfile_dir.c_str(), j, subnum);
    
    std::string pt1sfile_name_s = std::string(pt1sfile_name);
    std::string pt2sfile_name_s = std::string(pt2sfile_name);
    std::string hsfile_name_s = std::string(hsfile_name);
    
    string::size_type pos1 = pt2sfile_name_s.find(".root");
    string::size_type pos2 = pt1sfile_name_s.find(".root");
    string::size_type pos3 = hsfile_name_s.find(".root");

    string pt1sfile_calib(pt1sfile_name);
    pt1sfile_calib.replace(pos2,6,"_calib");
    pt1sfile_calib = pt1sfile_calib.substr(13);


    string pt2sfile_calib(pt2sfile_name);
    pt2sfile_calib.replace(pos1,6,"_calib");
    pt2sfile_calib = pt2sfile_calib.substr(13);

    
    string hsfile_calib(hsfile_name);
    hsfile_calib.replace(pos3,6,"_calib");
    hsfile_calib = hsfile_calib.substr(13);
    

    pt1sfile = pt1sfile_name_s;
    pt2sfile = pt2sfile_name_s;
    hsfile = hsfile_name_s;


  TString pt2s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt2sfile_calib.c_str());
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


  TString pt1s_calibname = TString::Format("%s%s.tsv", calibfile_dir.c_str(), pt1sfile_calib.c_str());
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
 
    tree1->Add(pt2sfile.c_str());
    tree1->SetBranchAddress("ADC", &adc_2s);
    tree2->Add(pt1sfile.c_str());
    tree2->SetBranchAddress("ADC", &adc_1s);
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
    

    loop_num = 16;

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
	
	cout << "Protype2s : number of events = " << tree1->GetEntries() << endl;
	cout << "Protype1s : number of events = " << tree2->GetEntries() << endl;
	cout << "Hodoscope : number of events = " << tree3->GetEntries() << endl;

	cout << "Select X1:" << "[" << omit_lowx1+1 << ", " <<16-omit_highx1 <<"]" << endl;
	cout << "Select Y1:" << "[" << omit_lowy1+1 << ", " <<16-omit_highy1 <<"]" << endl;
	cout << "Select Y2:" << "[" << omit_lowy2+1 << ", " <<16-omit_highy2 <<"]" << endl;
	cout << "Select X2:" << "[" << omit_lowx2+1 << ", " <<16-omit_highx2 <<"]" << endl;

	int less_evt = 0;


	int nent[3] = {};
	nent[0] = tree1->GetEntries();
	nent[1] = tree2->GetEntries();
	nent[2] = tree3->GetEntries();

	less_evt = TMath::MinElement(3, nent);

	cout << "less_evt = " << less_evt << endl;
	
	for(int evt = 0; evt < less_evt; evt++) {

	  tree1->GetEntry(evt);
	  tree2->GetEntry(evt);
	  tree3->GetEntry(evt);

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
	 
	  
	  for ( int i=0; i<64; i++ ) {
	    //  pe_h[i] = (double(adc_h[i]));
	    pe_h[i] = (double(adc_h[i]) - mean_ped[2][i]) / gain[2][i];
	    pe_h[i] = double(floor(pe_h[i]*10)/10);
	    
	    if ( i <= 15 ) {

	      if (pe_h[i] > 2.5
		  && i >=(omit_lowx1+0) && i<=(15-omit_highx1)) {

		HitCount_HSX1++;
		all_count_x1++;
		trig_HSX1 = true;
		if (pe_h[i] > max_HSX1) {
		  max_HSX1 = pe_h[i];
		  maxCh_HSX1 = 15-i;
		}
	      }
	    }
	    else if ( i>=16 && i <= 31 ) {

	      if (pe_h[i] > 2.5
		  && i >=(16+omit_lowy1) && i <= (31-omit_highy1)) {
		
		HitCount_HSY1++;
		all_count_y1++;
		trig_HSY1 = true;
		if (pe_h[i] > max_HSY1) {
		  max_HSY1 = pe_h[i];
		  maxCh_HSY1 = 31-i;
		}
	      }
	    }
	    else if ( i >=32 && i <= 47 ) {

	      if (pe_h[i] > 2.5
		  && i >=(32+omit_lowy2) && i <= (47-omit_highy2)) {
	
		HitCount_HSY2++;
		all_count_y2++;
		trig_HSY2 = true;
		if (pe_h[i] > max_HSY2) {
		  max_HSY2 = pe_h[i];
		  maxCh_HSY2 = i-32;
		}
	      }
	    }
	    else if ( i >= 48 && i <= 63 ) {

	      if (pe_h[i] > 2.5
		  && i >= (48+omit_lowx2) && i <= (63-omit_highx2)) {
		
		HitCount_HSX2++;
		all_count_x2++;
		trig_HSX2 = true;
		if (pe_h[i] > max_HSX2) {
		  max_HSX2 = pe_h[i];  
		  maxCh_HSX2 = i-48;
		}
	      }
	    }
	    
	  }
	  
	  if (trig_HSX1 && trig_HSY1 && trig_HSY2 && trig_HSX2) good_event = true;
	  if (HitCount_HSX1 == 1 && HitCount_HSY1 == 1 && HitCount_HSY2 == 1 && HitCount_HSX2 == 1) SingleHit = true;
	
	  if (trig_HSX1 && trig_HSY1) hNHit_Cell_U->Fill(x+1,y+1);	  
	  if (trig_HSY2 && trig_HSX2) hNHit_Cell_D->Fill(x+1,y+1);	  
	  if (good_event && SingleHit) hNHit_Cell_G->Fill(x+1,y+1);	  



	}
	  
      } // end of x direction loop
    } // end of y direction loop
  } // end of each file loop
  

  ///OUTPUT


 TCanvas*c6 = new TCanvas("c6","",800,600);
 c6->Clear();
 hNHit_Cell_G->Draw("COLZ");
 TString cellname = TString::Format("hodoprof/hodoscope_profile_good_event_run%d_%d.pdf",runnum,subnum);
 c6->SaveAs(cellname);

 TCanvas*c7 = new TCanvas("c7","",800,600);
 c7->Clear();
 hNHit_Cell_U->Draw("COLZ");
 cellname = TString::Format("hodoprof/hodoscope_profile_a_run%d_%d.pdf",runnum,subnum);
 c7->SaveAs(cellname);


 TCanvas*c8 = new TCanvas("c8","",800,600);
 c8->Clear();
 hNHit_Cell_D->Draw("COLZ");
 cellname = TString::Format("hodoprof/hodoscope_profile_b_run%d_%d.pdf",runnum,subnum);
 c8->SaveAs(cellname);

}

int main(int argc, char** argv)
{
  if (argc != 4) {
    cerr << "run_9cubes <1st run#> <subrun#> <# of files>"
	 << endl;
    return -1;
  }

  int runnum = atoi(argv[1]);
  int subnum = atoi(argv[2]);
  int file_count = atoi(argv[3]);

  run_hodoprof_check(runnum,subnum, file_count);

return 0;
}
