//modified by R.Fujita in 2018/10/30
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TPaveStats.h>
using namespace std;

void run_calib(string &type, int &nrun, int &nsub, int &file_count,
	       int &low, int &high, string &imgtype){

  gStyle->SetTitleBorderSize(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  Int_t fontid = 132;
  gStyle->SetStatFont(fontid);
  gStyle->SetLabelFont(fontid, "XYZ");
  gStyle->SetLabelFont(fontid, "");
  gStyle->SetTitleFont(fontid, "XYZ");
  gStyle->SetTitleFont(fontid, "");
  gStyle->SetTextFont(fontid);
  gStyle->SetFuncWidth(2);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPaintTextFormat("3.3f");

  TF1 *f1 = new TF1("gaus_ped","[0]*TMath::Gaus(x,[1],[2])");
  TF1 *f2 = new TF1("gaus_1pe","[0]*TMath::Gaus(x,[1],[2])");
  
  int adc[64];

  Double_t ped_const[64];
  Double_t pe_const[64];
  Double_t mean_ped[64];
  Double_t mean_1pe[64];
  Double_t gain[64];
  Double_t sigma_ped[64];
  Double_t sigma_1pe[64];
  Double_t chi2_ndf_ped[64];
  Double_t chi2_ndf_1pe[64];
 
  string rootfile_dir = "../tree_root/";
  string calibfile_dir = "calibfile/";
  
  for (int n=0; n<file_count; n++) { 
    TString inputfile_name_t = TString::Format("%s_%04d_%04d_tree2.root",
					       type.c_str(), nrun+n, nsub);
  std::string inputfile_name = std::string(inputfile_name_t);
  string::size_type pos1 = inputfile_name.find(".root");
  if(pos1 == string::npos) {
    cerr << inputfile_name << " is not a root file" << endl;
    return;
  }

  string::size_type pos2 = inputfile_name.find("hodo");
  string::size_type pos3 = inputfile_name.find("proto2s");
  string::size_type pos4 = inputfile_name.find("proto1s");
  if (pos2 == string::npos && pos3 == string::npos && pos4 == string::npos) {
    cerr << inputfile_name << " does not match calib config." << endl;
    return;
  }
  

  string calibfile_name(inputfile_name);
  calibfile_name.replace(pos1,6,"_calib");

  string inputfile = rootfile_dir + inputfile_name;
  
  TString calibdir;
  calibdir = TString::Format("%s%s", calibfile_dir.c_str(), inputfile_name.c_str());
  std::string calibdir_s = std::string(calibdir);
  calibdir_s = calibdir_s.substr(0, calibdir_s.length()-11);
  struct stat statBuf;

  TChain *tree1 = new TChain("tree");
  TString trname;

  tree1->Add(inputfile.c_str());
  tree1->SetBranchAddress("ADC", &adc);

 TH1I* h1[64];
 TString hname;  
 int binmax1[64];
 int binmax2[64];
 
 TCanvas*c1 = new TCanvas("c1","",800,600); //800,300
 c1->Clear();
 TCanvas*c2 = new TCanvas("c2","");
 c2->Clear();
 TCanvas*c3 = new TCanvas("c3","");
 c3->Clear();
 TCanvas*c4 = new TCanvas("c4","");
 c4->Clear();
 TCanvas*c5 = new TCanvas("c5","");
 c5->Clear();
 TCanvas*c6 = new TCanvas("c6","");
 c6->Clear();
 TCanvas*c7 = new TCanvas("c7","");
 c7->Clear();
 TCanvas*c8 = new TCanvas("c8","");
 c8->Clear();
 
 // c1->Divide(8,4);
 
for (int ch=0; ch < 64; ch++ ) {
  //  c1->cd(ch+1);
  hname = TString::Format("h1_%d",ch);  
  h1[ch] = new TH1I(hname, "", 4096, 0, 4096);
  for ( int evt=0; evt < tree1->GetEntries(); evt++ ) {
    tree1->GetEntry(evt);
    h1[ch]->Fill(adc[ch]);
  }
  
  c1->cd();
  //pedestal fit
  if (pos2 >= 0 && pos3 == string::npos
      && pos4 == string::npos) h1[ch]->GetXaxis()->SetRange(700,850);
  if (pos3 >= 0 && pos2 == string::npos
      && pos4 == string::npos) h1[ch]->GetXaxis()->SetRange(700,900);
  if (pos4 >= 0 && pos3 == string::npos
      && pos2 == string::npos) h1[ch]->GetXaxis()->SetRange(700,840);
  
  binmax1[ch] = h1[ch]->GetMaximumBin();
  ped_const[ch] = h1[ch]->Integral(binmax1[ch], binmax1[ch]); 
  f1->SetParameters(ped_const[ch], binmax1[ch], 5.);
  cout << ped_const[ch] << " " << binmax1[ch] << endl;
  h1[ch]->Fit("gaus_ped", "+", "", binmax1[ch]-12, binmax1[ch]+12);
  mean_ped[ch] = f1->GetParameter(1);
  sigma_ped[ch] = abs(f1->GetParameter(2));
  chi2_ndf_ped[ch] = f1->GetChisquare() / f1->GetNDF();

  
  //1p.e. peak fit
  if (pos2 >= 0 && pos3 == string::npos && pos4 == string::npos && ch >= 32) {
    h1[ch]->GetXaxis()->SetRange(binmax1[ch]+low+5,binmax1[ch]+high+5);
  } else {
    h1[ch]->GetXaxis()->SetRange(binmax1[ch]+low,binmax1[ch]+high);
  }
  binmax2[ch] = h1[ch]->GetMaximumBin();
  pe_const[ch] = h1[ch]->Integral(binmax2[ch], binmax2[ch]);
  f2->SetParameters(pe_const[ch], binmax2[ch], 7.);
  h1[ch]->Fit("gaus_1pe", "+", "", binmax2[ch]-8, binmax2[ch]+13);
  mean_1pe[ch] = f2->GetParameter(1);
  sigma_1pe[ch] = f2->GetParameter(2);
  chi2_ndf_1pe[ch] = f2->GetChisquare() / f2->GetNDF();
  gain[ch] = mean_1pe[ch] - mean_ped[ch];

  if (gain[ch] < 25 || gain[ch] > 50) {
    h1[ch]->GetXaxis()->SetRange(840,900);
    if (pos2 >= 0 && pos3 == string::npos && pos4 == string::npos) {
      f2->SetParameters(pe_const[ch], mean_ped[ch]+40, 7.);
      h1[ch]->Fit("gaus_1pe", "+", "", mean_ped[ch]+33, mean_ped[ch]+55);
    } else if ((pos3 >= 0 && pos2 == string::npos && pos4 == string::npos)
	       || (pos4 >= 0 && pos3 == string::npos && pos2 == string::npos)) {
      f2->SetParameters(pe_const[ch], mean_ped[ch]+30, 7.);
      h1[ch]->Fit("gaus_1pe", "+", "", mean_ped[ch]+23, mean_ped[ch]+45);
    }
    mean_1pe[ch] = f2->GetParameter(1);
    sigma_1pe[ch] = f2->GetParameter(2);
    chi2_ndf_1pe[ch] = f2->GetChisquare() / f2->GetNDF();
    gain[ch] = mean_1pe[ch] - mean_ped[ch];
  }

  
  h1[ch]->GetXaxis()->SetRange(700,1200);
  //  h1[ch]->GetXaxis()->SetRange(600,4000);
  
  h1[ch]->GetYaxis()->SetTitle("Entries");
  h1[ch]->GetXaxis()->SetTitle("ADC counts");
  h1[ch]->GetXaxis()->SetNdivisions(511);
  h1[ch]->GetYaxis()->SetNdivisions(511);
 
 
TString name1;
TString name2;
 
 if (stat(calibdir_s.c_str(),&statBuf) != 0)  mkdir(calibdir_s.c_str(), 0777);

 name1 = TString::Format("calib_ch%d (64ch readout)", ch);
 name2 = TString::Format("%s/%s_ch%d.%s", calibdir_s.c_str(), calibfile_name.c_str(),
			 ch,imgtype.c_str());

 h1[ch]->SetTitle(name1);
 gStyle->SetOptFit(1111);
 h1[ch]->Draw();
 gPad->SetLogy();
 c1->SaveAs(name2);
 }
 
 c2->cd();
TH1F *h2_2 = new TH1F("h2_2_ch32-63","",100,0.,100.);
h2_2->SetLineColor(kBlue);
TH1F *h2_1 = new TH1F("h2_1_ch0-31","",100,0.,100.);
h2_1->SetLineColor(kRed);
for (int n=32; n<64; n++) { h2_2->Fill(gain[n]); }
h2_2->GetXaxis()->SetTitle("1p.e. gain (ADC counts)");
h2_2->GetYaxis()->SetTitle("# of channels");
h2_2->GetYaxis()->SetTitleOffset(0.7);
h2_2->SetTitle("Gain distribution");
h2_2->Draw();
for (int n=0; n<32; n++) { h2_1->Fill(gain[n]); }
h2_1->Draw("sames");
c2->Update();
TPaveStats *st2_1 = (TPaveStats*)h2_1->FindObject("stats");
st2_1->SetTextColor(kRed);
TPaveStats *st2_2 = (TPaveStats*)h2_2->FindObject("stats");
st2_2->SetTextColor(kBlue);
st2_2->SetY1NDC(0.3);
st2_2->SetY2NDC(0.6);
c2->Modified();
TString namec2;
 namec2 = TString::Format("%s/%s_gain_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
c2->SaveAs(namec2);

 c3->cd();
TH1F *h3_2 = new TH1F("h3_2_ch32-63","",200,700.,900.);
h3_2->SetLineColor(kBlue);
TH1F *h3_1 = new TH1F("h3_1_ch0-31","",200,700.,900.);
h3_1->SetLineColor(kRed);
for (int n=32; n<64; n++) { h3_2->Fill(mean_ped[n]); }
h3_2->GetXaxis()->SetTitle("Pedestal (ADC counts)");
h3_2->GetYaxis()->SetTitle("# of channels");
h3_2->GetYaxis()->SetTitleOffset(0.7);
h3_2->SetTitle("Pedestal distribution");
h3_2->Draw();
for (int n=0; n<32; n++) { h3_1->Fill(mean_ped[n]); }
h3_1->Draw("sames");
c3->Update();
TPaveStats *st3_1 = (TPaveStats*)h3_1->FindObject("stats");
st3_1->SetTextColor(kRed);
TPaveStats *st3_2 = (TPaveStats*)h3_2->FindObject("stats");
st3_2->SetTextColor(kBlue);
st3_2->SetY1NDC(0.3);
st3_2->SetY2NDC(0.6);
c3->Modified();
TString namec3;
 namec3 = TString::Format("%s/%s_pedestal_dist.%s",calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
c3->SaveAs(namec3);

 c4->cd();
 TH1F *h4_1 = new TH1F("h4_1_ch0-31","",200,800.,1000.);
 h4_1->SetLineColor(kRed);
 TH1F *h4_2 = new TH1F("h4_2_ch32-63","",200,800.,1000.);
 h4_2->SetLineColor(kBlue);
 for (int n=32; n<64; n++) { h4_2->Fill(mean_1pe[n]); }
 h4_2->GetXaxis()->SetTitle("1p.e. (ADC counts)");
 h4_2->GetYaxis()->SetTitle("# of channels");
 h4_2->GetYaxis()->SetTitleOffset(0.7);
 h4_2->SetTitle("1p.e. distribution");
 h4_2->Draw();
 for (int n=0; n<32; n++) { h4_1->Fill(mean_1pe[n]); }
 h4_1->Draw("sames");
 c4->Update();
 TPaveStats *st4_1 = (TPaveStats*)h4_1->FindObject("stats");
 st4_1->SetTextColor(kRed);
 TPaveStats *st4_2 = (TPaveStats*)h4_2->FindObject("stats");
 st4_2->SetTextColor(kBlue);
 st4_2->SetY1NDC(0.3);
 st4_2->SetY2NDC(0.6);
 c4->Modified();
 TString namec4;
 namec4 = TString::Format("%s/%s_1pe_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
 c4->SaveAs(namec4);
 
 c5->cd();
TH1F *h5_1 = new TH1F("h5_1_ch0-31","",100,0.,20.);
h5_1->SetLineColor(kRed);
TH1F *h5_2 = new TH1F("h5_2_ch32-63","",100,0.,20.);
h5_2->SetLineColor(kBlue);
for (int n=32; n<64; n++) { h5_2->Fill(sigma_ped[n]); }
h5_2->GetXaxis()->SetTitle("Pedestal_sigma (ADC counts)");
h5_2->GetYaxis()->SetTitle("# of channels");
h5_2->GetYaxis()->SetTitleOffset(0.7);
h5_2->SetTitle("Pedestal_sigma distribution");
h5_2->Draw();
for (int n=0; n<32; n++) { h5_1->Fill(sigma_ped[n]); }
h5_1->Draw("sames");
 c5->Update();
 TPaveStats *st5_1 = (TPaveStats*)h5_1->FindObject("stats");
 st5_1->SetTextColor(kRed);
 TPaveStats *st5_2 = (TPaveStats*)h5_2->FindObject("stats");
 st5_2->SetTextColor(kBlue);
 st5_2->SetY1NDC(0.3);
 st5_2->SetY2NDC(0.6);
 c5->Modified();
 TString namec5;
 namec5 = TString::Format("%s/%s_pedsigma_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
 c5->SaveAs(namec5);

 c6->cd();
 TH1F *h6_1 = new TH1F("h6_1_ch0-31","",100,0.,20.);
 h6_1->SetLineColor(kRed);
 TH1F *h6_2 = new TH1F("h6_2_ch32-63","",100,0.,20.);
 h6_2->SetLineColor(kBlue);
 for (int n=32; n<64; n++) { h6_2->Fill(sigma_1pe[n]); }
 h6_2->GetXaxis()->SetTitle("1p.e._sigma (ADC counts)");
 h6_2->GetYaxis()->SetTitle("# of channels");
 h6_2->GetYaxis()->SetTitleOffset(0.7);
 h6_2->SetTitle("1p.e._sigma distribution");
 h6_2->Draw();
 for (int n=0; n<32; n++) { h6_1->Fill(sigma_1pe[n]); }
 h6_1->Draw("sames");
 c6->Update();
 TPaveStats *st6_1 = (TPaveStats*)h6_1->FindObject("stats");
 st6_1->SetTextColor(kRed);
 TPaveStats *st6_2 = (TPaveStats*)h6_2->FindObject("stats");
 st6_2->SetTextColor(kBlue);
 st6_2->SetY1NDC(0.3);
 st6_2->SetY2NDC(0.6);
 c6->Modified();
 TString namec6;
 namec6 = TString::Format("%s/%s_1pesigma_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
 c6->SaveAs(namec6);

 c7->cd();
 TH1F *h7_1 = new TH1F("h7_1_ch0-31","",150,0.,150.);
 h7_1->SetLineColor(kRed);
 TH1F *h7_2 = new TH1F("h7_2_ch32-63","",150,0.,150.);
 h7_2->SetLineColor(kBlue);
 for (int n=0; n<32; n++) { h7_1->Fill(chi2_ndf_ped[n]); }
 h7_1->GetXaxis()->SetTitle("Pedestal_chi2/ndf (ADC counts)");
 h7_1->GetYaxis()->SetTitle("# of channels");
 h7_1->GetYaxis()->SetTitleOffset(0.7);
 h7_1->SetTitle("Pedestal_chi2/ndf distribution");
 h7_1->Draw();
 for (int n=32; n<64; n++) { h7_2->Fill(chi2_ndf_ped[n]); }
 h7_2->Draw("sames");
 c7->Update();
 TPaveStats *st7_1 = (TPaveStats*)h7_1->FindObject("stats");
 st7_1->SetTextColor(kRed);
 TPaveStats *st7_2 = (TPaveStats*)h7_2->FindObject("stats");
 st7_2->SetTextColor(kBlue);
 st7_2->SetY1NDC(0.3);
 st7_2->SetY2NDC(0.6);
 c7->Modified();
 TString namec7;
 namec7 = TString::Format("%s/%s_pedchi2_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
 c7->SaveAs(namec7);

 c8->cd();
TH1F *h8_1 = new TH1F("h8_1_ch0-31","",100,0.,20.);
h8_1->SetLineColor(kRed);
TH1F *h8_2 = new TH1F("h8_2_ch32-63","",100,0.,20.);
h8_2->SetLineColor(kBlue);
for (int n=32; n<64; n++) { h8_2->Fill(chi2_ndf_1pe[n]); }
h8_2->GetXaxis()->SetTitle("1p.e._chi2/ndf (ADC counts)");
h8_2->GetYaxis()->SetTitle("# of channels");
h8_2->GetYaxis()->SetTitleOffset(0.7);
h8_2->SetTitle("1p.e._chi2/ndf distribution");
h8_2->Draw();
for (int n=0; n<32; n++) { h8_1->Fill(chi2_ndf_1pe[n]); }
h8_1->Draw("sames");
c8->Update();
TPaveStats *st8_1 = (TPaveStats*)h8_1->FindObject("stats");
st8_1->SetTextColor(kRed);
TPaveStats *st8_2 = (TPaveStats*)h8_2->FindObject("stats");
st8_2->SetTextColor(kBlue);
st8_2->SetY1NDC(0.3);
st8_2->SetY2NDC(0.6);
c8->Modified();
TString namec8;
 namec8 = TString::Format("%s/%s_1pechi2_dist.%s", calibdir_s.c_str(),
			  calibfile_name.c_str(),imgtype.c_str());
c8->SaveAs(namec8);

TString out;
out = TString::Format("%s%s.tsv", calibfile_dir.c_str(),calibfile_name.c_str());
ofstream outputfile(out);
 for (int n=0; n<64; n++) {
   outputfile << n << " " << mean_ped[n] << " " << sigma_ped[n] << " " << chi2_ndf_ped[n]
	      << " " << mean_1pe[n] << " " << sigma_1pe[n] << " " << chi2_ndf_1pe[n]
	      << " " << gain[n] << endl;
   delete h1[n];
 }
 delete h2_1; delete h2_2; delete h3_1; delete h3_2; delete h4_1; delete h4_2;
 delete h5_1; delete h5_2; delete h6_1; delete h6_2; delete h7_1; delete h7_2;
 delete h8_1; delete h8_2;
  }
}


int main(int argc, char** argv)
{
  if(argc != 8) {
    cerr << "run_calib <detector type> <nrun> <nsub> <# of files> <pedmax+LOW> <pedmax+HIGH>  <png, eps, pdf>" << endl;
    return -1;
  }

  string type = argv[1];
  int nrun = atoi(argv[2]);
  int nsub = atoi(argv[3]);
  int file_count = atoi(argv[4]);
  int low = atoi(argv[5]);
  int high = atoi(argv[6]);
  string imgtype = argv[7];

  run_calib(type, nrun, nsub, file_count, low, high, imgtype);
  return 0;
}
