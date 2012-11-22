#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TLegend.h>


// some constants
static const double BSMASS(5366.3), DSMASS(1968.49), KMASS(493.677),
  BDMASS(5279.53), DMASS(1869.62), PIMASS(139.57018), DSSTMASS(2112.34),
  KSTMASS(891.66), LBMASS(5620.2), LCMASS(2286.46), PMASS(938.27203);
  /*, RHOMASS(775.49), DSTMASS(2010.25), DSTMASS2(2460.1);*/


int kfactorp(TTree *intree, TTree *outtree, std::string mode);

void plotvar(TTree *tree, TH1 &hist, std::string var, int vidx=-1);


int main(int argc, char* argv[])
{
  // program arguments
  std::string mode;
  if (argc == 2) {
    mode = argv[1];
  } else {
    std::cout << "Not enough arguments!" << std::endl;
    return -1;
  }

  // read ntuple from file
  std::string fname = "smalltree-" + mode + ".root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("ftree"));

  // read generated 4-vectors
  fname = "eventtree-" + mode + ".root";
  TFile outfile(fname.c_str(), "read");
  TTree *outtree = dynamic_cast<TTree*>(outfile.Get("TwoBodyDecayGen_decaytree"));

  // B momentum distribution from MC
  TH1D hMC("hMC", "", 100, 0.0, 300.0);
  TH1D hgen("hgen", "", 100, 0.0, 300.0);

  hMC.SetLineColor(kAzure);
  hgen.SetLineColor(kRed);

  TLatex *label = new TLatex();
  label->SetTextSize(0.04);
  label->SetNDC(true);

  TLegend *legend = new TLegend(0.6, 0.65, 0.85, 0.45);
  legend->AddEntry(&hMC, "Monte Carlo", "l");
  legend->AddEntry(&hgen, "Generated", "l");
  legend->SetFillColor(0);
  legend->SetLineColor(0);

  plotvar(intree, hMC, "BsMom");
  plotvar(outtree, hgen, "particle_lvs", 0);
  hgen.Draw("hist");
  hMC.Draw("hist same");

  legend->SetHeader("Bs momentum");
  legend->Draw();

  fname = mode + "_Bs_mom_both.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  // reset histograms
  hMC.Reset("icesm");
  hgen.Reset("icesm");

  plotvar(intree, hMC, "DsMom");
  plotvar(outtree, hgen, "particle_lvs", 1);
  hMC.Draw("hist");
  hgen.Draw("hist same");

  legend->SetHeader("Dau1 momentum");
  legend->Draw();

  fname = mode + "_dau1_mom_both.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  // reset histograms
  hMC.Reset("icesm");
  hgen.Reset("icesm");

  plotvar(intree, hMC, "hMom");
  plotvar(outtree, hgen, "particle_lvs", 2);
  hMC.Draw("hist");
  hgen.Draw("hist same");

  legend->SetHeader("Dau2 momentum");
  legend->Draw();

  fname = mode + "_dau2_mom_both.png";
  gPad->Update();
  gPad->Print(fname.c_str());

  kfactorp(intree, outtree, mode);

  return 0;
}


int kfactorp(TTree *intree, TTree *outtree, std::string mode)
{
  intree->ResetBranchAddresses();
  outtree->ResetBranchAddresses();

  TLorentzVector *tru_BsMom = NULL;
  TLorentzVector *tru_hMom = NULL;
  TLorentzVector *tru_DsMom = NULL;
  TLorentzVector Bs_rec(0.0, 0.0, 0.0, 0.0);

  intree->SetBranchAddress("tru_BsMom", &tru_BsMom);
  intree->SetBranchAddress("tru_hMom" , &tru_hMom);
  intree->SetBranchAddress("tru_DsMom", &tru_DsMom);

  long ientries(intree->GetEntries());

  TH1D kfactorMCm("kfactorMCm", "", 1000, 0.9, 1.1);
  TH1D kfactorMCp("kfactorMCp", "", 1000, 0.9, 1.1);
  TH1D kfactorMCpm("kfactorMCpm", "", 1000, 0.9, 1.1);

  // TLatex *label = new TLatex();
  // label->SetTextSize(0.04);
  // label->SetNDC(true);

  double kfactorp(0.0), kfactorm(0.0), kfactorpm(0.0);

  for (unsigned i = 0; i < ientries; ++i) {
    intree->GetEntry(i);

    if (tru_BsMom->P() > 3E5) continue;

    tru_hMom->SetVectM(tru_hMom->Vect(), KMASS);
    Bs_rec = *tru_hMom + *tru_DsMom;

    kfactorp = Bs_rec.P() / tru_BsMom->P();
    kfactorm = tru_BsMom->M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactorMCm.Fill(kfactorm);
    kfactorMCp.Fill(kfactorp);
    kfactorMCpm.Fill(kfactorpm);
  }

  std::vector<TLorentzVector> *particle_lvs = NULL;
  TLorentzVector BsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector DsMom(0.0, 0.0, 0.0, 0.0);
  TLorentzVector hMom(0.0, 0.0, 0.0, 0.0);

  outtree->SetBranchAddress("particle_lvs", &particle_lvs);

  long oentries(outtree->GetEntries());

  // kfactor.Reset("ICESM");
  TH1D kfactortrum("kfactortrum", "", 1000, 0.85, 1.05);
  TH1D kfactortrup("kfactortrup", "", 1000, 0.85, 1.05);
  TH1D kfactortrupm("kfactortrupm", "", 1000, 0.85, 1.05);
  TH1D kfactortruf("kfactortruf", "", 1000, 0.85, 1.05);

  for (unsigned i = 0; i < oentries; ++i) {
    outtree->GetEntry(i);
    if (particle_lvs->size() < 3) continue;
    (*particle_lvs)[2].SetVectM((*particle_lvs)[2].Vect(), KMASS*1E-3);
    Bs_rec = (*particle_lvs)[1] + (*particle_lvs)[2];

    kfactorp = Bs_rec.P() / (*particle_lvs)[0].P();
    kfactorm = (*particle_lvs)[0].M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactortrum.Fill(kfactorm);
    kfactortrup.Fill(kfactorp);
    kfactortrupm.Fill(kfactorpm);

    // 10s of keV
    volatile int iBs[3] = {int((*particle_lvs)[0].X()*1E5), int((*particle_lvs)[0].Y()*1E5),
			   int((*particle_lvs)[0].Z()*1E5)};
    volatile int iDs[3] = {int((*particle_lvs)[1].X()*1E5), int((*particle_lvs)[1].Y()*1E5),
			   int((*particle_lvs)[1].Z()*1E5)};
    volatile int ih[3] = {int((*particle_lvs)[2].X()*1E5), int((*particle_lvs)[2].Y()*1E5),
			  int((*particle_lvs)[2].Z()*1E5)};

    float mBs((*particle_lvs)[0].M()*1E3),
      mDs((*particle_lvs)[1].M()*1E3),
      mh((*particle_lvs)[2].M()*1E3);

    // MeV
    double dBs[4] = {iBs[0]*1E-2, iBs[1]*1E-2, iBs[2]*1E-2, 0.};
    double dDs[4] = {iDs[0]*1E-2, iDs[1]*1E-2, iDs[2]*1E-2, 0.};
    double dh[4] = {ih[0]*1E-2, ih[1]*1E-2, ih[2]*1E-2, 0.};
    dBs[3] = std::sqrt(double(mBs) * double(mBs) + dBs[0] *dBs[0] +
		       dBs[1] * dBs[1] + dBs[2] * dBs[2]);
    dDs[3] = std::sqrt(double(mDs) * double(mDs) + dDs[0] *dDs[0] +
		       dDs[1] * dDs[1] + dDs[2] * dDs[2]);
    dh[3] = std::sqrt(double(mh) * double(mh) + dh[0] *dh[0] +
		       dh[1] * dh[1] + dh[2] * dh[2]);

    float fBs[4] = {float(dBs[0]), float(dBs[1]), float(dBs[2]), float(dBs[3])};
    float fDs[4] = {float(dDs[0]), float(dDs[1]), float(dDs[2]), float(dDs[3])};
    float fh[4] = {float(dh[0]), float(dh[1]), float(dh[2]), float(dh[3])};

    // GeV
    double ddBs[4] = {fBs[0]*1E-3, fBs[1]*1E-3, fBs[2]*1E-3, fBs[3]*1E-3};
    double ddDs[4] = {fDs[0]*1E-3, fDs[1]*1E-3, fDs[2]*1E-3, fDs[3]*1E-3};
    double ddh[4] =  {fh[0]*1E-3, fh[1]*1E-3, fh[2]*1E-3, fh[3]*1E-3};

    BsMom.SetXYZT(ddBs[0], ddBs[1], ddBs[2], ddBs[3]);
    DsMom.SetXYZT(ddDs[0], ddDs[1], ddDs[2], ddDs[3]);
    hMom.SetXYZT(ddh[0], ddh[1], ddh[2], ddh[3]);

    Bs_rec = DsMom + hMom;

    kfactorp = Bs_rec.P() / BsMom.P();
    kfactorm = BsMom.M() / Bs_rec.M();
    kfactorpm = kfactorp * kfactorm;

    kfactortruf.Fill(kfactorm);
  }

  kfactorMCm.SetLineColor(kAzure);
  kfactorMCp.SetLineColor(kAzure);
  kfactorMCpm.SetLineColor(kAzure);

  kfactortrum.SetLineColor(kRed);
  kfactortrup.SetLineColor(kRed);
  kfactortrupm.SetLineColor(kRed);

  kfactortruf.SetLineColor(kGreen);

  TLegend legend(0.2, 0.6, 0.45, 0.4);
  legend.AddEntry(&kfactorMCpm, "Monte Carlo", "l");
  legend.AddEntry(&kfactortrupm, "Generated", "l");
  legend.SetFillColor(0);
  legend.SetLineColor(0);

  std::string fname;
  fname = mode + "_kfactorm_both.png";
  kfactortrum.Draw("hist");
  kfactorMCm.Draw("hist same");
  legend.SetHeader("k-factor (m)");
  legend.Draw();
  gPad->Update();
  gPad->Print(fname.c_str());

  fname = mode + "_kfactorp_both.png";
  kfactortrup.Draw("hist");
  kfactorMCp.Draw("hist same");
  legend.SetHeader("k-factor (p)");
  legend.Draw();
  gPad->Update();
  gPad->Print(fname.c_str());

  fname = mode + "_kfactorpm_both.png";
  kfactortrupm.Draw("hist");
  kfactorMCpm.Draw("hist same");
  legend.SetHeader("k-factor (m/p)");
  legend.Draw();
  gPad->Update();
  gPad->Print(fname.c_str());

  fname = "hdump-" + mode + ".root";
  TFile hdump(fname.c_str(), "update");
  hdump.WriteTObject(&kfactorMCm, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactorMCp, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactorMCpm, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortrum, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortrup, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortrupm, NULL, "WriteDelete");
  hdump.WriteTObject(&kfactortruf, NULL, "WriteDelete");
  hdump.Close();

  return 0;
}


void plotvar(TTree *tree, TH1 &hist, std::string var, int vidx)
{
  tree->ResetBranchAddresses();
  TLorentzVector *lv = NULL;
  std::vector<TLorentzVector> *particle_lvs = NULL;

  if (vidx < 0) {
    tree->SetBranchAddress(var.c_str(), &lv);
  } else {
    tree->SetBranchAddress(var.c_str(), &particle_lvs);
  }

  long nentries(tree->GetEntries());
  for (unsigned i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    if (vidx < 0) {
      hist.Fill(1E-3 * lv->P());
    } else {
      hist.Fill((*particle_lvs)[vidx].P());
    }
  }

  return;
}
