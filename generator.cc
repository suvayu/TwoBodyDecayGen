#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TPad.h>

#include "TwoBodyDecayGen.hxx"


// some constants
static const double BSMASS(5366.3), DSMASS(1968.49), KMASS(493.677),
  BDMASS(5279.53), DMASS(1869.62), PIMASS(139.57018), DSSTMASS(2112.34),
  KSTMASS(891.66), LBMASS(5620.2), LCMASS(2286.46), PMASS(938.27203);
  /*, RHOMASS(775.49), DSTMASS(2010.25), DSTMASS2(2460.1);*/


void usage(char * prog)
{
  std::cout << "Usage: $ " << prog << " <nevents> <mode>"
    " # args are case sensitive" << std::endl;
}


int main(int argc, char* argv[])
{
  // program arguments
  if (argc > 3) {
    std::cout << "Too many arguments!" << std::endl;
    usage(argv[0]);
    return -1;
  }

  int nevents(100);
  std::string mode;
  if (argc == 3) {
    nevents = atol(argv[1]);
    mode = argv[2];
  } else {
    std::cout << "Not enough arguments!" << std::endl;
    usage(argv[0]);
    return -1;
  }

  // read ntuple from file
  std::string fname = "smalltree-" + mode + ".root";
  TFile infile(fname.c_str(), "read");
  TTree *intree = dynamic_cast<TTree*>(infile.Get("ftree"));

  // make dataset from histogram
  TH1D Bsmomp("Bsmomp", "", 100, 0.0, 300.0);
  intree->Draw("1E-3*tru_BsMom.P()>>Bsmomp");
  gPad->Print("Bs_mom_template.png");

  TH1D Bsmomn("Bsmomn", "", 100, 1.0, 6.0);
  intree->Draw("tru_BsMom.Eta()>>Bsmomn");
  gPad->Print("Bs_eta_template.png");

  // ROOT file dump
  fname = "eventtree-" + mode + ".root";
  TFile *file = new TFile(fname.c_str(), "recreate");
  file->cd();

  // generator config
  // double Bmasses[NDAUS] = {1.968, 0.494}; // DsK
  // if ("DsPi" == mode or "DsstPi" == mode) Bmasses[1] = PIMASS * 1E-3; // Pi

  // double Dsstmasses[NDAUS] = {0.0, 0.0}; // 0 for gamma
  // if ("DsstPi" == mode) {
  //   Bmasses[0] = DSSTMASS * 1E-3; // DsstPi
  //   Dsstmasses[0] = DSMASS * 1E-3; // Dsgamma
  // }

  // TwoBodyDecayGen stgen(DSSTMASS * 1E-3, Dsstmasses);
  // TwoBodyDecayGen generator(BSMASS * 1E-3, Bmasses, &stgen);
  std::vector<double> masses;
  masses.push_back(BSMASS * 1E-3);
  if ("DsK" == mode) {
    masses.push_back(DSMASS * 1E-3);
    masses.push_back(KMASS * 1E-3);
  } else if ("DsPi" == mode) {
    masses.push_back(DSMASS * 1E-3);
    masses.push_back(PIMASS * 1E-3);
  } else if ("DsstPi" == mode) {
    masses.push_back(DSSTMASS * 1E-3);
    masses.push_back(PIMASS * 1E-3);
    masses.push_back(DSMASS * 1E-3);
    masses.push_back(0.0);
  }

  std::vector<double> masses2;
  if ("DsstPi" == mode) {
    masses2.push_back(BSMASS * 1E-3);
    masses2.push_back(DSSTMASS * 1E-3);
    masses2.push_back(PIMASS * 1E-3);
    masses2.push_back(DSMASS * 1E-3);
    masses2.push_back(PIMASS * 1E-3);
  }

  TwoBodyDecayGen generator(&masses[0], masses.size());
  if ("DsstPi" == mode) {
    generator.add_decay_channel(&masses2[0], masses2.size(), 0.05);
  }
  generator.print();

  // generate, print summary and dump to ROOT file
  TTree* eventtree = generator.get_event_tree(nevents, &Bsmomp, &Bsmomn);
  eventtree->Print("all");
  file->WriteTObject(eventtree);
  file->Close();

  return 0;
}
