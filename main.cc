#include <iostream>
#include <cstdlib>

#include <TFile.h>
#include "TwoBodyDecayGen.hxx"


int main(int argc, char* argv[])
{
  if (argc > 2) {
    std::cout << "Too many arguments!" << std::endl;
    return -1;
  }

  unsigned long nevents(100);
  if (argc == 2) {
    nevents = atol(argv[1]);
  }

  double masses[NDAUS] = {1.968, 0.494};
  TwoBodyDecayGen generator(masses);
  TTree* eventtree = generator.get_event_tree(nevents);
  std::cout << eventtree->GetEntries() << " events generated" << std::endl;

  // eventtree->Scan("evt_wt");
  eventtree->Print("all");

  TFile file("eventtree.root", "recreate");
  file.cd();
  file.WriteTObject(eventtree);
  file.Close();

  return 0;
}
