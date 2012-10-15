#include <iostream>

#include <TFile.h>
#include "TwoBodyDecayGen.hxx"


int main()
{
  double masses[NDAUS] = {1.968, 0.494};
  TwoBodyDecayGen generator(masses);
  TTree* eventtree = generator.get_event_tree(100);
  std::cout << eventtree->GetEntries() << " events generated" << std::endl;

  // eventtree->Scan("evt_wt");
  eventtree->Print("all");

  TFile file("eventtree.root", "recreate");
  file.cd();
  file.WriteTObject(eventtree);
  file.Close();

  return 0;
}
