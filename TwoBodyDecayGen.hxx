/**
 * @file   TwoBodyDecayGen.hxx
 * @author Suvayu Ali <Suvayu.Ali@cernNOSPAM.ch>
 * @date   Wed Oct 10 22:39:13 2012
 *
 * @brief  This class implements an interface to generate 2-body decays
 *
 *
 */

#ifndef TWOBODYDECAYGEN_HXX
#define TWOBODYDECAYGEN_HXX

#include <string>
#include <vector>

#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>

#define NDAUS 2			/**< Number of daughters, fixed to 2 */


/**
 * This class defines a sequential 2-body decay tree and generates
 * events using the ROOT class TGenPhaseSpace.
 *
 * Note that this is a very simple phase space event generator and is
 * not aware of any resonances.  The decay tree is stored as an array
 * with particle masses and pointers to TwoBodyDecayGen instances for
 * subsequent decay vertices.  These pointers are NULL for leaf nodes
 * in the decay tree.  There are several constructors to instantiate a
 * TwoBodyDecayGen object; use of the constructor which takes an array
 * with particle masses and length of the decay tree (also the length
 * of the particle mass array) is recommended for simplicity of use.
 * Look at the constructor documentation for more details on the
 * format.
 *
 * @author Suvayu Ali <Suvayu.Ali@cernNOSPAM.ch>
 * @date 2012-11-05 Mon
 *
 */

class TwoBodyDecayGen {
public:

  typedef std::pair<std::vector<TwoBodyDecayGen*>, double> DauNode; /**< Node with B.F. */
  typedef std::vector<DauNode> DauNodeVec; /**< Vector of nodes (channels) */

  /**
   * Constructor 1
   *
   * @param mommass Mass of the mother in GeV/c²
   * @param dau1mass Mass of the first daughter in GeV/c²
   * @param dau2mass Mass of the second daughter in GeV/c²
   * @param dau1 Pointer to TwoBodyDecayGen object for first daughter
   * @param dau2 Pointer to TwoBodyDecayGen object for second daughter
   */
  TwoBodyDecayGen(double mommass, double dau1mass, double dau2mass,
		  TwoBodyDecayGen *dau1=NULL,
		  TwoBodyDecayGen *dau2=NULL);

  /**
   * Constructor 2
   *
   * @param mommass Mass of the mother in GeV/c²
   * @param daumasses Array of doubles with mass of the two daughters in GeV/c²
   * @param dau1 Pointer to TwoBodyDecayGen object for first daughter
   * @param dau2 Pointer to TwoBodyDecayGen object for second daughter
   */
  TwoBodyDecayGen(double mommass, double *daumasses,
		  TwoBodyDecayGen *dau1=NULL,
		  TwoBodyDecayGen *dau2=NULL);

  /**
   * Constructor that takes an array with particle masses for the
   * entire decay tree.
   *
   * No need to create the different decay vertices in the right order
   * any more.  Instead, pass an array with all the particle masses
   * and the number of particles (length of the array).  The array
   * should look something like these:
   *
   * 1. Bs → Ds*(Dsγ)ρ(ππ)
   *    double masses[7] = {Bs, Ds*, ρ, Ds, γ, π, π};
   * 2. Bs → Ds*(Dsγ)π
   *    double masses[5] = {Bs, Ds*, π, Ds, γ};
   * 3. Bs → Ds*(Dsπ)π
   *    double masses[5] = {Bs, Ds*, π, Ds, π};
   *
   * NB: Only one level of decay vertices have been tested at the
   * moment.
   *
   * @param masses Array of doubles with mass of all the particles in GeV/c²
   * @param nparts Number of particles in the decay tree (length of the array)
   */
  TwoBodyDecayGen(double *masses, unsigned nparts);

  ~TwoBodyDecayGen() {}

  /**
   * Add a new decay channel
   *
   * @param masses Array of doubles with mass of all the particles in GeV/c²
   * @param nparts Number of particles in the decay tree (lenght of the array)
   * @param brfr Branching fraction for the channel
   *
   * @return Status
   */
  bool add_decay_channel(double *masses, unsigned nparts,
			 double brfr=1.0);

  /**
   * Generate one event at a time
   *
   * @param momp Mother 4-momentum
   * @param particle_lvs std::vector used to return generated 4-momenta
   * @param ich Decay channel to generate (default: loop over all).
   *            Passing a -ve number loops over all channels
   *
   * @return Event weight
   */
  double generate(TLorentzVector &momp,
		  std::vector<TLorentzVector> &particle_lvs,
		  unsigned ich=-1);

  /**
   * Generate arbitrary number of events
   *
   * @param nevents Number of events to generate
   * @param hmomp Histogram template for 3-momentum of the mother mother particle
   *
   * @return Generated event tree
   */
  TTree* get_event_tree(unsigned nevents, TH1 *hmomp);

  /**
   * Print decay tree
   *
   * @param indent Spaces to indent (to denote decay level)
   */
  void print(unsigned indent=0);

private:

  TGenPhaseSpace _generator;	/**< Generator for the current decay vertex */
  double _mommass;		/**< Mother particle mass for the current decay vertex */
  double _daumasses[NDAUS];	/**< Array of the two daughter masses */
  DauNodeVec _dauchannels;	/**< Decay channels with B.F. (stored as pointers) */
};

#endif	// TWOBODYDECAYGEN_HXX
