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

// STL headers
#include <string>
#include <vector>
#include <deque>

// ROOT headers
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
  typedef std::pair<unsigned, double> chBFpair; /**< Channel id and B.F. pair */

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
   * <b>NB:</b> Only one level of decay vertices have been tested at
   * the moment.
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
   * @param nparts Number of particles in the decay tree (length of the array)
   * @param brfr Branching fraction for the channel
   *
   * @return Status
   */
  bool add_decay_channel(double *masses, unsigned nparts,
			 double brfr);

  /**
   * Return requested daughter decay node
   *
   * @param chid Decay channel id
   * @param dauid Daughter number (0 or 1)
   *
   * @return Pointer to the daughter decay node
   */
  TwoBodyDecayGen* get_daughter(unsigned chid, unsigned dauid);

  /**
   * Return BF for decay channel
   *
   * @param chid Decay channel id
   *
   * @return Branching fraction for decay channel
   */
  double get_brfr(unsigned chid);

  /**
   * Find leaf branches or decay nodes.
   *
   * This method traverses the decay tree and extracts the branching
   * fraction and the channel id from each node into a double-ended
   * queue.  It stops the queue everytime a leaf node is encountered.
   * The queue is then saved into a vector.
   *
   * <i>Implementation:</i> Follow each daughter node, until a leaf
   * node is found.  Return the depth of the recursive call, so in the
   * current function scope 0, increment it by 1 as we go up the tree.
   * Once a leaf node is found, branch off the search to the other
   * daughters one level above (recursion level = 1).
   *
   *       (recursion depth)
   *
   *       L3  - - -                     mother
   *                                    /      \
   *                                   /        \
   *       L2  - - -                  d1        d2 (leaf node)
   *                                 /  \      /  \
   *                                /    \    /    \
   *       L1  - - -   (leaf node) d1    d2  d3    d4
   *                              /  \
   *                             /    \
   *       L0  - - -            d1    d2
   *
   * <b>NB:</b> end points (d{2..4} @ L1 or d{1,2} @ L0) are stored as
   * NULL pointers
   *
   * @param brfrVec Vector with deque for each leaf branch / decay node
   * @param brfrQ Pointer to deque for each leaf branch / decay node
   *
   * @return Depth where leaf node was found (-ve numbers are invalid)
   */
  int find_leaf_nodes(std::vector<std::deque<chBFpair> > &brfrVec,
		       std::deque<chBFpair> &brfrQ);

  /**
   * Generate one event at a time
   *
   * @param momp Mother 4-momentum
   * @param particle_lvs std::vector used to return generated 4-momenta
   * @param chQ Queue with channels to generate
   *
   * @return Event weight
   */
  double generate(TLorentzVector &momp,
		  std::vector<TLorentzVector> &particle_lvs,
		  std::deque<chBFpair> chQ);

  /**
   * Return if the particle is in LHCb detector acceptance
   *
   * @param part_lv Particle 4-vector
   *
   * @return Inside LHCb acceptance or not
   */
  bool lv_in_LHCb(TLorentzVector &part_lv);

  /**
   * Generate arbitrary number of events
   *
   * @param nevents Number of events to generate
   * @param hmomp Template histogram for 3-momentum of the mother particle
   * @param hmomn Template histogram for pseudorapidity(η) of the mother particle
   *
   * @return Generated event tree
   */
  TTree* get_event_tree(unsigned nevents, TH1 *hmomp, TH1 *hmomn=NULL);

  /**
   * Print decay tree
   *
   * @param indent Spaces to indent (to denote decay level)
   */
  void print(unsigned indent=0);

private:

  /**
   * Print queue for debugging
   *
   * @param prefix String prefix in output
   * @param queue Queue to print
   */
  void _printQ(std::string prefix, std::deque<chBFpair> queue);

  /**
   * Print vector of queues for debugging
   *
   * @param prefix String prefix in output
   * @param queue Vector of queues to print
   */
  void _printQ(std::string prefix, std::vector<std::deque<chBFpair> > queue);

  static unsigned long long _count; /**< Debug message counter */
  TGenPhaseSpace _generator;	/**< Generator for the current decay vertex */
  double _mommass;		/**< Mother particle mass for the current decay vertex */
  double _daumasses[NDAUS];	/**< Array of the two daughter masses */
  DauNodeVec _dauchannels;	/**< Decay channels with BF (stored as pointers) */
};

#endif	// TWOBODYDECAYGEN_HXX
