// Bret Larget
// October 16-18, 24-25, 2012
// November 13-14, 17-19, 27-29 2012
// December 4, 2012
//
// Revisions by Karl Broman, 4-5 June 2014, for calling from R
//
// bp.C
//
// There is a bipartite graph with two types of nodes, PNodes (part nodes) and WNodes (whole nodes).
// Each node is active or inactive.
// Each edge connects a whole node and a part node.
// In a legal state, the following hold:
//   If a whole node is active, then each of its connected part nodes is active.
//   If all of the connected part nodes of a whole node are active, then the whole node is active.

// Here, we consider a larger state space where any set of whole nodes may be active or not.
// Each part node is active if any of its connected whole nodes are active; keep track of the number of active neighbors.
// A whole node is illegal if it is inactive, but all of its connected part nodes are active.
// The state is illegal if any whole nodes are illegal.
// We keep track of the number of illegal whole nodes.
//

// The MCMC chain samples from a function over the larger space that is reduced in the illegal region
// by a penalty of size (penalty-size)*(# of illegal nodes).  The consequences are that the chain spends more time
// in the legal region, the conditional distribution of the legal region matches the desired posterior, but the
// chain can potentially mix more freely among relative islands in the legal region by moving through illegal states.

// For each node, there is a vector of edges and a count of the number of active neighbors.
// The state contains the vectors of whole and part nodes and the vector of edges, plus a count of illegal whole nodes.
// The MCMC algorithm works as follows:
//   1. Pick a whole node at random; add it to the set of changed nodes.
//   2. Change the active state of this whole node.
//   3. Loop through the vector of connected part nodes:
//      4. If the node is not already part of the list of changed nodes, add it to the list.
//      5. Increment or decrement the number of active neighbors.
//      6. If the part node becomes active (or inactive) as the number of connected neighbors changes from 0 to 1 (or from 1 to 0),
//         then loop through connected whole nodes:
//         7. If the node is not already part of the list of changed nodes, add it to the list.
//         8. Modify the count of active neighbors of the whole node.
//   9. Calculate the new likelihood.
//  10. Loop through the changed whole nodes.
//      If the legal status of the whole node changed, then change the state count of illegal whole nodes.
//  11. Accept or reject.
//  12. If rejected, restore all of the changed nodes and the state to previous values.
#define VERSION 1.0

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <ctime>

#include "alias.h" // code to use the alias method to efficiently sample from discrete distributions (choose whole node nonuniformly)
#include <R.h>
#include <R_ext/Print.h>

using namespace std;

// public class, part of State
// contains input parameters to regulate the MCMC
class MCMCParameters
{
public:
  int numGenerations;
  int numBurnin;
  int subSampleRate;
  double penalty; // penalty for each illegal node, subtractd from the log-likelihood
  string initialState;
  MCMCParameters() {}
};

// Public class, part of State
// Exists as a simple convenience to pass into functions as a single entity
// Contains model parameters that do not change
class ModelNumbers
{
public:
  double alpha, logAlpha, logOneMinusAlpha;
  double gamma, logGamma, logOneMinusGamma;
  double p, logP, logOneMinusP;
  ModelNumbers() {}
};

// Public class, part of State
// Exists as a simple convenience to pass into functions as a single entity
// Contains model parameters that do change, so it is simple to copy to save and restore
class StateNumbers
{
public:
  int numIllegalNodes;
  int numActiveWholeNodes;
  int numActivePartNodes;
  int numInactiveWholeNodes;
  int numActiveOne;
  int numActiveZero;
  int numInactiveOne;
  int numInactiveZero;
  double logLikelihood;
  double logPrior;

  StateNumbers() : numIllegalNodes(0),
		   numActiveWholeNodes(0),
		   numActivePartNodes(0),
		   numInactiveWholeNodes(0),
		   numActiveOne(0),
		   numActiveZero(0),
		   numInactiveOne(0),
		   numInactiveZero(0),
		   logLikelihood(0.0),
		   logPrior(0.0)
    {}
  double calculateLogLikelihood(ModelNumbers m) {
    logLikelihood = numActiveOne*m.logGamma + numActiveZero*m.logOneMinusGamma + numInactiveOne*m.logAlpha + numInactiveZero*m.logOneMinusAlpha;
    return logLikelihood;
  }
  double calculateLogPrior(ModelNumbers m) {
    logPrior = numActiveWholeNodes*m.logP + numInactiveWholeNodes*m.logOneMinusP;
    return logPrior;
  }
};

class PNode;
class WNode;
class Node;

// Edge contains pointers to the nodes to which it connects
// In this program, edges never change after they are created
class Edge
{
private:
  PNode* pnode;
  WNode* wnode;
public:
  Edge(PNode* p,WNode* w) : pnode(p), wnode(w) {}
  PNode* getPNode() { return pnode; }
  WNode* getWNode() { return wnode; }
};

// Class Node is the base class for both whole nodes (WNode) and part nodes (PNode).
// It contains data and methods that are common to both types of nodes.
// protected items can be modified directly by derived classes (without get set methods).
// virtual functions set to 0 are defined separately in each derived class.

class Node
{
protected:
  string name;
  int number;
  vector<Edge*> edges;
  int numActiveNeighbors;
  int savedNumActiveNeighbors;
public:
  Node() {};
  ~Node() {
    edges.clear();
  }
  string getName() { return name; }
  int getNumber() { return number; }
  int getNumActiveNeighbors() { return numActiveNeighbors; }
  void addEdge(Edge* e) { edges.push_back(e); }
  int size() { return edges.size(); }
  virtual void restore() = 0;
  virtual bool isWNode() = 0;
};

// Class PNode inherits from class Node.
// The response is either 0 or 1 and is used in the likleihood calculation
// Keep track of the number of active neighbors: active if this is > 0
class PNode : public Node
{
private:
  int response;
public:
  PNode(string s,int n,int r) {
    name = s;
    number = n;
    response = r;
    numActiveNeighbors = 0;
  }
  bool isWNode() { return false; }
  bool active() { return ( numActiveNeighbors > 0 ); }
  int getResponse() { return response; }
  void save(set<Node*>& changedNodes) {
    if ( changedNodes.count(this) == 0 ) {
      changedNodes.insert(this);
      savedNumActiveNeighbors = numActiveNeighbors;
    }
  }
  void restore() { numActiveNeighbors = savedNumActiveNeighbors; }
  void initialAddActiveNeighbor(StateNumbers&);
  void mcmcFlipWNode(set<Node*>&,StateNumbers&,bool);
  void checkNumActiveNeighbors(ostream&);
};

// Class WNode inherits from Node
// The bool active is the object that is changed by the MCMC method.
// Illegal if not active, but number of active neighbors equals number of neighbors (all neighbors are active).
class WNode : public Node
{
private:
  bool active;
  bool savedActive;
  int responseSum; // sum of response for connected part nodes
public:
  WNode(string s,int n) {
    name = s;
    number = n;
    active = false;
    numActiveNeighbors = 0;
    responseSum = 0;
  }
  bool isWNode() { return true; }
  bool getActive() { return active; }
  string getName() { return name; }
  void save(set<Node*>& changedNodes) {
    if ( changedNodes.count(this) == 0 ) {
      changedNodes.insert(this);
      savedActive = active;
      savedNumActiveNeighbors = numActiveNeighbors;
    }
  }
  void restore() {
    active = savedActive;
    numActiveNeighbors = savedNumActiveNeighbors;
  }
  int countResponse();
  int getResponseSum() { return responseSum; }
  void initialSetActive(StateNumbers&);
  void initialAddActiveNeighbor() { numActiveNeighbors++; }
  void initiate(string,double,StateNumbers&);
  void addActiveNeighbor(set<Node*>& changedNodes) {
    save(changedNodes);
    numActiveNeighbors++;
  }
  void subtractActiveNeighbor(set<Node*>& changedNodes) {
    save(changedNodes);
    numActiveNeighbors--;
  }
  bool illegal() { return ( (active == false) && (numActiveNeighbors == edges.size()) ); }
  bool savedIllegal() { return ( (savedActive == false) && (savedNumActiveNeighbors == edges.size()) ); }
  void mcmcFlipWNode(set<Node*>&,StateNumbers&);
  bool checkNumActiveNeighbors(ostream&);
};

// Class State contains all of the objects that are updated by the MCMC and are needed to assess the likelihood
class State
{
private:
  MCMCParameters mcmcParameters;
  ModelNumbers modelNumbers;
  StateNumbers stateNumbers, savedStateNumbers;
  set<int> activeWholeNodeIndices; // used to keep track of the active whole nodes for faster updating of active counts
                                   // not in stateNumbers to avoid copying the set
  vector<int> activeCounts; // vector of length wnodes.size() with counts that each wnode is active when sampled
  vector<PNode*> pnodes;
  vector<WNode*> wnodes;
  vector<Edge*> edges;
  int sampleSize;
  string rootFile;
  Alias wDist;
public:
  State() { sampleSize = 0; }
  ~State() {
    for ( vector<PNode*>::iterator p=pnodes.begin(); p!=pnodes.end(); ++p )
      delete *p;
    pnodes.clear();
    for ( vector<WNode*>::iterator p=wnodes.begin(); p!=wnodes.end(); ++p )
      delete *p;
    wnodes.clear();
    for ( vector<Edge*>::iterator p=edges.begin(); p!=edges.end(); ++p )
      delete *p;
    edges.clear();
  }
  void setRootFile(string name) { rootFile = name; }
  string getRootFile() { return rootFile; }
  void setAlpha(double a) {
    modelNumbers.alpha=a;
    modelNumbers.logAlpha = log(a);
    modelNumbers.logOneMinusAlpha = log(1.0-a);
  }
  void setGamma(double b) {
    modelNumbers.gamma=b;
    modelNumbers.logGamma = log(b);
    modelNumbers.logOneMinusGamma = log(1.0-b);
  }
  void setP(double pp) {
    modelNumbers.p=pp;
    modelNumbers.logP = log(pp);
    modelNumbers.logOneMinusP = log(1.0-pp);
  }
  int getNumGenerations() { return mcmcParameters.numGenerations; }
  void setNumGenerations(int n) { mcmcParameters.numGenerations = n; }
  int getNumBurnin() { return mcmcParameters.numBurnin; }
  void setNumBurnin(int n) { mcmcParameters.numBurnin = n; }
  int getSubSampleRate() { return mcmcParameters.subSampleRate; }
  void setSubSampleRate(int n) { mcmcParameters.subSampleRate = n; }
  double getPenalty() { return mcmcParameters.penalty; }
  void setPenalty(double x) { mcmcParameters.penalty = x; }
  string getInitialState() { return mcmcParameters.initialState; }
  void setInitialState(string x) { mcmcParameters.initialState = x; }
  int getSampleSize() { return sampleSize; }
  void getCounts(vector<int>& counts,vector<string>& names,vector<int>& response, vector<int>& degree) {
    for ( int i=0; i< wnodes.size(); ++i ) {
      counts.push_back(activeCounts[i]);
      names.push_back(wnodes[i]->getName());
      response.push_back(wnodes[i]->getResponseSum());
      degree.push_back(wnodes[i]->size());
    }
  }
  void resizePNodes(int n) { pnodes.resize(n); }
  void resizeWNodes(int n) { wnodes.resize(n); }
  void addPNode(PNode* p,int i) { pnodes[i] = p; }
  void addWNode(WNode* w, int i) { wnodes[i] = w; }
  PNode* getPNode(int i) { return pnodes[i]; }
  WNode* getWNode(int i) { return wnodes[i]; }
  void addEdge(int w,int p) { 
    Edge* e = new Edge(pnodes[p],wnodes[w]);
    edges.push_back(e);
    wnodes[w]->addEdge(e);
    pnodes[p]->addEdge(e);
  }
  Edge* getEdge(int i) { return edges[i]; }
  double calculateInitialLogLikelihood();
  double calculateLogLikelihood() { return stateNumbers.calculateLogLikelihood(modelNumbers); }
  double getLogLikelihood() { return stateNumbers.logLikelihood; }
  double getSavedLogLikelihood() { return savedStateNumbers.logLikelihood; }
  double calculateLogPrior() { return stateNumbers.calculateLogPrior(modelNumbers); }
  double getLogPrior() { return stateNumbers.logPrior; }
  double getSavedLogPrior() { return savedStateNumbers.logPrior; }
  void save();
  void restore(set<Node*>&);
  void updateIllegal(set<Node*>&);
  void print(ostream&);
  void printAllIllegalNodes(ostream&);
  void initiateActivities();
  double mcmcFlipWNode(bool);
  void check(ostream&);
  void countResponse();
  void activateAllIllegalNodes();
  int getNumWholeNodes() { return wnodes.size(); }
  int getNumPartNodes() { return pnodes.size(); }
  int getNumEdges() { return edges.size(); }
  int getNumResponseOne() { return stateNumbers.numActiveOne + stateNumbers.numInactiveOne; }
  int getNumResponseZero() { return stateNumbers.numActiveZero + stateNumbers.numInactiveZero; }
  int getNumActiveOne() { return stateNumbers.numActiveOne; }
  int getNumActiveZero() { return stateNumbers.numActiveZero; }
  int getNumInactiveOne() { return stateNumbers.numInactiveOne; }
  int getNumInactiveZero() { return stateNumbers.numInactiveZero; }
  int getNumIllegalNodes() { return stateNumbers.numIllegalNodes; }
  int getNumActiveWholeNodes() { return stateNumbers.numActiveWholeNodes; }
  int getNumActivePartNodes() { return stateNumbers.numActivePartNodes; }
// use weight proportional to 6*(1 + ResponseSum)/(2 + Degree) + 1/sqrt(Degree)
// why this weight?  not sure, but it picks small degree node more often and high response nodes more often
  void initiateAlias() {
    int n = getNumWholeNodes();
    vector<double> prob(n,0);
    vector<int> index(n,0);
    double sum = 0;
    for ( int i=0; i<n; ++i ) {
      int m = wnodes[i]->size();
      index[i] = i;
      prob[i] = 6*(wnodes[i]->getResponseSum() + 1.0) / (2.0 + m) + 1.0 / sqrt( m );
      sum += prob[i];
    }
    for ( int i=0; i<n; ++i )
      prob[i] /= sum;
    
    wDist.initiate(prob,index);
  }
  void writeCounts(ofstream& f) {
    for ( int i=0; i<wnodes.size(); ++i )
      f << wnodes[i]->getName() << " " << setw(8) << activeCounts[i] << setw(12) << setprecision(6) << activeCounts[i] / (double)(mcmcParameters.numGenerations) << endl;
  }
};

// Save all of the state-level objects that can change
void State::save() { savedStateNumbers = stateNumbers; }

// Loop through the changed nodes and restore them.
// Then restore the state-level objects.
void State::restore(set<Node*>& changedNodes)
{
  for ( set<Node*>::iterator p=changedNodes.begin(); p!= changedNodes.end(); ++p )
    (*p)->restore();
  changedNodes.clear();
  stateNumbers = savedStateNumbers;
}

// The MCMC proposal to change the active status of nodes may have changed legal status of whole nodes.
// After all the active changes are processed, loop through the changed nodes and update the legal status.
// As a single proposal can touch nodes multiple times, it is far simpler to make all changes in all nodes first
// and then go back and update the legal status of the changed whole nodes.
// (When I tried to do this as the active states changed, I had a bug that led to the incorrect [negative] number of illegal states.)
void State::updateIllegal(set<Node*>& changedNodes)
{
  for ( set<Node*>::iterator p=changedNodes.begin(); p!= changedNodes.end(); ++p ) {
    if ( (*p)->isWNode() ) {
      if ( ((WNode*)(*p))->illegal() && !( ((WNode*)(*p))->savedIllegal() ) ) {
	stateNumbers.numIllegalNodes++;
      }
      else if ( !( ((WNode*)(*p))->illegal() ) && ((WNode*)(*p))->savedIllegal() ) {
	stateNumbers.numIllegalNodes--;
      }
    }
  }
}

// Only used for debugging.
void PNode::checkNumActiveNeighbors(ostream& f)
{
  int sum=0;
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    sum += ( (*e)->getWNode()->getActive() );
  if ( sum != numActiveNeighbors )
    f << "Node " << number << " mismatch: true = " << sum << ", stored = " << numActiveNeighbors << endl;
}

// Only used for debugging.
bool WNode::checkNumActiveNeighbors(ostream& f)
{
  int sum=0;
  for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
    sum += ( (*e)->getPNode()->active() );
  if ( sum != numActiveNeighbors )
    f << "Node " << number << " mismatch: true = " << sum << ", stored = " << numActiveNeighbors << endl;
  return illegal();
}

// Only used for debugging.
void State::check(ostream& f)
{
  f << "Check Whole Nodes:" << endl;
  int sumIllegal=0;
  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p )
    sumIllegal += ( (*p)->checkNumActiveNeighbors(f) == true ? 1 : 0 );
  if ( sumIllegal != stateNumbers.numIllegalNodes )
    f << "Illegal: true = " << sumIllegal << ", stored = " << stateNumbers.numIllegalNodes << endl;
  f << "Check Part Nodes:" << endl;
  for ( vector<PNode*>::iterator p=pnodes.begin(); p != pnodes.end(); ++p )
    (*p)->checkNumActiveNeighbors(f);
}

// Only used for debugging.
void State::printAllIllegalNodes(ostream& f)
{
  REprintf("Illegal nodes:");
  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p )
    if ( (*p)->illegal() )
      REprintf(" %d", (*p)->getNumber());
  REprintf("\n");
}

// Prints the changeable state-level information on a single line.
void State::print(ostream& f)
{
  f << setprecision(4) << setw(12) << stateNumbers.logLikelihood << " "
    << setprecision(4) << setw(12) << stateNumbers.logPrior << " "
    << setw(5) << stateNumbers.numActiveOne << " "
    << setw(5) << stateNumbers.numActiveZero << " "
    << setw(5) << stateNumbers.numInactiveOne << " "
    << setw(5) << stateNumbers.numInactiveZero << " "
    << setw(5) << stateNumbers.numIllegalNodes << " "
    << setw(5) << stateNumbers.numActiveWholeNodes << " "
    << setw(5) << stateNumbers.numActivePartNodes << endl;
}

// Set node to be active.
// Update the number of active neighbors for connected PNodes.
// Have them update the number of active neighbors for their connected neighbors if activity changes.
// This function is only called when the initial state is set, so no need to save states for MCMC.
void WNode::initialSetActive(StateNumbers& stateNumbers)
{
  active = true;
  stateNumbers.numActiveWholeNodes++;
  stateNumbers.numInactiveWholeNodes--;
  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++e )
    (*e)->getPNode()->initialAddActiveNeighbor(stateNumbers);
}

// Called when a connected whole node is made active at initiation.
void PNode::initialAddActiveNeighbor(StateNumbers& stateNumbers)
{
  ++numActiveNeighbors;
  if ( numActiveNeighbors == 1 ) {
    for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
      (*e)->getWNode()->initialAddActiveNeighbor();
    stateNumbers.numActivePartNodes++;
  }
}

// Loop through connected pnodes and count responses
int WNode::countResponse()
{
  responseSum = 0;
  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++e )
    responseSum += (*e)->getPNode()->getResponse();
  return responseSum;
}

void WNode::initiate(string initialState,double p,StateNumbers& stateNumbers)
{
//  if want to start with all whole nodes inactive, do nothing
  if ( initialState == "inactive" )
    return;

//  if want to start with a random start with probability p for each whole node to be active, do this
  if ( initialState == "random" ) {
    if ( unif_rand() < p )
      initialSetActive(stateNumbers);
    return;
  }

//  if want to start with all whole nodes active with positive response > 0.4, do this
  if ( initialState == "high" ) {
    if ( getResponseSum() / (double)(size()) > 0.4 )
      initialSetActive(stateNumbers);
    return;
  }
}

// Loop through wnodes, for each count the number of connected pnodes with response==1
void State::countResponse()
{
  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p )
    (*p)->countResponse();
}

// Activate all illegal nodes
void State::activateAllIllegalNodes()
{
  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p )
    if ( (*p)->illegal() )
      (*p)->initialSetActive(stateNumbers);
}

// Call the initiate routine for each whole node.
// Then loop through and activate all illegal nodes.
// After activation is set for all nodes, loop through and count the number of initial illegal nodes (should be zero!!!)
void State::initiateActivities()
{
  stateNumbers.numInactiveWholeNodes = getNumWholeNodes();

  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p )
    (*p)->initiate(mcmcParameters.initialState,modelNumbers.p,stateNumbers);

  activateAllIllegalNodes();

  activeCounts.resize(wnodes.size(), 0);

// loop through whole nodes and: (1) set initial counts for illegal nodes (ought to be 0 from previous function)
//                               (2) insert numbers of active nodes into set
  for ( vector<WNode*>::iterator p=wnodes.begin(); p != wnodes.end(); ++p ) {
    if ( (*p)->illegal() )
      stateNumbers.numIllegalNodes++;
    if ( (*p)->getActive() )
      activeWholeNodeIndices.insert( (*p)->getNumber() );
  }
}

// Loop through the part nodes.
// Each part node is categorized in a 2 by 2 manner:
//   active or inactive
//   response = 0 or 1
// Count the number of part nodes in each category.
// MCMC updates may change these counts.
// The counts and alpha/gamma are all that is needed to compute the likelihood.
double State::calculateInitialLogLikelihood()
{
  for ( vector<PNode*>::iterator p=pnodes.begin(); p!= pnodes.end(); ++p ) {
    if ( (*p)->active() ) {
      if ( (*p)->getResponse() == 1 ) // active, y=1
	stateNumbers.numActiveOne++;
      else // active, y=0
	stateNumbers.numActiveZero++;
    }
    else {
      if ((*p)->getResponse() == 1 ) // inactive, y=1
	stateNumbers.numInactiveOne++;
      else // inactive, y=0
	stateNumbers.numInactiveZero++;
    }
  }
  return calculateLogLikelihood();
}

// Function called on the single whole node whose active status is flopped with the flipWNode MCMC update.
// It then loops through all of its neighbors to continue the update.
void WNode::mcmcFlipWNode(set<Node*>& changedNodes,StateNumbers& stateNumbers)
{
  save(changedNodes);
  if ( active ) {
    active = false;
    stateNumbers.numActiveWholeNodes--;
    stateNumbers.numInactiveWholeNodes++;
  }
  else {
    active = true;
    stateNumbers.numActiveWholeNodes++;
    stateNumbers.numInactiveWholeNodes--;
  }
  for ( vector<Edge*>::iterator e=edges.begin(); e!= edges.end(); ++e )
    (*e)->getPNode()->mcmcFlipWNode(changedNodes,stateNumbers,active);
}

// Function called when the part node is connected to the whole node that is flipped with the flipWNode MCMC update.
// Change number of active neighbors.
// If this changes the active status (by moving from 0 to 1 or vice versa), then loop through connected whole nodes.
void PNode::mcmcFlipWNode(set<Node*>& changedNodes,StateNumbers& stateNumbers,bool activeCallingNode)
{
  save(changedNodes);
  if ( activeCallingNode ) {
    numActiveNeighbors++;
    if ( numActiveNeighbors == 1 ) {
      stateNumbers.numActivePartNodes++;
      if ( response == 1 ) {
	stateNumbers.numActiveOne++;
	stateNumbers.numInactiveOne--;
      }
      else {
	stateNumbers.numActiveZero++;
	stateNumbers.numInactiveZero--;
      }
      for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
	(*e)->getWNode()->addActiveNeighbor(changedNodes);
    }
  }
  else {
    numActiveNeighbors--;
    if ( numActiveNeighbors == 0 ) {
      stateNumbers.numActivePartNodes--;
      if ( response == 1 ) {
	stateNumbers.numActiveOne--;
	stateNumbers.numInactiveOne++;
      }
      else {
	stateNumbers.numActiveZero--;
	stateNumbers.numInactiveZero++;
      }
      for ( vector<Edge*>::iterator e=edges.begin(); e!=edges.end(); ++e )
	(*e)->getWNode()->subtractActiveNeighbor(changedNodes);
    }
  }
}

// Run under modified space where modified log-likelihood = regular log-likelihood - penalty*(number of illegal states)
double State::mcmcFlipWNode(bool postBurnin)
{
  save();
  set<Node*> changedNodes;
// use alias method to select wnode
  int i = wDist.pick();
  (wnodes[i])->mcmcFlipWNode(changedNodes,stateNumbers);
  updateIllegal(changedNodes);
  double newLogLikelihood = calculateLogLikelihood();
  double newLogPrior = calculateLogPrior();
  double acceptanceProbability = exp( (newLogLikelihood + newLogPrior - mcmcParameters.penalty*stateNumbers.numIllegalNodes)
				      - (savedStateNumbers.logLikelihood + savedStateNumbers.logPrior - mcmcParameters.penalty*savedStateNumbers.numIllegalNodes) );
  if ( acceptanceProbability > 1 )
    acceptanceProbability = 1.0;
  if ( unif_rand() > acceptanceProbability ) { // reject proposal
    restore(changedNodes);
  }
  else { // accept proposal
    if ( activeWholeNodeIndices.count(i) == 1 )
      activeWholeNodeIndices.erase(i);
    else
      activeWholeNodeIndices.insert(i);
  }
// activeCounts method of storing is much faster (1 or 2 clicks versus 20+)
// change storage from whole nodes and looping through all to set of active node indices and counts with the state
  if ( postBurnin && stateNumbers.numIllegalNodes == 0 ) {
    sampleSize++;
    for ( set<int>::iterator p=activeWholeNodeIndices.begin(); p!= activeWholeNodeIndices.end(); ++p )
      activeCounts[*p]++;
  }

  return acceptanceProbability;
}

// Long messy function that does many things.
// Might be nice to separate into component functions.
// 1. Parse the command line; exit if failure.
// 2. Initialize random number generator.
// 3. Initialize state variables.
// 4. Read in whole nodes and initialize them in state.
// 5. Read in part nodes and initialize them in state.
// 6. Read in the edge file and initialize corresponding nodes in state.
// 7. Select a random starting state by selecting random active whole nodes.
// 8. Update the number of active neighbors for all nodes and update the number of illegal nodes.
// 9. Calculate the initial log-likelihood.
void initialize(int n_whole, char **whole_names,
                int n_part, char **part_names, int *part_activity,
                int n_edge, char **edge_whole, char **edge_part,
                double alpha, double gamma, double p,
                int nburn, int ngen, int sub,
                double penalty, char *initial,
                State &state)

{
  int numGen=0, numBurnin=0, subSampleRate=1, i;
  string initialState = "inactive";

  initialState = (string)initial;
  if ( (initialState != "inactive") && (initialState != "random") && (initialState != "high") ) {
    error("Error: flag --initial must be followed by one of these arguments [ inactive | random | high ]");
  }

  numGen = ngen;
  numBurnin = nburn;
  subSampleRate = sub;

  if ( alpha<=0 || alpha >=gamma)
    error("alpha must be > 0 and < gamma");
  if ( gamma<=alpha || gamma >= 1 )
    error("gamma must be >= alpha and < 1");
  if ( p<= 0 || p >= 1 )
    error("p must be > 0 and < 1");
  if ( numGen < 0 || numBurnin < 0 || subSampleRate < 0)
    error("numGen, numBurnin, and subSampleRate must be >= 0");

// Set alpha and gamma and p
  REprintf("Setting alpha, gamma, and p....");
  state.setAlpha(alpha);
  state.setGamma(gamma);
  state.setP(p);
  REprintf("done.\n");

// Set MCMC parameters
  state.setNumGenerations(numGen);
  state.setNumBurnin(numBurnin);
  state.setSubSampleRate(subSampleRate);
  state.setPenalty(penalty);
  state.setInitialState(initialState);

// Convert the whole nodes
  REprintf("Converting whole nodes....");
  map<string,int> wnodeNames;
  for(i=0; i<n_whole; i++) {
    string name = (string)(whole_names[i]);
    if(wnodeNames.count(name) > 0)
      error("Error: whole node %d contains a name already used.\n", i+1);
    wnodeNames[name] = i;
  }
  REprintf("done.\n");

// Convert the part nodes
  REprintf("Reading in part nodes....");
  map<string,int> pnodeNames;
  map<string,int> responses;
  for(i=0; i<n_part; i++) {
    string name = (string)(part_names[i]);
    if(pnodeNames.count(name) > 0)
      error("Error: part node %d contains a name already used.\n", i+1);
    pnodeNames[name] = i;
    responses[name] = part_activity[i];
  }
  REprintf("done.\n");

// Initialize nodes in State
  // 11/29/2012 fixed bug where edges were based on alphabetical order and not the order read in
  //            wnodes and pnodes now stored in the order read in, not the map order
  //            number should now match index in state
  REprintf("Initializing nodes in state....");
  state.resizePNodes(pnodeNames.size());
  state.resizeWNodes(wnodeNames.size());

  for ( map<string,int>::iterator p=wnodeNames.begin(); p !=wnodeNames.end(); ++p )
    state.addWNode(new WNode(p->first,p->second), p->second);

  map<string,int>::iterator r=responses.begin();
  for ( map<string,int>::iterator p=pnodeNames.begin(); p !=pnodeNames.end(); ++p) {
    state.addPNode(new PNode(p->first,p->second,r->second), p->second);
    ++r;
  }
  REprintf("done.\n");

// Convert edges and initialize edges in State
  REprintf("Converting and initializing edges....");
  for(i=0; i<n_edge; i++) {
    string wname = (string)(edge_whole[i]);
    string pname = (string)(edge_part[i]);
    if((wnodeNames.count(wname) == 0) || (pnodeNames.count(pname) == 0)) {
      error("Error: edge %d (%s - %s) contains an invalid node name.\n",
            i+1, edge_whole[i], edge_part[i]);
    }
    state.addEdge(wnodeNames[wname], pnodeNames[pname]);
  }
  REprintf("done.\n");

// Initialize activities
  REprintf("Initiating activities at random....");
  state.countResponse();
  state.initiateActivities();
  REprintf("done.\n");

  state.calculateInitialLogLikelihood();
  state.calculateLogPrior();

// Initialize alias for mcmcFlipWNode
  state.initiateAlias();

  REprintf("State Summary:\n");
  REprintf("   Number of whole nodes = %d\n", state.getNumWholeNodes());
  REprintf("   Number of part nodes  = %d\n", state.getNumPartNodes());
  REprintf("    Number Response == 1 = %d\n", state.getNumResponseOne());
  REprintf("    Number Response == 0 = %d\n", state.getNumResponseZero());
  REprintf("   Number of edges       = %d\n", state.getNumEdges());
}

void bp(int n_whole, char **whole_names,
        int n_part, char **part_names, int *part_activity,
        int n_edge, char **edge_whole, char **edge_part,
        double alpha, double gamma, double p,
        int nburn, int ngen, int sub,
        double penalty, char *initial,
        char **resultName, int resultName_size,
        double *resultProb,
        int *resultCount, int *resultSample,
        int *resultDegree, int *resultResponse,
        int n_samples, double **SamplesDouble, int **SamplesInt)
{
  State state;
  int cursample=0;

  time_t initTime;
  time(&initTime);

  initialize(n_whole, whole_names,
             n_part, part_names, part_activity,
             n_edge, edge_whole, edge_part,
             alpha,gamma,p,nburn,ngen,sub,penalty,initial,
             state);

  const int window=100000;
  time_t start,current;
  time(&start);

// burnin
  REprintf("Begin burnin runs.\n");
  for ( int k=1; k<=state.getNumBurnin(); ++k ) {

    R_CheckUserInterrupt(); /* check for ^C */

    double acceptanceProbability = state.mcmcFlipWNode(false);
    if ( k % state.getSubSampleRate() == 0 ) {
      if(cursample >= n_samples) {
        warning("cursample exceeds available slots\n");
      } else {
        SamplesDouble[0][cursample] = acceptanceProbability;
        SamplesDouble[1][cursample] = state.getLogLikelihood();
        SamplesDouble[2][cursample] = state.getLogPrior();

        SamplesInt[0][cursample] = state.getNumActiveOne();
        SamplesInt[1][cursample] = state.getNumActiveZero();
        SamplesInt[2][cursample] = state.getNumInactiveOne();
        SamplesInt[3][cursample] = state.getNumInactiveZero();
        SamplesInt[4][cursample] = state.getNumIllegalNodes();
        SamplesInt[5][cursample] = state.getNumActiveWholeNodes();
        SamplesInt[6][cursample] = state.getNumActivePartNodes();
      }
      cursample++;
    }
    if ( k % window == 0 ) {
      time(&current);
      double secondsPerGeneration = difftime(current,start)/k;
      double remainBurnMinutes = (state.getNumBurnin()-k) * secondsPerGeneration / 60;
      double remainSampleMinutes = state.getNumGenerations() * secondsPerGeneration / 60;
      REprintf("%5d", (int)((100.0*k)/state.getNumBurnin()));
	    REprintf("%% complete, estimated burnin time remaining = %d", (int)remainBurnMinutes);
	    REprintf(" minutes, total = %d", (int)(remainBurnMinutes+remainSampleMinutes) + 1);
      REprintf(" minutes.\n");
    }
  }

// sample
  REprintf("Begin sample runs.\n");
  time(&start);
  for ( int k=1; k<=state.getNumGenerations(); ++k ) {

    R_CheckUserInterrupt(); /* check for ^C */

    double acceptanceProbability = state.mcmcFlipWNode(true);
    if ( k % state.getSubSampleRate() == 0 ) {
      if(cursample >= n_samples) {
        warning("cursample exceeds available slots\n");
      } else {
        SamplesDouble[0][cursample] = acceptanceProbability;
        SamplesDouble[1][cursample] = state.getLogLikelihood();
        SamplesDouble[2][cursample] = state.getLogPrior();

        SamplesInt[0][cursample] = state.getNumActiveOne();
        SamplesInt[1][cursample] = state.getNumActiveZero();
        SamplesInt[2][cursample] = state.getNumInactiveOne();
        SamplesInt[3][cursample] = state.getNumInactiveZero();
        SamplesInt[4][cursample] = state.getNumIllegalNodes();
        SamplesInt[5][cursample] = state.getNumActiveWholeNodes();
        SamplesInt[6][cursample] = state.getNumActivePartNodes();
      }
      cursample++;
    }
    if ( k % window == 0 ) {
      time(&current);
      double secondsPerGeneration = difftime(current,start)/k;
      int remainSeconds = (int)((state.getNumGenerations()-k) * secondsPerGeneration);
      int remainMinutes = remainSeconds / 60 + 1;
      REprintf("%5d", (int)((100.0*k)/state.getNumGenerations()));
	    REprintf("%% complete, estimated time remaining = ");
      if ( remainMinutes > 1 )
        REprintf("%d minutes.\n", remainMinutes);
      else
        REprintf("%d seconds.\n", remainSeconds);
    }
  }
// save output
  vector<int> counts;
  vector<string> names;
  vector<int> response;
  vector<int> degree;
  state.getCounts(counts,names,response,degree);

  for ( int i=0; i<counts.size(); ++i ) {
    strncpy(resultName[i], names[i].c_str(), resultName_size);

    resultProb[i] = counts[i]/(double)(state.getSampleSize());
    resultCount[i] = counts[i];
    resultSample[i] = state.getSampleSize();
    resultDegree[i] = degree[i];
    resultResponse[i] = response[i];
  }

  time_t endTime;
  time(&endTime);
}

extern "C" {
  void R_bp(int *n_whole, char **whole_names,
            int *n_part, char **part_names, int *part_activity,
            int *n_edge, char **edge_whole, char **edge_part,
            double *alpha, double *gamma, double *p,
            int *nburn, int *ngen, int *sub, 
            double *penalty, char **initial,
            char **resultName, int *resultName_size,
            double *resultProb,
            int *resultCount, int *resultSample,
            int *resultDegree, int *resultResponse,
            int *n_samples, double *samplesDouble, int *samplesInt)
  {
    double *SamplesDouble[3];
    int *SamplesInt[7];
    int i;

    GetRNGstate();

    SamplesDouble[0] = samplesDouble;
    for(i=1; i<3; i++)
      SamplesDouble[i] = SamplesDouble[i-1] + *n_samples;
    SamplesInt[0] = samplesInt;
    for(i=1; i<7; i++)
      SamplesInt[i] = SamplesInt[i-1] + *n_samples;

    bp(*n_whole, whole_names,
       *n_part, part_names, part_activity,
       *n_edge, edge_whole, edge_part,
       *alpha, *gamma, *p,
       *nburn, *ngen, *sub,
       *penalty, *initial,
       resultName, *resultName_size,
       resultProb, resultCount, 
       resultSample, resultDegree, resultResponse,
       *n_samples, SamplesDouble, SamplesInt);

    PutRNGstate();
  }
}
