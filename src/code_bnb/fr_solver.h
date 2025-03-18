#ifndef FR_SOLVER_H
#define FR_SOLVER_H

#include <vector>
#include <set>
#include <map>
#include <queue>
#include <cmath>

#include "selfish_routing.h"

struct BBNode
{
  BBNode(std::map<int,std::vector<std::pair<int,int> > > _edgesToUse, std::set<int> _completedDemands, std::map<int,std::vector<std::pair<int,int> > > _edgesToNotUse, int _frIndexToBranchOn, std::vector<std::vector<double>> _initialSrFlows, double _lowerBound, double _upperBound, int _counter) : edgesToUse(_edgesToUse), completedDemands(_completedDemands), edgesToNotUse(_edgesToNotUse), frIndexToBranchOn(_frIndexToBranchOn), initialSrFlows(_initialSrFlows), lowerBound(_lowerBound), upperBound(_upperBound), counter(_counter) {}

  void print() const;

  // A note on node processing order.
  // The initial lower bound can only change if it becomes infeasible.
  // This will only happen when all out of node are marked False.
  bool operator < (const BBNode& rhs) const
  {
    //if (lowerBound > rhs.lowerBound + 1)
    //{
    //  return true;
    //}
    if (counter < rhs.counter)
    {
      return false;
    }
    return true;
  }

  std::map<int,std::vector<std::pair<int,int> > > edgesToUse;
  std::set<int> completedDemands;
  std::map<int,std::vector<std::pair<int,int> > > edgesToNotUse;
  int frIndexToBranchOn;

  std::vector<std::vector<double>> initialSrFlows;

  double lowerBound;
  double upperBound;

  int counter;
};

class FRSolver
{
  public:
    FRSolver(NetworkDesignProblem problem);

    bool solveBranchAndBound();

  private:
    double primalHeuristic(const BBNode node, std::pair<int,int>& branchEdge);
    double calculateLowerBound(const BBNode node);

    double upperBound;
    double lowerBound;

    NetworkDesignProblem problem;
    SelfishRoutingSolver selfishRoutingSolver;

    std::priority_queue<BBNode> nodeQueue;
    double minLeafNodeLowerBound;
    int nodeCounter;
    int nodesProcessed;

    double extraMilliseconds;
};

#endif
