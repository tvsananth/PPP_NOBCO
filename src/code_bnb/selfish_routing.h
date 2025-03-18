#ifndef SELFISH_H
#define SELFISH_H

#include <iostream>
#include <fstream>
#include <stdio.h>

#include <vector>
#include <set>
#include <map>

#include <cmath>

enum DistanceType
{
  GRADIENT = 1,
  FREE_TRAVEL = 2
};

enum FunctionType
{
  USER_OPTIMAL = 1,
  SYSTEM_OPTIMAL = 2
};

class NetworkDesignProblem
{
  public:
    NetworkDesignProblem(std::string nodeFileName, std::string edgeFileName, std::string outputLogFileName);
 
    std::string ologFileName; 

    int numNodes;
    int numEdges;

    std::set<int> demandNodes;
    std::set<int> exitNodes;
    std::vector<int> nodeDemands;

    std::vector<int> frDemandNodes;
    std::set<int> frEntryNodes;

    double alpha;
    double beta;
    std::vector<std::vector<double> > edgeFreeTravelTimes;
    std::vector<std::vector<double> > edgeCapacities;
    std::vector<std::vector<int> > numLanes;

    bool useInitialSRValues;
    bool useInitialFRValues;
    std::vector<std::vector<double> > initialSrFlows;
    std::set<std::pair<int,int> > initialFrEdges;
 
    std::map<int,std::pair<int,int> > edgeIndices;
};

class SelfishRoutingSolver
{
  public:
    SelfishRoutingSolver(NetworkDesignProblem _problem, bool shouldPrint);

    double solve(std::set<std::pair<int,int> > frEdges, FunctionType functionType, std::vector<std::vector<double>> initialSrFlows);
    double calculateShortestPathFromNodeToSetOfNodes(int demandNode, const std::set<int>& finishNodes, std::vector<int>& shortestPathByNodes, const std::vector<std::pair<int,int> > edgesNotToUse, DistanceType distanceType, std::set<std::pair<int,int> > frEdges) const;
    std::vector<std::vector<double> > getFlows() { return flows; }

  private:
    void findInitialSolution();
    double calculateFunctionValue(const std::vector<std::vector<double> >& flows, FunctionType functionType) const;
    double calculateTotalEvacuationTime(const std::vector<std::vector<double> >& flows) const;
    bool checkFlowConservation() const;

    // methods for gradient descent
    void calculateGradient(FunctionType functionType);
    void resetDescentFlows();
    void updateDescentFlows(const std::vector<int>& shortestPathByNodes);
    double gradientDescentTolerance = 0.0001;

    std::vector<std::vector<double> > flows;
    std::vector<std::vector<double> > gradient;
    std::vector<std::vector<double> > descentFlows;

    NetworkDesignProblem problem;
    bool shouldPrint;
};

#endif
