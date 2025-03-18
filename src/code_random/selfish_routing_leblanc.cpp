#include "selfish_routing_leblanc.h"

#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <vector>
#include <algorithm>
#include <map>

#include <cmath>

NetworkDesignProblem::NetworkDesignProblem(std::string nodeFileName, std::string edgeFileName)
{
  useInitialFRValues = false;
  useInitialSRValues = false;

  // initialize nodes
  std::ifstream nodeIfs(nodeFileName);
  std::string line;
  int nodeIndex = -1;
  while (std::getline(nodeIfs, line))
  {
    // skip header
    if (nodeIndex == -1)
    {
      nodeIndex = 0;
      continue;
    }
  
    std::istringstream iss(line);
    std::vector<std::string> row;
    std::string word;
    while (std::getline(iss, word, ','))
    {
      row.push_back(word);
    }
    int demand = std::stoi(row[2]);
    bool exit = (std::stoi(row[3]) == 1);
    bool injured = (std::stoi(row[4]) == 1);

    nodeDemands.push_back(demand);
    if (demand > 0)
    {
      demandNodes.insert(nodeIndex);
    }
    if (exit)
    {
      exitNodes.insert(nodeIndex);
      frEntryNodes.insert(nodeIndex);
    }
    if (injured)
    {
      frDemandNodes.push_back(nodeIndex);
    }

    nodeIndex = nodeIndex + 1;
  }
  numNodes = nodeDemands.size();

  // initialize edges
  edgeCapacities.resize(numNodes);
  edgeFreeTravelTimes.resize(numNodes);
  numLanes.resize(numNodes);
  initialSrFlows.resize(numNodes);
  for (int i=0; i<numNodes; ++i)
  {
    edgeCapacities[i].resize(numNodes);
    edgeFreeTravelTimes[i].resize(numNodes);
    numLanes[i].resize(numNodes);
    initialSrFlows[i].resize(numNodes);
  }

  std::ifstream edgeIfs(edgeFileName);
  int edgeIndex = -1;
  while (std::getline(edgeIfs, line))
  {
    // skip header
    if (edgeIndex == -1)
    {
      edgeIndex = 0;
      continue;
    }

    std::istringstream iss(line);
    std::vector<std::string> row;
    std::string word;
    while (std::getline(iss, word, ','))
    {
      row.push_back(word);
    }
    int from = std::stoi(row[0]);
    int to = std::stoi(row[1]);

    // make forward edges
    edgeFreeTravelTimes[from][to] = std::stod(row[2]);
    edgeCapacities[from][to] = std::stod(row[3]);
    numLanes[from][to] = std::stoi(row[4]);
    edgeIndices[edgeIndex] = std::make_pair(from,to);
    numEdges = edgeIndex;

    // make reverse edges
    //edgeFreeTravelTimes[to][from] = std::stod(row[2]);
    //edgeCapacities[to][from] = std::stod(row[3]);
    //numLanes[to][from] = std::stoi(row[4]);
    //edgeIndex = edgeIndex + 1;

    // initial SR flows if there
    double initialSrFlow = std::stod(row[5]);
    if (initialSrFlow > 0.00001)
    {
      initialSrFlows[from][to] = initialSrFlow;
      //useInitialSRValues = true;
    }

    // initial FR choices if there
    int FR = std::stoi(row[6]);
    if (FR > 0)
    {
      int min = std::min(from,to);
      int max = std::max(from, to);
      initialFrEdges.insert(std::make_pair(min,max));
      useInitialFRValues = true;
    }
 
    edgeIndex = edgeIndex + 1;
  }

  // initialize alpha and beta
  alpha = 0.15;
  beta = 4;

  std::cout << "num nodes: " << numNodes << std::endl;
  std::cout << "demands:" << std::endl;
  std::cout << "exit nodes:" << std::endl;
  for (int exitNode : exitNodes)
  {
    std::cout << exitNode << "," << std::endl;
  }
  std::cout << "demand nodes:" << std::endl;
  for (int demandNode : demandNodes)
  {
    std::cout << demandNode << "[" << nodeDemands[demandNode] << "]," << std::endl;
  }
  std::cout << "edge capacities:" << std::endl;
  for (int row=0; row<numNodes; ++row)
  {
    for (int col=0; col<numNodes; ++col)
    {
      if (edgeCapacities[row][col] > 0)
      {
        std::cout << row << "," << col << ": " << edgeCapacities[row][col] << std::endl;
      }
    }
  }
  std::cout << "edge free flow times:" << std::endl;
  for (int row=0; row<numNodes; ++row)
  {
    for (int col=0; col<numNodes; ++col)
    {
      if (edgeFreeTravelTimes[row][col] > 0)
      {
        std::cout << row << "," << col << ": " << edgeFreeTravelTimes[row][col] << std::endl;
      }
    }
  }
  std::cout << "num lanes:" << std::endl;
  for (int row=0; row<numNodes; ++row)
  {
    for (int col=0; col<numNodes; ++col)
    {
      if (numLanes[row][col] > 0)
      {
        std::cout << row << "," << col << ": " << numLanes[row][col] << std::endl;
      }
    }
  }
}

SelfishRoutingSolver::SelfishRoutingSolver(NetworkDesignProblem _problem, bool _shouldPrint) : problem(_problem), shouldPrint(_shouldPrint)
{
  // initialize flows
  flows.resize(_problem.numNodes);
  for (int i=0; i<_problem.numNodes; ++i)
  {
    flows[i].resize(_problem.numNodes);
  }

  // initialize gradient
  gradient.resize(_problem.numNodes);
  for (int i=0; i<_problem.numNodes; ++i)
  {
    gradient[i].resize(_problem.numNodes);
  }

  // initialize descent flows
  descentFlows.resize(_problem.numNodes);
  for (int i=0; i<_problem.numNodes; ++i)
  {
    descentFlows[i].resize(_problem.numNodes);
  }
}

double SelfishRoutingSolver::calculateFunctionValue(const std::vector<std::vector<double> >& flowValues, FunctionType functionType) const
{
  double returnValue = 0.0;
  for (int row=0; row<problem.numNodes; ++row)
  {
    for (int col=0; col<problem.numNodes; ++col)
    {
      double flow = flowValues[row][col];
      if (flow > 0.00001)
      {
        double t0 = problem.edgeFreeTravelTimes[row][col];
        double c0 = problem.edgeCapacities[row][col];
        double arcValue = t0 * flow + (problem.alpha*t0*std::pow(flow,(problem.beta+1)) / (std::pow(c0,problem.beta)*(problem.beta + 1)));
        if (functionType == FunctionType::SYSTEM_OPTIMAL)
        {
          arcValue = t0 * flow + (problem.alpha*t0*std::pow(flow,(problem.beta+1)) / std::pow(c0,problem.beta));
        }
        returnValue = returnValue + arcValue;
      }
    }
  }

  return returnValue;
}

double SelfishRoutingSolver::calculateTotalEvacuationTime(const std::vector<std::vector<double> >& flowValues) const
{
  double returnValue = 0.0;
  for (int row=0; row<problem.numNodes; ++row)
  {
    for (int col=0; col<problem.numNodes; ++col)
    {
      double flow = flowValues[row][col];
      if (flow > 0.00001)
      {
        double t0 = problem.edgeFreeTravelTimes[row][col];
        double c0 = problem.edgeCapacities[row][col];
        double arcValue = t0 * flow + (problem.alpha*t0*std::pow(flow,problem.beta+1) / std::pow(c0,problem.beta));
        returnValue = returnValue + arcValue;
        //std::cout << "t0: " << t0 << " c0: " << c0 << " flow0: " << flow << '\n';
      }
    }
  }

  return returnValue;
}

void SelfishRoutingSolver::findInitialSolution()
{
  if (problem.useInitialSRValues)
  {
    for (int row=0; row<problem.numNodes; ++row)
    {
      for (int col=0; col<problem.numNodes; ++col)
      {
        flows[row][col] = problem.initialSrFlows[row][col];
      }
    }
  }
  else
  {
    for (int row=0; row<problem.numNodes; ++row)
    {
      for (int col=0; col<problem.numNodes; ++col)
      {
        gradient[row][col] = problem.edgeFreeTravelTimes[row][col];
      }
    }

    for (int demandNode : problem.demandNodes)
    {
      std::vector<std::pair<int,int> > emptyEdgesNotToUse;
      std::vector<int> shortestPathByNodes;
      std::set<std::pair<int,int> > emptyFrEdges;
      calculateShortestPathFromNodeToSetOfNodes(demandNode, problem.exitNodes, shortestPathByNodes, emptyEdgesNotToUse, DistanceType::FREE_TRAVEL, emptyFrEdges);
      for (int index=0; index<(shortestPathByNodes.size()-1); ++index)
      {
        int fromNodeIndex = shortestPathByNodes[index];
        int toNodeIndex = shortestPathByNodes[index+1];
        flows[fromNodeIndex][toNodeIndex] = flows[fromNodeIndex][toNodeIndex] + problem.nodeDemands[demandNode];
      }
    }
  }

  if (shouldPrint)
  {
    for (int row=0; row<problem.numNodes; ++row)
    {
      for (int col=0; col<problem.numNodes; ++col)
      {
        if (flows[row][col] > 0)
        {
          std::cout << row << "," << col << ": " << flows[row][col] << std::endl;
        }
      }
    }
  }
}

void SelfishRoutingSolver::calculateGradient(FunctionType functionType)
{
  for (int row=0; row<problem.numNodes; ++row)
  {
    for (int col=0; col<problem.numNodes; ++col)
    {
      double flow = flows[row][col];
      double t0 = problem.edgeFreeTravelTimes[row][col];
      double c0 = problem.edgeCapacities[row][col];
      double gradientValue = t0 + (problem.alpha*t0*std::pow(flow,problem.beta) / std::pow(c0,problem.beta));
      if (functionType == FunctionType::SYSTEM_OPTIMAL)
      {
        gradientValue = t0 + ((problem.beta+1)*problem.alpha*t0*std::pow(flow,problem.beta) / std::pow(c0,problem.beta));
      }
      gradient[row][col] = gradientValue;
    }
  }
}

void SelfishRoutingSolver::resetDescentFlows()
{
  for (int row=0; row<problem.numNodes; ++row)
  {
    for (int col=0; col<problem.numNodes; ++col)
    {
      descentFlows[row][col] = 0.0;
    }
  }
}

double SelfishRoutingSolver::calculateShortestPathFromNodeToSetOfNodes(int fromNode, const std::set<int>& finishNodes, std::vector<int>& shortestPathByNodes, std::vector<std::pair<int,int> > edgesNotToUse, DistanceType distanceType, std::set<std::pair<int,int> > frEdges) const
{
  // Dijkstra from fromNode until first exitNode is reached
  // 1. mark all nodes unvisited and create set of unvisited nodes
  std::vector<int> bestPredecessors(problem.numNodes, -1);
  std::vector<bool> marked(problem.numNodes, false);
  std::set<int> unvisitedNodes;
  for (int nodeIndex=0; nodeIndex<problem.numNodes; ++nodeIndex)
  {
    unvisitedNodes.insert(nodeIndex);
  }

  // 2. Assign a distance value to all nodes, 0 for source fromNode, infinity for others
  double infinity = std::numeric_limits<double>::max();
  std::vector<double> distances(problem.numNodes, infinity);
  distances[fromNode] = 0.0;
  int currentNode = fromNode;

  // special case when path is fully given
  if (finishNodes.find(fromNode) != finishNodes.end())
  {
    std::cout << "shortest path trivial" << std::endl;
    return 0;
  }

  // 3. Begin loop of processing nodes
  while (true)
  {
    for (int toNode : unvisitedNodes)
    {
      // check if edge allowed
      std::pair<int,int> edge = std::make_pair(currentNode,toNode);
      if ((problem.edgeCapacities[currentNode][toNode] > 0.0001) && (std::find(edgesNotToUse.begin(), edgesNotToUse.end(),edge) == edgesNotToUse.end()))
      {
        // update distances
        double possibleDistance = distances[currentNode] + gradient[currentNode][toNode];
        if (distanceType == DistanceType::FREE_TRAVEL)
        {
          // mark FR edges as 0 distance
          if (frEdges.find(edge) != frEdges.end())
          {
            possibleDistance = distances[currentNode];
          }
          else
          {
            possibleDistance = distances[currentNode] + problem.edgeFreeTravelTimes[currentNode][toNode];
          }
        }

        if (possibleDistance < distances[toNode])
        {
          distances[toNode] = possibleDistance;
          bestPredecessors[toNode] = currentNode;
        }
      }
    }
    marked[currentNode] = true;
    unvisitedNodes.erase(currentNode);
    //std::cout << "node removed: " << currentNode << " with d: " << distances[currentNode] << std::endl;

    // 4. Choose new current node
    int oldCurrentNode = currentNode;
    double newBestDistance = infinity;
    for (int nextNode : unvisitedNodes)
    {
      if (distances[nextNode] < newBestDistance)
      {
        newBestDistance = distances[nextNode];
        currentNode = nextNode;
      }
    }

    if (finishNodes.find(currentNode) != finishNodes.end())
    {
      // Found closest exit node, construct shortest path using stored bestPredecessors
      int currNodeOnPath = currentNode;
      while (currNodeOnPath != fromNode)
      {
        shortestPathByNodes.push_back(currNodeOnPath);
        int nextNodeOnPath = bestPredecessors[currNodeOnPath];
        currNodeOnPath = nextNodeOnPath;
      }
      shortestPathByNodes.push_back(currNodeOnPath);

      std::reverse(shortestPathByNodes.begin(), shortestPathByNodes.end());
      return distances[shortestPathByNodes.back()];
    }

    // handle unconnected case
    if (oldCurrentNode == currentNode)
    {
      std::cout << "no shortest path, disconnected" << std::endl;
      return infinity;
    }
  }
}

void SelfishRoutingSolver::updateDescentFlows(const std::vector<int>& shortestPathByNodes)
{
  int demandNode = shortestPathByNodes[0];
  for (int index=0; index<(shortestPathByNodes.size()-1); ++index)
  {
    int fromNodeIndex = shortestPathByNodes[index];
    int toNodeIndex = shortestPathByNodes[index+1];
    descentFlows[fromNodeIndex][toNodeIndex] = descentFlows[fromNodeIndex][toNodeIndex] + problem.nodeDemands[demandNode];
  }
}

void SelfishRoutingSolver::solve(std::set<std::pair<int,int> > frEdges, FunctionType functionType, double &pcur_time, double &pcur_time_evac, std::vector<double> &pcur_sol_sr_double)
{
  // print to send to ananth for debugging
  /*
  for (int index=0; index<problem.numEdges; ++index)
  {
    auto pair = problem.edgeIndices[index];
    if (frEdges.find(pair) != frEdges.end())
    {
      std::cout << "1" << std::endl;
    }
    else
    {
      std::cout << "0" << std::endl;
    }
  }
  */

  // Decrease capacities by 1 lane for FR edges. Redo at the end
  for (auto frEdge : frEdges)
  {
    double edgeCapacity = problem.edgeCapacities[frEdge.first][frEdge.second];
    int numberOfLanes = problem.numLanes[frEdge.first][frEdge.second];
    problem.edgeCapacities[frEdge.first][frEdge.second] = edgeCapacity * (numberOfLanes - 1) / (numberOfLanes * 1.0);
 
    double edgeCapacityReverse = problem.edgeCapacities[frEdge.second][frEdge.first];
    int numberOfLanesReverse = problem.numLanes[frEdge.second][frEdge.first];
    problem.edgeCapacities[frEdge.second][frEdge.first] = edgeCapacityReverse * (numberOfLanesReverse - 1) / (numberOfLanesReverse * 1.0);
  }

  //for (auto frEdge : frEdges)
  //{
  //  std::cout << frEdge.first << " " << frEdge.second << '\n';
  //  std::cout << problem.edgeCapacities[frEdge.first][frEdge.second] << '\n';
  //  std::cout << frEdge.second << " " << frEdge.first << '\n';
  //  std::cout << problem.edgeCapacities[frEdge.second][frEdge.first] << '\n';
  // }

  if (shouldPrint)
  {
    std::cout << "find initial solution" << std::endl;
  }
  auto beginFindInitialSolution = std::chrono::high_resolution_clock::now();
  findInitialSolution();
  double functionValue = calculateFunctionValue(flows, functionType);
  double totalEvacuationTime = calculateTotalEvacuationTime(flows);
  auto endFindInitialSolution = std::chrono::high_resolution_clock::now();
  auto initialSolutionDuration = std::chrono::duration_cast<std::chrono::seconds>(endFindInitialSolution - beginFindInitialSolution);

  if (shouldPrint)
  {
    std::cout << "start gradient descent" << std::endl;
  }
  auto beginGradientDescent = std::chrono::high_resolution_clock::now();
  int iteration = 0;
  bool shouldContinue = true;
  while (shouldContinue)
  {
    if (shouldPrint)
    {
      std::cout << "iter: " << iteration << " val: " << functionValue << std::endl;
    }
    checkFlowConservation();

    // Calculate gradient, which will be edge costs for shortest paths
    calculateGradient(functionType);
    /*
    std::cout << "gradient: " << std::endl;
    for (int row=0; row<problem.numNodes; ++row)
    {
      for (int col=0; col<problem.numNodes; ++col)
      {
        if (gradient[row][col] > 0.001)
        {
          std::cout << row << "," << col << ": " << gradient[row][col] << std::endl;
        }
      }
    }
    */

    // For each demand node, find shortest path to exit, to get optimal directions
    resetDescentFlows();
    for (int demandNode : problem.demandNodes)
    {
      // use Dijkstra to calculate shortest paths
      std::vector<std::pair<int,int> > emptyEdgesNotToUse;
      std::vector<int> shortestPathByNodes;
      std::set<std::pair<int,int> > emptyFrEdges;
      calculateShortestPathFromNodeToSetOfNodes(demandNode, problem.exitNodes, shortestPathByNodes, emptyEdgesNotToUse, DistanceType::GRADIENT, emptyFrEdges);
      updateDescentFlows(shortestPathByNodes);
      /*
      std::cout << "descent flows: " << std::endl;
      for (int row=0; row<problem.numNodes; ++row)
      {
        for (int col=0; col<problem.numNodes; ++col)
        {
          if (descentFlows[row][col] > 0.001)
          {
            std::cout << row << "," << col << ": " << descentFlows[row][col] << std::endl;
          }
        }
      }
      */
    }

    // Bolzano search (bisecting method) to find best step size problem.alpha in [0,1]
    double lowerBound = 0;
    std::vector<std::vector<double> > lowerBoundFlows(flows);
    double lowerBoundValue = functionValue;

    double upperBound = 1;
    std::vector<std::vector<double> > upperBoundFlows(descentFlows);
    double upperBoundValue = calculateFunctionValue(upperBoundFlows, functionType);

    bool continueSearch = true;
    double bestFunctionValue = functionValue;
    std::vector<std::vector<double> > bestFlows(flows);
    while (continueSearch)
    {
      // step size is midpoint
      double midpointBound = (upperBound + lowerBound) / 2.0;

      // get new flow values and new function value
      std::vector<std::vector<double> > midpointFlows(flows);
      for (int row=0; row<problem.numNodes; ++row)
      {
        for (int col=0; col<problem.numNodes; ++col)
        {
          // direction is d1=y1-x1, so min f(x1 + a*d1)
          midpointFlows[row][col] = midpointFlows[row][col] + ((descentFlows[row][col] - midpointFlows[row][col]) * midpointBound);
        }
      }
      double midpointFunctionValue = calculateFunctionValue(midpointFlows, functionType);
      //std::cout << "step: " << midpointBound << " val: " << midpointFunctionValue << std::endl;

      // replace with worst of lower and upper bound
      if (upperBoundValue < lowerBoundValue)
      {
        lowerBound = midpointBound;
        lowerBoundValue = midpointFunctionValue;
        lowerBoundFlows = midpointFlows;
      }
      else
      {
        upperBound = midpointBound;
        upperBoundValue = midpointFunctionValue;
        upperBoundFlows = midpointFlows;
      }

      // end once the window is small enough
      if (upperBound - lowerBound < 0.00001)
      {
        bestFunctionValue = midpointFunctionValue;
        bestFlows = midpointFlows;
        continueSearch = false;
      }
    }

    // check termination condition, else set new bests and iterate
    if (100.0 * (functionValue - bestFunctionValue) / functionValue < gradientDescentTolerance)
    {
      shouldContinue = false;
    }
    else
    {
      functionValue = bestFunctionValue;
      flows = bestFlows;
      iteration = iteration + 1;
    }
  }

  if (shouldPrint)
  {
    std::cout << "optimal flows:" << std::endl;
    for (int row=0; row<problem.numNodes; ++row)
    {
      for (int col=0; col<problem.numNodes; ++col)
      {
        if (flows[row][col] > 0.00001)
        {
          std::cout << row << "," << col << ": " << flows[row][col] << std::endl;
        }
      }
    }
  }

  checkFlowConservation();

  // Look at total evacuation times or equilibrium times?
  // We want to avoid choosing that road for FR next time.
  double optimalFlowsValue = calculateFunctionValue(flows, functionType);
  double optimalTotalEvacuationTime = calculateTotalEvacuationTime(flows);
  if (shouldPrint)
  {
    std::cout << "optimal flows value: " << optimalFlowsValue << std::endl;
    std::cout << "total evacuation time: " << optimalTotalEvacuationTime << std::endl;
  }

  auto endGradientDescent = std::chrono::high_resolution_clock::now();
  auto gradientDescentDuration = std::chrono::duration_cast<std::chrono::seconds>(endGradientDescent - beginGradientDescent);

  if (shouldPrint)
  {
    std::cout << "initial solution took: " << initialSolutionDuration.count() << " seconds" << std::endl;
    std::cout << "gradient descent took: " << gradientDescentDuration.count() << " seconds" << std::endl;
    std::cout << "total solve time took: " << initialSolutionDuration.count() + gradientDescentDuration.count() << " seconds" << std::endl;
  }

  // Restore capacities of FR lanes
  for (auto frEdge : frEdges)
  {
    double edgeCapacity = problem.edgeCapacities[frEdge.first][frEdge.second];
    int numberOfLanes = problem.numLanes[frEdge.first][frEdge.second];
    problem.edgeCapacities[frEdge.first][frEdge.second] = edgeCapacity * numberOfLanes / ((numberOfLanes - 1) * 1.0);
 
    double edgeCapacityReverse = problem.edgeCapacities[frEdge.second][frEdge.first];
    int numberOfLanesReverse = problem.numLanes[frEdge.second][frEdge.first];
    problem.edgeCapacities[frEdge.second][frEdge.first] = edgeCapacityReverse * numberOfLanesReverse / ((numberOfLanesReverse - 1) * 1.0);
  }

  if (shouldPrint)
  {
    std::cout << " pcur: " << pcur_time << " pcur_time: " << pcur_time_evac << '\n';
  }
  pcur_time = optimalFlowsValue;
  pcur_time_evac = optimalTotalEvacuationTime;
  if (shouldPrint)
  {  
    std::cout << " pcur: " << pcur_time << " pcur_time: " << pcur_time_evac << " " << problem.numEdges << '\n';
  }
  for (int i_r1=0; i_r1 <= problem.numEdges; i_r1++)
     {
      (pcur_sol_sr_double)[i_r1] = flows[problem.edgeIndices[i_r1].first][problem.edgeIndices[i_r1].second];
     } 
  //return optimalTotalEvacuationTime;
}

bool SelfishRoutingSolver::checkFlowConservation() const
{
  // check flow
  bool error = false;
  for (int node=0; node<problem.numNodes; ++node)
  {
    double inFlow = 0;
    double outFlow = 0;
    for (int otherNode=0; otherNode<problem.numNodes; ++otherNode)
    {
      if (node != otherNode)
      {
        inFlow += flows[otherNode][node];
        outFlow += flows[node][otherNode];
      }
    }

    double netFlow = inFlow - outFlow;
    // demand node, exit nodes, other nodes
    if (problem.demandNodes.find(node) != problem.demandNodes.end())
    {
      if ((netFlow > (-1 * problem.nodeDemands[node] + 0.1)) || (netFlow < (-1 * problem.nodeDemands[node] - 0.1)))
      {
        error = true;
        std::cout << "ERROR at node: " << node << " inFlow: " << inFlow << " outFlow: " << outFlow << std::endl;
      }
    }
    else if (problem.exitNodes.find(node) != problem.exitNodes.end())
    {
      if (netFlow < -0.1)
      {
        error = true;
        std::cout << "ERROR at node: " << node << " inFlow: " << inFlow << " outFlow: " << outFlow << std::endl;
      }
    }
    else
    {
      if ((netFlow > 0.1) || netFlow < -0.1)
      {
        error = true;
        std::cout << "ERROR at node: " << node << " inFlow: " << inFlow << " outFlow: " << outFlow << std::endl;
      }
    }
  }

  return error;
}
