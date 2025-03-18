#include "fr_solver.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <iostream>
#include <fstream>

void BBNode::print() const
{
  std::cout << "node: " << counter << std::endl;
  std::cout << "lower bound: " << lowerBound << std::endl;
  std::cout << "upper bound: " << upperBound << std::endl;
  std::cout << "index to branch: " << frIndexToBranchOn << std::endl;

  std::cout << "edges to use: ";
  for (auto indexEdge : edgesToUse)
  {
    std::cout << indexEdge.first << ": ";
    for (auto edge : indexEdge.second)
    {
      std::cout << edge.first << "->" << edge.second << ", ";
    }
  }
  std::cout << std::endl;
 
  std::cout << "edges to not use: ";
  for (auto indexEdge : edgesToNotUse)
  {
    std::cout << indexEdge.first << ": ";
    for (auto edge : indexEdge.second)
    {
      std::cout << edge.first << "->" << edge.second << ", ";
    }
  }
  std::cout << std::endl;
}

FRSolver::FRSolver(NetworkDesignProblem _problem) : problem(_problem), selfishRoutingSolver(_problem, false), nodeCounter(1), nodesProcessed(0), upperBound(std::numeric_limits<double>::max()), lowerBound(0), extraMilliseconds(0) {}

bool FRSolver::solveBranchAndBound()
{
  auto startBranchAndBound = std::chrono::high_resolution_clock::now();

  // Create root node
  std::map<int,std::vector<std::pair<int,int> > > emptyMap;
  std::set<int> emptySet;
  std::vector<std::vector<double>> emptyInitialSrFlows;
  BBNode rootNode(emptyMap, emptySet, emptyMap, 0, emptyInitialSrFlows, lowerBound, upperBound, nodeCounter);
  nodeQueue.push(rootNode);
  std::cout << "STATS lb[" << lowerBound << "] ub[" << upperBound << "] time[0] numNodes[" << nodesProcessed << "]" << std::endl;

  std::ofstream file_out_bnb(problem.ologFileName, std::ios::out | std::ios::binary);
  double time_duration_count = 0;
  // Process until node queue empty or problem solved
  std::set<int> openNodes;
  openNodes.insert(nodeCounter);
  std::map<int,double> openNodeLowerBounds;
  while (!nodeQueue.empty() && (lowerBound != upperBound))
  {
    BBNode node = nodeQueue.top();
    nodeQueue.pop();
    std::cout << "Processing Node: " << node.counter << std::endl;
    ++nodesProcessed;
    node.print();

    // Calculate UB with Primal Heuristic
    std::pair<int,int> edgeToBranchOn;
    double nodeUpperBound = primalHeuristic(node, edgeToBranchOn);

    // TODO(akarahal) Want to store the solution, so that we can 
    // Also use UB solution to guide children nodes by storing info
    // Loop over FR edges, add 0.5 times the distance to a dictionary
    // Update node to store this information, so creating next nodes passes this dict
    // Update shortest path (or before/after calculation) to add/subtract the amounts
    if ((nodeUpperBound < upperBound) && nodeUpperBound > 0)
    {
      upperBound = nodeUpperBound;
    }

    if (edgeToBranchOn == std::make_pair(-1,-1))
    {
      std::cout << "infeasible no new node created" << std::endl;
    }
    else
    {
      // evaluate lower bound and check
      //node.lowerBound = calculateLowerBound(node);
      openNodeLowerBounds[node.counter] = node.lowerBound;
      if (node.counter == 1)
      {
        lowerBound = node.lowerBound;
      }

      if (node.lowerBound > upperBound)
      {
        std::cout << "lower bound: " << node.lowerBound << std::endl;
        std::cout << "lower bound greater than upper bound, no new node created" << std::endl;
      }
      else
      {
        // branch to create new nodes
        int nextFrIndexToBranchOn = (node.frIndexToBranchOn == problem.frDemandNodes.size() - 1) ? 0 : node.frIndexToBranchOn + 1;
        int numChecked = 0;
        bool allDemandComplete = false;
        while (node.completedDemands.find(nextFrIndexToBranchOn) != node.completedDemands.end())
        {
          nextFrIndexToBranchOn = (nextFrIndexToBranchOn == problem.frDemandNodes.size() - 1) ? 0 : nextFrIndexToBranchOn + 1;

          numChecked = numChecked + 1;
          if (numChecked > problem.frDemandNodes.size())
          {
            allDemandComplete = true;
            break;
          }
        }

        if (!allDemandComplete)
        {
          // One node uses the edge
          nodeCounter = nodeCounter + 1;
          std::map<int,std::vector<std::pair<int,int> > > updatedEdgesToUse1(node.edgesToUse);
          std::set<int> updatedCompletedDemands1(node.completedDemands);
          std::map<int,std::vector<std::pair<int,int> > > updatedEdgesToNotUse1(node.edgesToNotUse);
          if (updatedEdgesToUse1.find(node.frIndexToBranchOn) != updatedEdgesToUse1.end())
          {
            updatedEdgesToUse1[node.frIndexToBranchOn].push_back(edgeToBranchOn);
          }
          else
          {
            std::vector<std::pair<int,int> > initializeVector;
            initializeVector.push_back(edgeToBranchOn);
            updatedEdgesToUse1[node.frIndexToBranchOn] = initializeVector;
          }
          // if we have branched all the way to an FR entry node, mark as completed
          if (problem.frEntryNodes.find(edgeToBranchOn.first) != problem.frEntryNodes.end())
          {
            updatedCompletedDemands1.insert(node.frIndexToBranchOn);
          }
          BBNode nodeWithBranchingEdge(updatedEdgesToUse1, updatedCompletedDemands1, updatedEdgesToNotUse1, nextFrIndexToBranchOn, selfishRoutingSolver.getFlows(), node.lowerBound, nodeUpperBound, nodeCounter);
          nodeQueue.push(nodeWithBranchingEdge);
          openNodes.insert(nodeCounter);
          openNodeLowerBounds[nodeCounter] = node.lowerBound;

          // One node does not use the edge
          nodeCounter = nodeCounter + 1;
          std::map<int,std::vector<std::pair<int,int> > > updatedEdgesToUse2(node.edgesToUse);
          std::set<int> updatedCompletedDemands2(node.completedDemands);
          std::map<int,std::vector<std::pair<int,int> > > updatedEdgesToNotUse2(node.edgesToNotUse);
          if (updatedEdgesToNotUse2.find(node.frIndexToBranchOn) != updatedEdgesToNotUse2.end())
          {
            updatedEdgesToNotUse2[node.frIndexToBranchOn].push_back(edgeToBranchOn);
          }
          else
          {
            std::vector<std::pair<int,int> > initializeVector;
            initializeVector.push_back(edgeToBranchOn);
            updatedEdgesToNotUse2[node.frIndexToBranchOn] = initializeVector;
          }
          BBNode nodeWithoutBranchingEdge(updatedEdgesToUse2, updatedCompletedDemands2, updatedEdgesToNotUse2, nextFrIndexToBranchOn, selfishRoutingSolver.getFlows(), node.lowerBound, nodeUpperBound, nodeCounter);
          nodeQueue.push(nodeWithoutBranchingEdge);
          openNodes.insert(nodeCounter);
          openNodeLowerBounds[nodeCounter] = node.lowerBound;
        }
        else
        {
          std::cout << "all demands completed, no new node created" << std::endl;
        }
      }
    }

    // Pop to process the next node
    openNodes.erase(node.counter);

    // Print stats
    auto currentTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startBranchAndBound);
    lowerBound = std::numeric_limits<double>::max(); 
    for (auto nodeCount : openNodes)
    {
      //std::cout << "open node: " << nodeCount << " has lb: " << openNodeLowerBounds[nodeCount] << std::endl;
      lowerBound = std::min(lowerBound, openNodeLowerBounds[nodeCount]);
    }
    std::cout << "STATS lb[" << lowerBound << "] ub[" << upperBound << "] time[" << duration.count() - extraMilliseconds / 1000.0 << "] numNodes[" << nodesProcessed << "]" << std::endl;

    time_duration_count = duration.count();
    file_out_bnb.write((char *) &lowerBound, sizeof(double));
    file_out_bnb.write((char *) &upperBound, sizeof(double));
    file_out_bnb.write((char *) &time_duration_count, sizeof(double));

  }

  file_out_bnb.close();

  auto endBranchAndBound = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(endBranchAndBound - startBranchAndBound);

  std::cout << "branch and bound took: " << duration.count() - extraMilliseconds / 1000.0 << " seconds" << std::endl;
  return true;
}

double FRSolver::primalHeuristic(const BBNode node, std::pair<int,int>& branchEdge)
{
  // include all edges to use
  std::set<std::pair<int,int> > frEdges;
  for (auto edgeToUse : node.edgesToUse)
  {
    for (auto edge : edgeToUse.second)
    {
      int min = std::min(edge.first, edge.second);
      int max = std::max(edge.first, edge.second);
      auto forwardEdge = std::make_pair(min, max);
      frEdges.insert(forwardEdge);
    }
  }

  // For each FR demand node
  for (int frIndex=0; frIndex<problem.frDemandNodes.size(); ++frIndex)
  {
    int frDemandNode = problem.frDemandNodes[frIndex];

    // Include branch-node-(+)-decisions
    int updatedFrDemandNode = frDemandNode;
    if (node.edgesToUse.find(frIndex) != node.edgesToUse.end())
    {
      updatedFrDemandNode = node.edgesToUse.at(frIndex).back().first;
    }
 
    // Loop over all FR entry points, to find shortest one
    std::set<int> frDemandNodeAsSet;
    frDemandNodeAsSet.insert(updatedFrDemandNode);
    std::vector<int> bestShortestPathByNodes;
    double bestShortestPathDistance = std::numeric_limits<double>::max(); 
    for (int frEntryNode : problem.frEntryNodes)
    {
      // Get shortest path without branch-node-(-)-decisions
      std::vector<int> shortestPathByNodes;
      std::vector<std::pair<int,int> > frEdgesToNotUse;
      if (node.edgesToNotUse.find(frIndex) != node.edgesToNotUse.end())
      {
        frEdgesToNotUse = node.edgesToNotUse.at(frIndex);
      }
      double shortestPathDistance = selfishRoutingSolver.calculateShortestPathFromNodeToSetOfNodes(frEntryNode, frDemandNodeAsSet, shortestPathByNodes, frEdgesToNotUse, DistanceType::FREE_TRAVEL, frEdges);

      // Get shortest route length
      //std::cout << "fr entry node: " << frEntryNode << " gives distance: " << shortestPathDistance << std::endl;
      if (shortestPathDistance < bestShortestPathDistance)
      {
        bestShortestPathDistance = shortestPathDistance;
        bestShortestPathByNodes = shortestPathByNodes;
      }
    }

    // fully specified path is the best but shouldn't be branched on
    if (bestShortestPathByNodes.empty())
    {
      // handle disconnected case (using -1 sentinel value)
      if (bestShortestPathDistance >= std::numeric_limits<double>::max())
      {
        branchEdge = std::make_pair(-1,-1);
        std::cout << "disconnected fr demand node " << updatedFrDemandNode << std::endl; 
        return -1;
      }

      if (node.frIndexToBranchOn == frIndex)
      {
        branchEdge = node.edgesToUse.at(frIndex).back();
      }
    }
    else
    {
      // branching rule: choose for current index the first edge from demand to exit
      for (int index=0; index<bestShortestPathByNodes.size()-1; ++index)
      {
        auto edge = std::make_pair(bestShortestPathByNodes[index],bestShortestPathByNodes[index+1]);
        int min = std::min(edge.first, edge.second);
        int max = std::max(edge.first, edge.second);
        auto forwardEdge = std::make_pair(min, max);
        frEdges.insert(forwardEdge);
        if ((index == bestShortestPathByNodes.size()-2) && (node.frIndexToBranchOn == frIndex))
        {
          branchEdge = edge;
        }
      }
    }

    std::cout << "fr route for " << frDemandNode << std::endl;
    for (auto node : bestShortestPathByNodes)
    {
      std::cout << node << ",";
    }
    std::cout << std::endl;
    //std::cout << "branch edge:" << branchEdge.first << "," << branchEdge.second << std::endl;
  }
 
  /*
  // Using these FR choices, run selfish routing
  double systemOptimalValue = selfishRoutingSolver.solve(frEdges, FunctionType::SYSTEM_OPTIMAL, node.initialSrFlows);
  std::cout << "system optimal value: " << systemOptimalValue << std::endl;

  // Checking FIR, so need both solutions if it might be best
  if ((systemOptimalValue < upperBound) && systemOptimalValue > 0)
  {
    // compute user optimal
    auto startTime = std::chrono::high_resolution_clock::now();
    double userOptimalValue = selfishRoutingSolver.solve(frEdges, FunctionType::USER_OPTIMAL, node.initialSrFlows);
    std::cout << "user optimal value: " << userOptimalValue << std::endl;

    // go back to system optimal if needed to be stored
    systemOptimalValue = selfishRoutingSolver.solve(frEdges, FunctionType::SYSTEM_OPTIMAL, node.initialSrFlows);
    std::cout << "system optimal value: " << systemOptimalValue << std::endl;

    // remove this time from clock
    auto endTime = std::chrono::high_resolution_clock::now();
    auto userDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "user duration: " << userDuration.count() << std::endl;
    extraMilliseconds += userDuration.count();
  }
  */

  double optimalValue = selfishRoutingSolver.solve(frEdges, FunctionType::USER_OPTIMAL, node.initialSrFlows);
  std::cout << "optimal value: " << optimalValue << std::endl;
  return optimalValue;

  //return systemOptimalValue;
}

double FRSolver::calculateLowerBound(const BBNode node)
{
  // Update the graph with branch-node-(+)-decisions
  // Like Leblanc, use optimal non-selfish routing objective function as the lower bound
  std::set<std::pair<int,int> > frEdges;
  for (auto edgeToUse : node.edgesToUse)
  {
    for (auto edge : edgeToUse.second)
    {
      int min = std::min(edge.first, edge.second);
      int max = std::max(edge.first, edge.second);
      auto forwardEdge = std::make_pair(min, max);
      frEdges.insert(forwardEdge);
    }
  }

  std::vector<std::vector<double>> emptyInitialFlows;
  double nodeLowerBound = selfishRoutingSolver.solve(frEdges, FunctionType::SYSTEM_OPTIMAL, emptyInitialFlows);

  return nodeLowerBound;
}
