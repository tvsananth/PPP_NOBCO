#include "selfish_routing.h"
#include "fr_solver.h"

void usage()
{
  std::cout << std::endl;
  std::cout << "USAGE: leblanc_solver nodeFileName edgeFileName SR/FR" << std::endl;
}

int main(int argc, char **argv)
{
  srand(time(NULL));

  if (argc == 5)
  {
    std::string nodeFileName = argv[1];
    std::string edgeFileName = argv[2];
    std::string solveMethod = argv[3];
    std::string outputLogFileName  = argv[4];
    NetworkDesignProblem problem(nodeFileName, edgeFileName, outputLogFileName);
    SelfishRoutingSolver solver(problem, true);
 
    if (solveMethod == "SR")
    {
      std::vector<std::vector<double>> emptyInitialSrFlows;
      std::cout << "solve user optimal" << std::endl;
      solver.solve(problem.initialFrEdges, FunctionType::USER_OPTIMAL, emptyInitialSrFlows);
      std::cout << "solve system optimal" << std::endl;
      solver.solve(problem.initialFrEdges, FunctionType::SYSTEM_OPTIMAL, emptyInitialSrFlows);
      std::cout << "done solving" << std::endl;
    }
    else
    {
      FRSolver frSolver(problem);
      frSolver.solveBranchAndBound();
      std::cout << "done solving" << std::endl;
    }
  }
  else
  {
    usage();
    return -1;
  }

  return 0;
};
