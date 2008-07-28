#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include "SVMSupport.hpp"
#include <ga/make_ga.h>
#include <eo>
#include <es.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <eo/apply.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef eoReal<eoMinimizingFitness> Indi;

double
safe_svm_evaluator(GRNModel& m, const std::list<std::string>& testingSet)
{
  int pipes[2];
  pid_t pid;

  pipe(pipes);
  if ((pid = fork()) == 0)
  {
    close(pipes[0]);

    m.trainSVMs();
    double result = m.testSVMs(testingSet);
    write(pipes[1], &result, sizeof(double));

    exit(0);
  }
  else
  {
    close(pipes[1]);

    struct timeval tv;
    // Give it 10 minutes...
    tv.tv_sec = 600;
    tv.tv_usec = 0;
    fd_set s;
    FD_ZERO(&s);
    FD_SET(pipes[0], &s);
    double result;
    if (select(pipes[0] + 1, &s, NULL, &s, &tv) <= 0)
    {
      kill(pid, SIGKILL);
      result = 1E60;
    }
    else
    {
      if (read(pipes[0], &result, sizeof(double)) != sizeof(double))
        result = 1E60;
    }
    close(pipes[0]);
    int status;
    wait(&status);

    return result;
  }
  return 1E60;
}

class EvaluateSVMFit
  : public eoEvalFunc<Indi>
{
public:
  EvaluateSVMFit(GRNModel& aM, const std::list<std::string>& aTestingSet,
                 std::list<double>& aLastRun)
    : mM(aM), mTestingSet(aTestingSet), mLastRun(aLastRun)
  {
  }
  
  virtual void operator() (Indi & aIndi)
  {
    if (aIndi.invalid())
    {
      double gamma = exp(aIndi[0]), C = exp(aIndi[1]), p = exp(aIndi[2]);
      mM.setSVMParameters(gamma, C, p);

      double fitness;
      if (!mLastRun.empty())
      {
        fitness = mLastRun.front();
        mLastRun.pop_front();
      }
      else
        fitness = safe_svm_evaluator(mM, mTestingSet);
      aIndi.fitness(fitness);
      std::cout << "SVM Result: log2 gamma (" << aIndi[0] << ") "
                   "log2 C (" << aIndi[1] << ") log2 p (" << aIndi[2]
                << ") Result (" << fitness << ")"
                << std::endl;
    }
  }

private:
  GRNModel& mM;
  const std::list<std::string>& mTestingSet;
  std::list<double>& mLastRun;
};

int
main(int argc, char** argv)
{
  const unsigned int T_SIZE = 3; // size for tournament selection
  const unsigned int VEC_SIZE = 3; // Number of object variables in genotypes
  const unsigned int POP_SIZE = 10; // Size of population
  const unsigned int MAX_GEN = 100; // Maximum number of generation before STOP
  const unsigned int MIN_GEN = 10;  // Minimum number of generation before ...
  const unsigned int STEADY_GEN = 10; // stop after STEADY_GEN gen. without improvement
  const float P_CROSS = 0.8;	// Crossover probability
  const float P_MUT = 0.5;	// mutation probability
  const double EPSILON = 0.01;	// range for real uniform mutation
  double SIGMA = 0.3;       	// std dev. for normal mutation
  // some parameters for chosing among different operators
  const double hypercubeRate = 0.5;     // relative weight for hypercube Xover
  const double segmentRate = 0.5;  // relative weight for segment Xover
  const double uniformMutRate = 0.5;  // relative weight for uniform mutation
  const double detMutRate = 0.5;      // relative weight for det-uniform mutation
  const double normalMutRate = 0.5;   // relative weight for normal mutation
  const unsigned int SEED = 42;	// seed for random number generator

  po::options_description desc;
  std::string matrixdir, model, trainingset, testingset, lastrun;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "The directory set up by SOFT2Matrix")
    ("model", po::value<std::string>(&model),
     "The gene regulatory network model")
    ("trainingset", po::value<std::string>(&trainingset),
     "The list of arrays which have been selected for inclusion in the training set")
    ("testingset", po::value<std::string>(&testingset),
     "The list of arrays which have been selected for inclusion in the testing set")
    ("lastrun", po::value<std::string>(&lastrun),
     "The output file from the last run, to re-use scores from (optional)")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("matrixdir"))
      wrong = "matrixdir";
    else if (!vm.count("model"))
      wrong = "model";
    else if (!vm.count("trainingset"))
      wrong = "trainingset";
    else if (!vm.count("testingset"))
      wrong = "testingset";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }
  
  if (!fs::is_directory(matrixdir))
  {
    std::cout << "Matrix directory doesn't exist."
              << std::endl;
    return 1;
  }

  if (!fs::is_regular(model))
  {
    std::cout << "Model file doesn't exist or not regular file."
              << std::endl;
    return 1;
  }

  if (!fs::is_regular(trainingset))
  {
    std::cout << "Training set file doesn't exist or not regular file."
              << std::endl;
    return 1;
  }

  if (!fs::is_regular(testingset))
  {
    std::cout << "Testing set file doesn't exist or not regular file."
              << std::endl;
    return 1;
  }

  std::list<double> lastRun;
  if (fs::is_regular(lastrun))
  {
    std::ifstream flastrun(lastrun.c_str());

    static const boost::regex prev(".*SVM Result: .* Result \\(([^\\)]+)\\).*");
    while (flastrun.good())
    {
      std::string l;
      std::getline(flastrun, l);
      boost::smatch m;
      if (!boost::regex_match(l, m, prev))
      {
        continue;
      }

      lastRun.push_back(strtod(m[1].str().c_str(), NULL));
    }
  }

  ExpressionMatrixProcessor emp(matrixdir);
  GRNModel m(model, emp, 50);

  std::list<std::string> trainingArrays, testingArrays;
  m.loadArraySet(trainingset, trainingArrays);
  m.loadArraySet(testingset, testingArrays);
  m.loadSVMTrainingData(trainingArrays);

  // We seed it just so we can restart if need be.
  rng.reseed(SEED);

  EvaluateSVMFit eval(m, trainingArrays, lastRun);
  eoUniformGenerator<double> uGen(-15.0, 15.0);
  eoInitFixedLength<Indi> random(VEC_SIZE, uGen);
  eoPop<Indi> pop(POP_SIZE, random);
  apply<Indi>(eval, pop);
  pop.sort();
  std::cout << "Initial Population" << std::endl;
  std::cout << pop;
  eoDetTournamentSelect<Indi> selectOne(T_SIZE);
  eoSelectPerc<Indi> select(selectOne);// by default rate==1
  eoGenerationalReplacement<Indi> replace;
  eoSegmentCrossover<Indi> xoverS;
  eoHypercubeCrossover<Indi> xoverA;
  eoPropCombinedQuadOp<Indi> xover(xoverS, segmentRate);
  xover.add(xoverA, hypercubeRate, true);
  
  eoUniformMutation<Indi>  mutationU(EPSILON);
  eoDetUniformMutation<Indi>  mutationD(EPSILON);
  eoNormalMutation<Indi>  mutationN(SIGMA);
  eoPropCombinedMonOp<Indi> mutation(mutationU, uniformMutRate);
  mutation.add(mutationD, detMutRate);
  mutation.add(mutationN, normalMutRate, true);

  eoGenContinue<Indi> genCont(MAX_GEN);
  eoSteadyFitContinue<Indi> steadyCont(MIN_GEN, STEADY_GEN);
  eoFitContinue<Indi> fitCont(0);
  eoCombinedContinue<Indi> continuator(genCont);
  continuator.add(steadyCont);
  continuator.add(fitCont);

  eoSGATransform<Indi> transform(xover, P_CROSS, mutation, P_MUT);
  eoEasyEA<Indi> gga(continuator, eval, select, transform, replace);
  gga(pop);

  pop.sort();
  std::cout << "Final Population:"
            << std::endl << pop << std::endl;

  return 0;
}
