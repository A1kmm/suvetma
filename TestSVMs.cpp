#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include "SVMSupport.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class ResultSaver
{
public:
  ResultSaver(uint32_t nGenes, const std::string& aFile)
    : mnGenes(nGenes)
  {
    mData = new double[nGenes];
    mOutput = fopen(aFile.c_str(), "w");
  }

  ~ResultSaver()
  {
    delete [] mData;
    fclose(mOutput);
  }

  void startRow(uint32_t aArray)
  {
    for (uint32_t i = 0; i < mnGenes; i++)
      mData[i] = std::numeric_limits<double>::quiet_NaN();
  }

  void result(uint32_t aGene, double aResult)
  {
    mData[aGene] = aResult;
  }

  void endRow(uint32_t aArray)
  {
    fwrite(mData, mnGenes * sizeof(double), 1, mOutput);
  }

private:
  double* mData;
  FILE* mOutput;
  uint32_t mnGenes;
};

int
main(int argc, char** argv)
{
  po::options_description desc;
  std::string matrixdir, model, svmdir, testingset, output;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "The directory set up by SOFT2Matrix")
    ("model", po::value<std::string>(&model),
     "The gene regulatory network model")
    ("svmdir", po::value<std::string>(&svmdir),
     "The directory containing the support vector machines")
    ("testingset", po::value<std::string>(&testingset),
     "The list of arrays which have been selected for inclusion in the "
     "testing set")
    ("output", po::value<std::string>(&output),
     "The file to write the per-array data into")
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
    else if (!vm.count("svmdir"))
      wrong = "svmdir";
    else if (!vm.count("testingset"))
      wrong = "testingset";
    else if (!vm.count("output"))
      wrong = "output";
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

  if (!fs::is_directory(svmdir))
  {
    std::cout << "SVM directory doesn't exist."
              << std::endl;
    return 1;
  }

  if (!fs::is_regular(testingset))
  {
    std::cout << "Testing set file doesn't exist or not regular file."
              << std::endl;
    return 1;
  }

  ExpressionMatrixProcessor emp(matrixdir);
  GRNModel m(model, emp);
  ResultSaver rs(emp.getNumGenes(), output);
  m.loadSVMs(svmdir);

  std::list<std::string> testingSet;
  m.loadArraySet(testingset, testingSet);
  
  m.testSVMs(testingSet, rs);
}
