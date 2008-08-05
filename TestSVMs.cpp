#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include "SVMSupport.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

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
  m.loadSVMs(svmdir);

  std::list<std::string> testingSet;
  m.loadArraySet(testingset, testingSet);
  
  std::vector<std::pair<double, uint32_t> > results;
  m.testSVMs(testingSet, results);

  std::ofstream ofile(output.c_str());

  std::vector<std::pair<double, uint32_t> >::iterator i;
  std::list<std::string>::iterator j;
  for (i = results.begin(), j = testingSet.begin();
       i != results.end();
       i++, j++)
    ofile << *j << "\t" << (*i).first << "\t" << (*i).second
          << std::endl;
}
