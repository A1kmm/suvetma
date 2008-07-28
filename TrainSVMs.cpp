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
  std::string matrixdir, model, svmdir, trainingset;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "The directory set up by SOFT2Matrix")
    ("model", po::value<std::string>(&model),
     "The gene regulatory network model")
    ("svmdir", po::value<std::string>(&svmdir),
     "The directory containing the support vector machines")
    ("trainingset", po::value<std::string>(&trainingset),
     "The list of arrays which have been selected for inclusion in the training set")
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
    else if (!vm.count("trainingset"))
      wrong = "trainingset";
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
    fs::path p(svmdir);
    try
    {
      fs::create_directory(p);
    }
    catch (std::exception e)
    {
      std::cout << "SVM directory doesn't exist and couldn't be created: "
                << e.what() << std::endl;
      return 1;
    }
  }

  if (!fs::is_regular(trainingset))
  {
    std::cout << "Training set file doesn't exist or not regular file."
              << std::endl;
    return 1;
  }

  ExpressionMatrixProcessor emp(matrixdir);
  GRNModel m(model, emp, 50);

  std::list<std::string> trainingArrays;
  m.loadArraySet(trainingset, trainingArrays);
  m.loadSVMTrainingData(trainingArrays);
  m.trainSVMs();
  m.saveSVMs(svmdir);
  
  return 0;
}
