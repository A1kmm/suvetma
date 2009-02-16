/*
    Train the support vector machines.
    Copyright (C) 2008-2009  Andrew Miller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
  double loggamma, logC, nu;
  bool dontReplace;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "The directory set up by SOFT2Matrix")
    ("model", po::value<std::string>(&model),
     "The gene regulatory network model")
    ("svmdir", po::value<std::string>(&svmdir),
     "The directory containing the support vector machines")
    ("trainingset", po::value<std::string>(&trainingset),
     "The list of arrays which have been selected for inclusion in the training set")
    ("gamma", po::value<double>(&loggamma),
     "The value of the RBF parameter gamma, as a base-e logarithm of the value")
    ("C", po::value<double>(&logC),
     "The value of the SVM parameter C, as a base-e logarithm of the value")
    ("nu", po::value<double>(&nu),
     "The value of the SVM parameter nu, as a base-e logarithm of the value")
    ("dont-replace", "Indicates that existing SVMs shouldn't be replaced")
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

  dontReplace = (vm.count("dont-replace") != 0);

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
  GRNModel m(model, emp);

  std::list<std::string> trainingArrays;
  m.loadArraySet(trainingset, trainingArrays);
  m.setSVMParameters(exp(loggamma), exp(logC), nu);
  m.loadSVMTrainingData(trainingArrays);
  m.trainAndSaveSVMs(svmdir, dontReplace);
  
  return 0;
}
