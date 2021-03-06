/*
    Sign test the fits.
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
#include <string>
#include <cstdio>
#include <sys/types.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class SignTestFits
{
public:
  SignTestFits
  (
   const std::string& aControl,
   const std::string& aModel
  )
    : mfControl(NULL), mfModel(NULL), mnGreater(0), mnTotal(0)
  {
    mfControl = fopen(aControl.c_str(), "r");
    mfModel = fopen(aModel.c_str(), "r");
    mBuffer1 = new double[kBufferSize];
    mBuffer2 = new double[kBufferSize];

    processFitData();

    std::cout << "Of " << mnTotal << " trials, the control had greater error "
              << "in " << mnGreater << std::endl;
  }

  ~SignTestFits()
  {
    if (mfControl)
      fclose(mfControl);
    if (mfModel)
      fclose(mfModel);
    if (mBuffer1)
      delete [] mBuffer1;
    if (mBuffer2)
      delete [] mBuffer2;
  }

  void
  processFitData()
  {
    do
    {
      mBufferUtilisation = fread(mBuffer1, sizeof(double), kBufferSize, mfControl);
      if (fread(mBuffer2, sizeof(double), mBufferUtilisation, mfModel) != mBufferUtilisation)
      {
        std::cerr << "Warning: read mismatch" << std::endl;
      }
      updateSignStatistics();
    }
    while (mBufferUtilisation == kBufferSize);
  }

private:
  FILE * mfControl, * mfModel;
  double * mBuffer1, * mBuffer2;
  static const uint32_t kBufferSize = 512000;
  uint32_t mBufferUtilisation;
  uint32_t mnGreater, mnTotal;

  void
  updateSignStatistics()
  {
    for (uint32_t i = 0; i < mBufferUtilisation; i++)
    {
      if (isfinite(mBuffer1[i]) && isfinite(mBuffer2[i]) && mBuffer1[i] != mBuffer2[i])
      {
        mnGreater += (mBuffer1[i] >= mBuffer2[i]);
        mnTotal++;
      }
    }
  }
};

int
main(int argc, char** argv)
{
  po::options_description desc;
  std::string controlMatrix, modelMatrix;

  desc.add_options()
    ("controlmatrix", po::value<std::string>(&controlMatrix),
     "The error matrix for the model expected to fit poorly")
    ("modelmatrix", po::value<std::string>(&modelMatrix),
     "The error matrix for the model expected to fit well")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("controlmatrix"))
      wrong = "controlmatrix";
    else if (!vm.count("modelmatrix"))
      wrong = "modelmatrix";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!fs::is_regular(controlMatrix))
  {
    std::cout << "Control matrix file not found." << std::endl;
    return 1;
  }

  if (!fs::is_regular(modelMatrix))
  {
    std::cout << "Model matrix file not found." << std::endl;
    return 1;
  }

  SignTestFits stf(controlMatrix, modelMatrix);
}
