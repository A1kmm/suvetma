/*
    Perform sign tests, a gene at a time.
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
#include <list>
#include <fstream>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class SignTestFits
{
public:
  SignTestFits
  (
   const std::string& aControl,
   const std::string& aModel,
   const std::string& aGeneList
  )
    : mfControl(NULL), mfModel(NULL)
  {
    loadGeneList(aGeneList);

    mfControl = fopen(aControl.c_str(), "r");
    mfModel = fopen(aModel.c_str(), "r");
    mBuffer1 = new double[mNGenes];
    mBuffer2 = new double[mNGenes];
    mControlWorse = new uint32_t[mNGenes];
    mTotal = new uint32_t[mNGenes];
    memset(mControlWorse, 0, sizeof(uint32_t) * mNGenes);
    memset(mTotal, 0, sizeof(uint32_t) * mNGenes);

    processFitData();
    dumpFitData();
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
    if (mControlWorse)
      delete [] mControlWorse;
    if (mTotal)
      delete [] mTotal;
  }

  void
  processFitData()
  {
    while (true)
    {
      if (fread(mBuffer1, sizeof(double), mNGenes, mfControl) !=
          mNGenes)
        break;

      if (fread(mBuffer2, sizeof(double), mNGenes, mfModel) !=
          mNGenes)
        break;

      updateSignStatistics();
    }
  }

  void
  dumpFitData()
  {
    std::cout << "\"gene\",\"total\",\"control.better\"" << std::endl;
    uint32_t v = 0;
    for (std::list<std::string>::iterator i = mGenes.begin();
         i != mGenes.end();
         i++, v++)
      std::cout << "\"" << *i << "\",\"" << mTotal[v] << "\",\""
                << mControlWorse[v] << "\"" << std::endl;
  }

private:
  FILE * mfControl, * mfModel;
  double * mBuffer1, * mBuffer2;
  uint32_t * mControlWorse, * mTotal;
  uint32_t mNGenes;
  std::list<std::string> mGenes;

  void
  updateSignStatistics()
  {
    for (uint32_t i = 0; i < mNGenes; i++)
    {
      if (finite(mBuffer1[i]) && finite(mBuffer2[i]))
      {
        mTotal[i]++;
        mControlWorse[i] += (mBuffer1[i] >= mBuffer2[i]);
      }
    }
  }

  void
  loadGeneList(const std::string& aFileName)
  {
    std::ifstream gl(aFileName.c_str());

    while (gl.good())
    {
      std::string l;
      std::getline(gl, l);

      if (!gl.good())
        break;

      mGenes.push_back(l);
    }

    // The + 1 is a workaround for a bug in the other tools due to the trailing newline.
    mNGenes = mGenes.size();
  }
};

int
main(int argc, char** argv)
{
  po::options_description desc;
  std::string controlMatrix, modelMatrix, geneList;

  desc.add_options()
    ("controlmatrix", po::value<std::string>(&controlMatrix),
     "The error matrix for the model expected to fit poorly")
    ("modelmatrix", po::value<std::string>(&modelMatrix),
     "The error matrix for the model expected to fit well")
    ("genelist", po::value<std::string>(&geneList),
     "The file containing the list of genes")
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
    else if (!vm.count("genelist"))
      wrong = "genelist";
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

  if (!fs::is_regular(geneList))
  {
    std::cout << "Gene list not found." << std::endl;
  }

  SignTestFits stf(controlMatrix, modelMatrix, geneList);
}
