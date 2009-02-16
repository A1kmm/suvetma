/*
    SVM tool support classes
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

#include "SVMSupport.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/tokenizer.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

ExpressionMatrixProcessor::ExpressionMatrixProcessor
(
 const std::string& aMatrixDir
)
  : mDataFile(NULL), mRow(NULL)
{
  fs::path md(aMatrixDir);

  {
    fs::path arraylist(md);
    arraylist /= "arrays";

    std::ifstream af(arraylist.string().c_str());
    uint32_t index = 0;
    
    while (af.good())
    {
      std::string array;
      std::getline(af, array);

      if (!af.good())
        break;
      
      mArrayIndices.insert(std::pair<std::string, uint32_t>(array, index++));
    }

    mnArrays = index;
  }

  {
    fs::path genelist(md);
    genelist /= "genes";

    std::ifstream gf(genelist.string().c_str());
    uint32_t index = 0;

    while (gf.good())
    {
      std::string gene;
      std::getline(gf, gene);

      if (!gf.good())
        break;

      mGeneIndices.insert(std::pair<std::string, uint32_t>(gene, index++));
    }

    mnGenes = index;
  }

  {
    fs::path datafile(md);
    datafile /= "data";

    mDataFile = fopen(datafile.string().c_str(), "r");
  }

  mRow = new double[mnGenes];
}

ExpressionMatrixProcessor::~ExpressionMatrixProcessor()
{
  if (mDataFile != NULL)
    fclose(mDataFile);
  if (mRow)
    delete [] mRow;
}

uint32_t
ExpressionMatrixProcessor::getIndexOfGene(const std::string& aGene)
{
  return mGeneIndices[aGene];
}

uint32_t
ExpressionMatrixProcessor::getIndexOfArray(const std::string& aArray)
{
  return mArrayIndices[aArray];
}

void
ExpressionMatrixProcessor::setArray
(
 uint32_t aArray
)
{
  fseek(mDataFile, aArray * mnGenes * sizeof(double), SEEK_SET);
  fread(mRow, mnGenes * sizeof(double), 1, mDataFile);
}

double
ExpressionMatrixProcessor::getDataPoint
(
 uint32_t aGene
)
{
  return mRow[aGene];
}

SupportVectorMachine::SupportVectorMachine
(
 ExpressionMatrixProcessor& aEMP,
 const std::string& aRegulatedGene,
 uint32_t aNumRegulators
)
  : mEMP(aEMP), mRegulatedGene(aEMP.getIndexOfGene(aRegulatedGene)),
    mNumRegulators(aNumRegulators), mModel(NULL), mGamma(0.1), mC(0.1),
    mNu(0.1), mTestNodes(NULL), mRegulatedGeneName(aRegulatedGene)
{
  mRegulatingGenes.reserve(aNumRegulators);
  mProblem.y = NULL;
  mProblem.x = NULL;

  mTestNodes = new svm_node[mNumRegulators + 1];
}

SupportVectorMachine::~SupportVectorMachine()
{
  if (mModel)
    svm_destroy_model(mModel);

  if (mProblem.y)
    delete [] mProblem.y;

  if (mProblem.x)
  {
    if (mProblem.x[0])
      delete [] mProblem.x[0];
    delete [] mProblem.x;
  }

  if (mTestNodes)
  {
    delete mTestNodes;
  }
}

void
SupportVectorMachine::addRegulatingGene(uint32_t aRegulatingGene)
{
  mRegulatingGenes.push_back(aRegulatingGene);
}

void
SupportVectorMachine::setupProblem(uint32_t aSize)
{
  // Handle this case specially so we can make the simplifying assumption that
  // mNumRegulators > 0 later...
  if (mNumRegulators == 0)
  {
    mNumFinite = 0;
    mSum = 0.0;
    return;
  }

  mProblem.l = 0;
  mProblem.y = new double[aSize];
  mProblem.x = new svm_node*[aSize];

  double* yp = mProblem.y;
  svm_node** xp = mProblem.x;
  svm_node* p = new svm_node[aSize * (mNumRegulators + 1)];

  for (uint32_t i = 0; i < aSize; i++)
  {
    xp[i] = p;
    p += (mNumRegulators + 1);
  }

  mYp = mProblem.y;
  mXp = mProblem.x;
}

void
SupportVectorMachine::loadTrainingRow()
{
  if (mNumRegulators == 0)
  {
    double rgl = mEMP.getDataPoint(mRegulatedGene);
    if (!isfinite(rgl))
      return;
    
    mNumFinite++;
    mSum += rgl;

    return;
  }

  bool containsNaNs = false;

  for (std::vector<uint32_t>::iterator j = mRegulatingGenes.begin();
       j != mRegulatingGenes.end(); j++)
    if (!isfinite(mEMP.getDataPoint(*j)))
      containsNaNs = true;

  // We just ignore the whole array if there are NaNs...
  if (containsNaNs)
    return;

  if (!isfinite(mEMP.getDataPoint(mRegulatedGene)))
    return;

  mProblem.l++;

  *mYp++ = mEMP.getDataPoint(mRegulatedGene);
  svm_node* p;
  p = *mXp++;

  uint32_t k = 1;
  for (std::vector<uint32_t>::iterator j = mRegulatingGenes.begin();
       j != mRegulatingGenes.end(); p++, j++, k++)
  {
    p->index = k;
    p->value = mEMP.getDataPoint(*j);
  }
  p->index = -1;
}

void
SupportVectorMachine::train()
{
  // Handle this case specially so we can make the simplifying assumption that
  // mNumRegulators > 0 later...
  if (mNumRegulators == 0)
  {
    mAverage = mSum / mNumFinite;
    return;
  }

  mParameter.svm_type = NU_SVR;
  mParameter.kernel_type = RBF;
  mParameter.degree = 3;
  mParameter.gamma = mGamma;
  mParameter.coef0 = 0;
  mParameter.p = 0;
  mParameter.cache_size = 100;
  mParameter.C = mC;
  mParameter.eps = 1E-3;
  mParameter.nu = mNu;
  mParameter.shrinking = 1;
  mParameter.probability = 0;
  mParameter.nr_weight = 0;
  mParameter.weight_label = NULL;
  mParameter.weight = NULL;

  if (mModel != NULL)
    svm_destroy_model(mModel);

  mModel = svm_train(&mProblem, &mParameter);
}

void
GRNModel::trainAndSaveSVMs(const std::string& aSVMDir, bool aDontReplace)
{
  fs::path dir(aSVMDir);
  
  for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
       i != mSVMs.end();
       i++)
  {
    (*i)->train();
    fs::path targ(dir);
    if (fs::is_regular(targ) && aDontReplace)
      continue;
    targ /= (*i)->getRegulatedGeneName();
    (*i)->save(targ.string());
  }
}

void
SupportVectorMachine::save(const std::string& aFilename)
{
  svm_save_model(aFilename.c_str(), mModel);
}

void
SupportVectorMachine::load(const std::string& aFilename)
{
  if (mModel)
    svm_destroy_model(mModel);

  mModel = svm_load_model(aFilename.c_str());

  assert(mModel);
}

double
SupportVectorMachine::testOnRow()
{
  // We just ignore the whole array if there are NaNs...
  svm_node* p = mTestNodes;
  uint32_t k = 1;
  for (std::vector<uint32_t>::iterator j = mRegulatingGenes.begin();
       j != mRegulatingGenes.end(); j++)
  {
    double v(mEMP.getDataPoint(*j));
    if (!isfinite(v))
      return std::numeric_limits<double>::quiet_NaN();
    else
    {
      p->index = k++;
      p->value = v;
      p++;
    }
  }

  p->index = -1;

  double answer = mEMP.getDataPoint(mRegulatedGene);
  if (!isfinite(answer))
    return std::numeric_limits<double>::quiet_NaN();

  double x = (svm_predict(mModel, mTestNodes) - answer);
  return x * x;
}

GRNModel::GRNModel(const std::string& aModel,
                   ExpressionMatrixProcessor& aEMP,
                   uint32_t aGeneLimit)
  : mEMP(aEMP)
{
  std::ifstream m(aModel.c_str());

  bool unlimitedGenes = (aGeneLimit == 0);

  std::string l;
  std::getline(m, l);
  if (l != "VERTICES")
  {
    std::cerr << "Expected GRN model to start with VERTICES. Not loaded."
              << std::endl;
    return;
  }
  
  const static boost::regex vertexline("VERTEX ([0-9]+) (.*)");

  while (m.good())
  {
    std::getline(m, l);
    if (l == "ENDVERTICES")
      break;

    boost::smatch res;
    if (!boost::regex_match(l, res, vertexline))
      continue;

    uint32_t vno = strtoul(res[1].str().c_str(), NULL, 10);
    mHGNCByVertex[vno] = res[2];
  }

  const static boost::regex edgeline("EDGES ([0-9]+) \\((.*)\\)");

  while (m.good() && (unlimitedGenes || aGeneLimit-- > 0))
  {
    std::getline(m, l);
    boost::smatch res;
    if (!boost::regex_match(l, res, edgeline))
    {
      aGeneLimit++;
      continue;
    }

    uint32_t g(strtoul(res[1].str().c_str(), NULL, 10));
    std::string targGene(mHGNCByVertex[g]);

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(" ");
    std::string rgl(res[2]);
    tokenizer t(rgl, sep);
    std::vector<uint32_t> regulators;
    for (tokenizer::iterator i = t.begin(); i != t.end(); i++)
    {
      std::string g(*i);

      regulators.push_back(aEMP.getIndexOfGene
                           (mHGNCByVertex[strtoul(g.c_str(), NULL, 10)]));
    }
    
    if (targGene == "")
    {
      aGeneLimit++;
      continue;
    }

    SupportVectorMachine* svm(new SupportVectorMachine(aEMP, targGene,
                                                       regulators.size()));
    for (std::vector<uint32_t>::iterator i = regulators.begin();
         i != regulators.end();
         i++)
      svm->addRegulatingGene(*i);

    mSVMs.push_back(svm);
  }
}

GRNModel::~GRNModel()
{
  for (
       std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
       i != mSVMs.end();
       i++
      )
    delete *i;
}

void
GRNModel::saveSVMs(const std::string& aFilename)
{
  fs::path dir(aFilename);

  for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
       i != mSVMs.end();
       i++)
  {
    fs::path targ(dir);
    targ /= (*i)->getRegulatedGeneName();
    (*i)->save(targ.string());
  }
}

void
GRNModel::loadSVMs(const std::string& aFilename)
{
  fs::path dir(aFilename);

  for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
       i != mSVMs.end();
       i++)
  {
    fs::path targ(dir);
    targ /= (*i)->getRegulatedGeneName();
    (*i)->load(targ.string());
  }
}
