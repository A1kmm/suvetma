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

#include <string>
#include <inttypes.h>
#include <map>
#include <vector>
#include <list>
#include <svm.h>
#include <fstream>
#include <math.h>

class ExpressionMatrixProcessor
{
public:
  ExpressionMatrixProcessor(const std::string& aMatrixDir);
  ~ExpressionMatrixProcessor();

  uint32_t getIndexOfGene(const std::string& aGene);
  uint32_t getIndexOfArray(const std::string& aArray);
  void setArray(uint32_t aArray);
  double getDataPoint(uint32_t aGene);
  uint32_t getNumGenes() const { return mnGenes; }
  uint32_t getNumArrays() const { return mnArrays; }

private:
  std::map<std::string, uint32_t> mArrayIndices, mGeneIndices;
  uint32_t mnGenes, mnArrays;
  FILE* mDataFile;
  double* mRow;
};

class SupportVectorMachine
{
public:
  SupportVectorMachine(ExpressionMatrixProcessor& aEMP,
                       const std::string& aRegulatedGene,
                       uint32_t aNumRegulators);
  ~SupportVectorMachine();

  void addRegulatingGene(uint32_t aRegulatingGene);

  void setupProblem(uint32_t aSize);
  void loadTrainingRow();

  void train();
  double testOnRow();
  void save(const std::string& aFilename);
  void load(const std::string& aFilename);

  const std::string& getRegulatedGeneName()
  {
    return mRegulatedGeneName;
  }

  uint32_t getRegulatedGene()
  {
    return mRegulatedGene;
  }

  void
  setParameters(double aGamma, double aC, double aNu)
  {
    mGamma = aGamma;
    mC = aC;
    mNu = aNu;
  }

private:
  ExpressionMatrixProcessor& mEMP;
  std::vector<uint32_t> mRegulatingGenes;
  uint32_t mRegulatedGene, mNumRegulators;
  struct svm_model* mModel;
  struct svm_problem mProblem;
  struct svm_parameter mParameter;
  double mAverage, mSum, mGamma, mC, mNu;
  uint32_t mNumFinite;
  double * mYp;
  svm_node ** mXp, *mTestNodes;
  std::string mRegulatedGeneName;
};

class GRNModel
{
public:
  GRNModel(const std::string& aModel,
           ExpressionMatrixProcessor& aEMP,
           uint32_t aGeneLimit = 0);
  ~GRNModel();

  template<class Container> void loadSVMTrainingData(const Container& aTrainingArrays)
  {
    uint32_t nTraining(aTrainingArrays.size());

    for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
         i != mSVMs.end();
         i++)
      (*i)->setupProblem(nTraining);

    for (typename Container::const_iterator i = aTrainingArrays.begin();
         i != aTrainingArrays.end();
         i++)
    {
      mEMP.setArray(mEMP.getIndexOfArray(*i));
      
      for (std::list<SupportVectorMachine*>::iterator j = mSVMs.begin();
           j != mSVMs.end();
           j++)
        (*j)->loadTrainingRow();
    }
  }

  void setSVMParameters(double aGamma, double aC, double aNu)
  {
    for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
         i != mSVMs.end();
         i++)
      (*i)->setParameters(aGamma, aC, aNu);
  }

  void trainSVMs()
  {
    for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
         i != mSVMs.end();
         i++)
      (*i)->train();
  }

  void trainAndSaveSVMs(const std::string& aSVMDir, bool aDontReplace = false);

  template<class Container> double testSVMs(const Container& aTestingArrays)
  {
    double testScore = 0.0;

    for (typename Container::const_iterator i = aTestingArrays.begin();
         i != aTestingArrays.end();
         i++)
    {
      mEMP.setArray(mEMP.getIndexOfArray(*i));
      
      for (std::list<SupportVectorMachine*>::iterator j = mSVMs.begin();
           j != mSVMs.end();
           j++)
      {
        double r((*j)->testOnRow());
        if (isfinite(r))
          testScore += r;
      }
    }

    return testScore;
  }

  template<class Container, class Listener>
  void testSVMs(const Container& aTestingArrays,
                Listener& aResults)
  {
    for (typename Container::const_iterator i = aTestingArrays.begin();
         i != aTestingArrays.end();
         i++)
    {
      uint32_t idx = mEMP.getIndexOfArray(*i);
      mEMP.setArray(idx);
      aResults.startRow(idx);

      double arraySum = 0.0;
      uint32_t arrayCount = 0;

      for (std::list<SupportVectorMachine*>::iterator j = mSVMs.begin();
           j != mSVMs.end();
           j++)
        aResults.result((*j)->getRegulatedGene(),
                        (*j)->testOnRow());

      aResults.endRow(idx);
    }
  }

  void saveSVMs(const std::string& aSVMDir);
  void loadSVMs(const std::string& aSVMDir);

  template<class Container> void loadArraySet(const std::string& aFilename,
                                              Container& aArrays)
  {
    std::ifstream s(aFilename.c_str());
    
    while (s.good())
    {
      std::string l;
      std::getline(s, l);
      
      if (l != "")
        aArrays.push_back(l);
    }
  }

private:
  ExpressionMatrixProcessor& mEMP;
  std::map<uint32_t, std::string> mHGNCByVertex;
  std::list<SupportVectorMachine*> mSVMs;
};
