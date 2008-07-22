#include <string>
#include <inttypes.h>
#include <map>
#include <vector>
#include <list>
#include <svm.h>
#include <fstream>

class ExpressionMatrixProcessor
{
public:
  ExpressionMatrixProcessor(const std::string& aMatrixDir);
  ~ExpressionMatrixProcessor();

  uint32_t getIndexOfGene(const std::string& aGene);
  uint32_t getIndexOfArray(const std::string& aArray);
  void setArray(uint32_t aArray);
  double getDataPoint(uint32_t aGene);

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
  void save(const std::string& aFilename);
  void load(const std::string& aFilename);

  const std::string& getRegulatedGeneName()
  {
    return mRegulatedGeneName;
  }

private:
  ExpressionMatrixProcessor& mEMP;
  std::vector<uint32_t> mRegulatingGenes;
  uint32_t mRegulatedGene, mNumRegulators;
  struct svm_model* mModel;
  struct svm_problem mProblem;
  struct svm_parameter mParameter;
  double mAverage, mGamma, mC, mp;
  uint32_t mNumFinite;
  double * mYp;
  svm_node ** mXp;
  std::string mRegulatedGeneName;
};

class GRNModel
{
public:
  GRNModel(const std::string& aModel,
           ExpressionMatrixProcessor& aEMP);
  ~GRNModel();

  template<class Container> void trainSVMs(const Container& aTrainingArrays)
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

    for (std::list<SupportVectorMachine*>::iterator i = mSVMs.begin();
         i != mSVMs.end();
         i++)
      (*i)->train();
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
      
      aArrays.push_back(l);
    }
  }


private:
  ExpressionMatrixProcessor& mEMP;
  std::map<uint32_t, std::string> mHGNCByVertex;
  std::list<SupportVectorMachine*> mSVMs;
};
