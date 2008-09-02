#include <inttypes.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <list>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class LevelBinner
{
public:
  LevelBinner()
  {
    uint32_t i = 0;
    double f = kLowestBinTop;

    for (i = 0; i < kBinCount - 1; i++, f *= kBinFactor)
      mLevels.push_back(f);
  }

  uint32_t
  getBinByValue(double v)
  {
    return (std::lower_bound(mLevels.begin(), mLevels.end(), v) -
            mLevels.begin());
  }

  std::string
  getBinName(uint32_t bin)
  {
    char buf[40];
    if (bin == 0)
    {
      snprintf(buf, 40, "<%g", kLowestBinTop);
      return buf;
    }
    else if (bin >= (kBinCount - 1))
    {
      snprintf(buf, 40, ">%g", kLowestBinTop * pow(kBinFactor, kBinCount - 2));
      return buf;
    }

    double lower = kLowestBinTop, upper = kLowestBinTop * kBinFactor;

    while (--bin)
    {
      lower *= kBinFactor;
      upper *= kBinFactor;
    }

    snprintf(buf, 40, "%g-%g", lower, upper);
    return buf;
  }

  static const double kLowestBinTop = 0.25;
  static const double kBinFactor = 2;
  static const uint32_t kBinCount = 18;

private:
  std::vector<double> mLevels;
};

LevelBinner gLevelBinner;

class GeneProfileBuilder
{
public:
  GeneProfileBuilder(const std::string& aMatrixDir)
    : mBinnedData(NULL)
  {
    std::list<std::string> arrays;
    fs::path matrixpath(aMatrixDir);
    fs::path arrayfile(matrixpath);
    arrayfile /= "arrays";
    readList(arrays, arrayfile);
    fs::path genefile(matrixpath);
    genefile /= "genes";
    readList(mGenes, genefile);
    
    uint32_t nGenes = mGenes.size();
    double* row = new double[nGenes];
    
    matrixpath /= "data";
    FILE* matrix = fopen(matrixpath.string().c_str(), "r");
    
    mBinnedData = new uint32_t[mGenes.size() * LevelBinner::kBinCount];
    memset(mBinnedData, 0, sizeof(uint32_t) * mGenes.size() *
           LevelBinner::kBinCount);

    for (std::list<std::string>::iterator i(arrays.begin());
         i != arrays.end(); i++)
    {
      fread(row, nGenes * sizeof(double), 1, matrix);
      processRow(row);
    }

    delete [] row;
  }

  void
  writeProfile(std::ostream& aDest)
  {
    aDest << "\"Gene\"";
    for (uint32_t i = 0; i < LevelBinner::kBinCount; i++)
      aDest << ",\"" << gLevelBinner.getBinName(i) << "\"";

    aDest << std::endl;

    uint32_t* p = mBinnedData;
    for (std::list<std::string>::iterator i(mGenes.begin());
         i != mGenes.end(); i++)
    {
      aDest << "\"" << *i << "\"";
      for (uint32_t j = 0; j < LevelBinner::kBinCount; j++)
        aDest << "," << *p++;
      aDest << std::endl;
    }
  }

  ~GeneProfileBuilder()
  {
    if (mBinnedData)
      delete [] mBinnedData;
  }

private:
  template<typename Container>
  void
  readList(Container& aStorage, const fs::path& aPath)
  {
    std::ifstream s(aPath.string().c_str());
    while (s.good())
    {
      std::string l;
      std::getline(s, l);
      aStorage.push_back(l);
    }
  }

  void
  processRow(double* aData)
  {
    uint32_t* p = mBinnedData;
    for (uint32_t i = 0, l = mGenes.size(); i < l; i++)
    {
      if (finite(aData[i]))
        p[gLevelBinner.getBinByValue(aData[i])]++;
      p += LevelBinner::kBinCount;
    }
  }

  uint32_t* mBinnedData;
  std::list<std::string> mGenes;
};

int
main(int argc, char** argv)
{
  po::options_description desc;
  std::string matrixdir, outputfile;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "The directory set up by SOFT2Matrix")
    ("output", po::value<std::string>(&outputfile),
     "The file to store gene expression levels into")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("matrixdir"))
      wrong = "matrixdir";
    if (!vm.count("output"))
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

  GeneProfileBuilder gpb(matrixdir);
  std::ofstream s(outputfile.c_str());

  gpb.writeProfile(s);
}
