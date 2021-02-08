#ifndef O2_ITS_TESTDATAREADER_FLOW
#define O2_ITS_TESTDATAREADER_FLOW

#include <vector>
#include <deque>
#include <memory>
#include "Rtypes.h"  // for Digitizer::Class, Double_t, ClassDef, etc
#include "TObject.h" // for TObject
#include "TGaxis.h"

#include "TFile.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTReconstruction/RawPixelReader.h"

#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <fstream>
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTReconstruction/Clusterer.h"
#include "RootInclude.h"

#include "ITSBase/GeometryTGeo.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSMFTReconstruction/DigitPixelReader.h"
#include "DataFormatsITSMFT/Digit.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTReconstruction/GBTWord.h"
#include "ITSMFTReconstruction/PayLoadCont.h"
#include "ITSMFTReconstruction/PixelData.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include <chrono>
#include <thread>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSMFTReconstruction/Clusterer.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITSMFTReconstruction/PixelData.h"

using namespace o2::framework;
using o2::itsmft::SegmentationAlpide;
namespace o2
{
namespace its
{

class TestDataReader : public Task
{
  using ChipPixelData = o2::itsmft::ChipPixelData;
  using PixelReader = o2::itsmft::PixelReader;
  using Cluster = o2::itsmft::Cluster;
  using PixelData = o2::itsmft::PixelData;

 public:
  TestDataReader() = default;
  ~TestDataReader() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  //-----------------Cluster Fuctions

  void resetColumn(int* buff) //AID cluster
  {
    std::memset(buff, -1, sizeof(int) * SegmentationAlpide::NRows);
  }
  void swapColumnBuffers()
  {
    int* tmp = mCurr;
    mCurr = mPrev;
    mPrev = tmp;
  }
  void expandPreCluster(UInt_t ip, UShort_t row, int preClusIndex)
  {
    auto& firstIndex = mPreClusterHeads[mPreClusterIndices[preClusIndex]];
    mPixels.emplace_back(firstIndex, ip);
    firstIndex = mPixels.size() - 1;
    mCurr[row] = preClusIndex;
  }

  void addNewPrecluster(UInt_t ip, UShort_t row)
  {
    mPreClusterHeads.push_back(mPixels.size());
    // new head does not point yet (-1) on other pixels, store just the entry of the pixel in the ChipPixelData
    mPixels.emplace_back(-1, ip);
    int lastIndex = mPreClusterIndices.size();
    mPreClusterIndices.push_back(lastIndex);
    mCurr[row] = lastIndex; // store index of the new precluster in the current column buffer
  }
  void adjustBoundingBox(const o2::itsmft::PixelData pix, UShort_t& rMin, UShort_t& rMax,
                         UShort_t& cMin, UShort_t& cMax) const
  {
    if (pix.getRowDirect() < rMin) {
      rMin = pix.getRowDirect();
    }
    if (pix.getRowDirect() > rMax) {
      rMax = pix.getRowDirect();
    }
    if (pix.getCol() < cMin) {
      cMin = pix.getCol();
    }
    if (pix.getCol() > cMax) {
      cMax = pix.getCol();
    }
  }
  //-----------------------------------
  void setNChips(int n)
  {
    mChips.resize(n);
  }
  std::vector<std::string> GetFName(std::string folder);

 private:
  void initChip(UInt_t first);
  void updateChip(UInt_t ip);
  void finishChip(std::vector<Cluster>* fullClus);

  std::unique_ptr<TFile> mFile = nullptr;
  //o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS> mRawReader; //AID
  std::unique_ptr<o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS>> mRawReader; //AID

  //  o2::itsmft::ChipPixelData mChipData;
  //  std::size_t rofEntry = 0, nrofdig = 0;
  //  std::unique_ptr<TFile> outFileDig;
  //  std::unique_ptr<TTree> outTreeDig; // output tree with digits
  //  std::unique_ptr<TTree> outTreeROF; // output tree with ROF records
  std::vector<ChipPixelData> mChips;
  std::vector<o2::itsmft::Digit> mDigits;
  //std::vector<o2::itsmft::QcDigit> mDigits; //qcDigit
  std::vector<o2::itsmft::Digit> mMultiDigits;

  ChipPixelData* mChipData = nullptr;
  std::string mInputName = "Split9.bin";
  int mIndexPush;
  //  int mPixelSize;
  std::vector<int> mNDigits;
  std::vector<std::string> mFolderNames;
  std::vector<std::string> mNowFolderNames;
  std::vector<std::vector<std::string>> mFileNames;
  std::vector<std::vector<std::string>> mNowFileNames;
  std::vector<std::string> mDiffFolderName;
  std::vector<std::string> mDiffFileNamePush;
  std::vector<std::vector<std::string>> mDiffFileNames;
  std::vector<std::string> NewNextFold;
  std::string mWorkDir;
  std ::string mRunType;

  int mResetCommand;
  std::string mRunID;
  int mNEvent;
  int mEventPerPush;
  int mEventRegistered;
  int mTotalPixelSize;
  static constexpr int sNError = 11;
  //			unsigned int Error[sNError];
  std::array<unsigned int, sNError> mErrors;
  std::vector<std::array<unsigned int, sNError>> mErrorsVec;
  std::vector<std::array<unsigned int, sNError>> mErrorsVecTest;
  //  int pos;
  int j;
  int mFileDone;
  int mFileID;
  int mEpNumber;
  int mRunNumber;
  int mTrackError;
  int mIndexPushEx;
  int mFileRemain;
  int mFileInfo;
  //Immediate Injection Variables//

  int mNewFileInj;
  //  int mNewFileInjAction;
  std::vector<o2::itsmft::Digit> mDigitsTest;
  std::vector<o2::itsmft::Digit> mMultiDigitsTest;
  //gsl::span<const o2::itsmft::QcDigit> DigitsArrayInp; ?? 15.02.2020
  gsl::span<o2::itsmft::Digit> DigitsArrayInp;
  int mMaxPixelSize;

  //-----------------------------Cluster stuff
  o2::itsmft::ROFRecord mROFRef; // ROF reference
  bool doCluster = kTRUE;        //kFALSE;
  std::vector<Cluster> fullClus;
  std::vector<Cluster>* fullClusPtr = nullptr;
  std::vector<Cluster> compClus;
  int mMaxBCSeparationToMask = 6000. / o2::constants::lhc::LHCBunchSpacingNS + 10;
  std::vector<ChipPixelData> mChipsOld; // previously processed chips data (for masking)
  int* mCurr;                           // pointer on the 1st row of currently processed mColumnsX
  int* mPrev;                           // pointer on the 1st row of previously processed mColumnsX
  std::vector<std::pair<int, UInt_t>> mPixels;
  std::vector<int> mPreClusterHeads; // index of precluster head in the mPixels
  std::vector<int> mPreClusterIndices;
  UShort_t mCol = 0xffff;    ///< Column being processed
  bool mNoLeftColumn = true; ///< flag that there is no column on the left to check
  int mColumn1[SegmentationAlpide::NRows + 2];
  int mColumn2[SegmentationAlpide::NRows + 2];
  std::array<PixelData, Cluster::kMaxPatternBits * 2> mPixArrBuff; //! temporary buffer for pattern calc.
  int mClustersCount = 0;                                          ///< number of clusters in the output container
  o2::itsmft::GeometryTGeo* mGeometry = nullptr;                   //! ITS OR MFT upgrade geometry

  static const std::string sRunTypeFileName;
  static const std::string sFakeRateDefConfig;
  static const std::string sThresholdDefConfig;
  static const std::string sClusterDefConfig;
};

/// create a processor spec
/// read simulated ITS digits from a root file
framework::DataProcessorSpec getTestDataReaderSpec();

} // namespace its
} // namespace o2

#endif /* O2_ITS_RAWPIXELREADER_FLOW */
