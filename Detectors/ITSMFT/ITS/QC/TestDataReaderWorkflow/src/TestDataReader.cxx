#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <string>
#include "TTree.h"

#include "Framework/ControlService.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTReconstruction/GBTWord.h"
#include "ITSMFTReconstruction/PayLoadCont.h"
#include "ITSMFTReconstruction/PixelData.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSMFTReconstruction/RawPixelReader.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "ITSQCDataReaderWorkflow/TestDataReader.h"
#include "DetectorsBase/GeometryManager.h"
#include <TCanvas.h>
#include <iostream>
#include <dirent.h>
#include <cstdio>
#include <algorithm>
#include <iterator>
#include <thread>

using namespace o2::framework;
using namespace o2::itsmft;
using namespace std;
using Segmentation = o2::itsmft::SegmentationAlpide;

namespace o2
{
namespace its
{

const std::string TestDataReader::sRunTypeFileName = "Config/RunType.dat";
const std::string TestDataReader::sFakeRateDefConfig = "Config/ConfigFakeRate.dat";
const std::string TestDataReader::sThresholdDefConfig = "Config/ConfigThreshold.dat";
const std::string TestDataReader::sClusterDefConfig = "Config/ConfigCluster.dat";

void TestDataReader::init(InitContext& ic)
{
  mIndexPush = 0;

  std::ifstream RunFileType(sRunTypeFileName);
  RunFileType >> mRunType;

  LOG(DEBUG) << "OLD CONFIG: RunFileType = " << mRunType;

  if (mRunType == "FakeHitRate") {
    std::ifstream EventPush(sFakeRateDefConfig);
    EventPush >> mEventPerPush >> mTrackError >> mWorkDir;
    doCluster = kFALSE;
  }

  if (mRunType == "ThresholdScan") {
    std::ifstream EventPush(sThresholdDefConfig);
    EventPush >> mEventPerPush >> mTrackError >> mWorkDir;
    doCluster = kFALSE;
  }

  if (mRunType == "Cluster") {
    std::ifstream EventPush(sClusterDefConfig);
    EventPush >> mEventPerPush >> mTrackError >> mWorkDir;
    doCluster = kTRUE;
  }

  LOG(DEBUG) << "OLD CONFIG: EventPerPush = " << mEventPerPush << "   TrackError = " << mTrackError << "  WorkDir = " << mWorkDir;
  LOG(DEBUG) << "DONE Reset Histogram Decision";

  o2::base::GeometryManager::loadGeometry();
  mFolderNames = GetFName(mWorkDir);

  cout << "NFolder = " << mFolderNames.size() << endl;
  for (int i = 0; i < mFolderNames.size(); i++) {

    cout << "FDN = " << mFolderNames[i] << endl;

    mFileNames.push_back(GetFName(mFolderNames[i]));

    cout << "FDN File Size = " << mFileNames[i].size() << endl;

    for (int j = 0; j < mFileNames[i].size(); j++) {

      cout << "FDN File = " << mFileNames[i][j] << endl;
    }
  }
  for (int i = 0; i < sNError; i++) {
    mErrors[i] = 0;
  }

  //		mEventPerPush = 3000;
  mEventRegistered = 0;
  mTotalPixelSize = 0;

  //	GetFileName("infile");

  //

  const Int_t numOfChips = o2::itsmft::ChipMappingITS::getNChips();
  LOG(DEBUG) << "numOfChips = " << numOfChips;
  setNChips(numOfChips);
  j = 0;
  mNDigits.clear();
  mFileDone = 1;
  mFileRemain = 0;
  mNewFileInj = 1;
  mMaxPixelSize = 58700095;
  fullClusPtr = &fullClus;
  mGeometry = GeometryTGeo::Instance();
  mGeometry->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L)); // make sure T2L matrices are loaded
}

void TestDataReader::run(ProcessingContext& pc)
{
  // Keep checking new folders (and files in the folders)
  // If found, process them one by one

  //Defining all local variables
  // int j = 0;
  int NEventPre;
  int NEvent;
  double PercentDone = 0;
  int ErrorDetcted;
  Int_t ChipWasReaded = 0;
  cout << "----------------------------------------------------------" << endl
       << endl;

  cout << "New Cycle" << endl;
  if (mFolderNames.size() > 0)
    cout << "Old Folder Size = " << mFolderNames.size() << " " << mFolderNames[0] << endl;
  else
    cout << "Old Folder Size = " << mFolderNames.size() << endl;
  // Get folders in working directory, put all files in a vector
  mNowFolderNames = GetFName(mWorkDir);

  if (mNowFolderNames.size() > 0)
    cout << "Now Folder Size = " << mNowFolderNames.size() << " " << mNowFolderNames[0] << endl;
  else
    cout << "Now NFolder = " << mNowFolderNames.size() << endl;
  for (int i = 0; i < mNowFolderNames.size(); i++) {
    mNowFileNames.push_back(GetFName(mNowFolderNames[i]));
  }

  // Check for new folders comparing with the previous cycle
  std::set_difference(mNowFolderNames.begin(), mNowFolderNames.end(), mFolderNames.begin(), mFolderNames.end(), std::inserter(mDiffFolderName, mDiffFolderName.begin()));

  cout << "Difference Size Between New and Initial Runs = " << mDiffFolderName.size() << endl;

  // No new folder
  if (mDiffFolderName.size() == 0) {
    cout << "No New Run -- No Need to Reset" << endl;
    mResetCommand = 0;
    pc.outputs().snapshot(Output{"ITS", "TEST", 0, Lifetime::Timeframe}, mResetCommand);
  }

  // New folders found, send the reset signal and reload configuration
  if (mDiffFolderName.size() > 0) {
    cout << "New Run Started -- Reset All Histograms" << endl;
    mResetCommand = 1;
    pc.outputs().snapshot(Output{"ITS", "TEST", 0, Lifetime::Timeframe}, mResetCommand);
    for (int i = 0; i < sNError; i++) {
      mErrors[i] = 0;
    }
    mErrorsVec.clear();
    mResetCommand = 0;
    mNewFileInj = 1;
    mIndexPush = 0;
    std::ifstream RunFileType(sRunTypeFileName);
    RunFileType >> mRunType;

    LOG(INFO) << "NEW CONFIG: RunFileType = " << mRunType;

    if (mRunType == "FakeHitRate") {
      std::ifstream EventPush(sFakeRateDefConfig);
      EventPush >> mEventPerPush >> mTrackError >> mWorkDir;
    }

    if (mRunType == "ThresholdScan") {
      std::ifstream EventPush(sThresholdDefConfig);
      EventPush >> mEventPerPush >> mTrackError >> mWorkDir;
    }

    LOG(INFO) << "NEW CONFIG: EventPerPush = " << mEventPerPush << "   TrackError = " << mTrackError << "  mWorkDir = " << mWorkDir;
    LOG(INFO) << "DONE Reset Histogram Decision";
  }

  LOG(DEBUG) << "Start Creating New Now Vector";

  LOG(INFO) << "Get IN LOOP";
  // Checking for new files, extract them into new vector
  for (int i = 0; i < mFolderNames.size(); i++) {
    std::set_difference(mNowFileNames[i].begin(), mNowFileNames[i].end(), mFileNames[i].begin(), mFileNames[i].end(), std::inserter(mDiffFileNamePush, mDiffFileNamePush.begin()));
    mDiffFileNames.push_back(mDiffFileNamePush);
    cout << "Difference File Size Between New and Initial Runs " << mDiffFileNames[i].size() << endl;

    mDiffFileNamePush.clear();
  }

  LOG(INFO) << "DONE GRABING Existing";

  //Getting the new files from new folders that does not exist in the previous cycle
  for (int i = mFolderNames.size(); i < mNowFolderNames.size(); i++) {
    mDiffFileNames.push_back(mNowFileNames[i]);
    cout << "New File Size Between New and Initial Runs " << mDiffFileNames[i].size() << endl;
  }

  LOG(INFO) << "Total New Files = " << mDiffFileNames.size();

  LOG(DEBUG) << "DONE Creating Difference";

  LOG(INFO) << "mDiffFileNames Size = " << mDiffFileNames.size();

  LOG(INFO) << "Start Loop";

  //Start Decoding New Files by loop through the new file vector

  for (int i = 0; i < mNowFolderNames.size(); i++) {

    LOG(INFO) << "i = " << i << "    mDiffFileNames[i].size() = " << mDiffFileNames[i].size();

    //Getting the folder name ID

    int pos = mNowFolderNames[i].find_last_of("/");

    if (pos != string::npos)
      mRunID = mNowFolderNames[i].substr(pos + 1);

    LOG(INFO) << "FileDone = " << mFileDone << endl;

    //Reading files one by one

    // if (mDiffFileNames[i].size() > 0 && mFileDone == 1)  //AID
    if (mDiffFileNames[i].size() > 0) {

      mFileRemain = mDiffFileNames[i].size();
      //mFileDone = 0;
      cout << "RunID = " << mRunID << endl;
      cout << "File Location = " << mDiffFileNames[i][0] << endl;
      cout << i << " Remaining = " << mFileRemain << endl;
      //Getting the RunID
      size_t last_index1 = mRunID.find_last_not_of("0123456789");
      string RunIDS = mRunID.substr(last_index1 + 1);
      //Getting the FileID
      string FileIDS;

      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      pos = mDiffFileNames[i][0].find_last_of("/");
      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 1.5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      if (pos != string::npos)
        FileIDS = mDiffFileNames[i][0].substr(pos + 1);
      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 1.666 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      cout << mDiffFileNames[i][0].substr(pos + 8, 1) << endl;
      mEpNumber = stoi(mDiffFileNames[i][0].substr(pos + 8, 1));
      cout << "ep: " << mEpNumber << endl;
      cout << "Before FileIDS = " << FileIDS << endl;

      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      size_t last_index2 = FileIDS.find_last_not_of("0123456789");
      FileIDS = FileIDS.substr(last_index2 + 1);
      //Extracting the RunID and File to integer in order to make it possible to inject to the QC
      mRunNumber = std::stoi(RunIDS);
      mFileID = std::stoi(FileIDS);

      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      ofstream fout(Form("ErrorData/ErrorLogRun%d_File%d.dat", mRunNumber, mFileID));
      fout << " START OF ERROR REPORT For Run " << mRunNumber << "  File " << mFileID << endl;
      //Reading the first file (Because it is one by one)
      mInputName = mDiffFileNames[i][0];

      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      //mEventRegistered = 0;
      //LOG(INFO) << "mRunNumber = " << mRunNumber;

      //Inject fake thing digits for the QC to update immediately
      //mDigitsTest is the fake digit for updating the QC immediately on the QC GUI

      cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      cout << "[AID] 6" << endl;

      // if (mNewFileInj == 1) {
      if (mFileDone == 1) { //AID
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Debug Point 6 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cout << "New File Injected, Now Updating the Canvas and Light" << endl;
        mDigitsTest.emplace_back(0, 0, 0, 0);
        // mMultiDigitsTest.push_back(mDigitsTest[0]);
        mErrorsVecTest.push_back(mErrors);
        mFileDone = 2;
        mFileInfo = mFileDone + mFileRemain * 10;
        pc.outputs().snapshot(Output{"ITS", "Run", 0, Lifetime::Timeframe}, mRunNumber);
        pc.outputs().snapshot(Output{"ITS", "EP", 0, Lifetime::Timeframe}, mEpNumber);
        pc.outputs().snapshot(Output{"ITS", "File", 0, Lifetime::Timeframe}, mFileID);
        pc.outputs().snapshot(Output{"ITS", "Error", 0, Lifetime::Timeframe}, mErrorsVecTest[0]);
        cout << " @@@@@@@@@@@@@@@@@@@@@@@@@22 I'm sending mFileInfo " << mFileInfo << endl;
        pc.outputs().snapshot(Output{"ITS", "Finish", 0, Lifetime::Timeframe}, mFileInfo);
        pc.outputs().snapshot(Output{"ITS", "DIGITS", 0, Lifetime::Timeframe}, mMultiDigitsTest);
        //pc.outputs().snapshot(Output{ "ITS", "CLUSTERS", 0, Lifetime::Timeframe }, fullClus);

        mNewFileInj = 0;
        mEventRegistered = 0;
        mErrorsVecTest.clear();
        mDigitsTest.clear();
        mMultiDigitsTest.clear();
        //-----------------AID Improtant
        LOG(INFO) << "[AID] Opened new FILE!!!! inpName = " << mInputName;
        mRawReader = std::make_unique<o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS>>(); //AID how to solve problem with new readers ??
        mRawReader->setPadding128(true);                                                         // payload GBT words are padded to 16B
        mRawReader->setVerbosity(0);
        mRawReader->setMinTriggersToCache(1000000);
        mRawReader->openInput(mInputName);
        //---------------------------------

        if (mFolderNames.size() < mNowFolderNames.size())
          mFileNames.push_back(NewNextFold);
        cout << "Done!!! You should see the Canvas Updated " << endl;
        break;
      }

      //DONE Injection//
      /*
 // AID Commented
      o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS> mRawReader;
      mRawReader.setPadding128(true); // payload GBT words are padded to 16B
      mRawReader.setVerbosity(0);
      mRawReader.setMinTriggersToCache(1025);

      mRawReader.openInput(mInputName);
*/
      int Index = 0;
      int IndexMax = -1;
      int NChip = 0;
      int NChipMax = -1;
      int TimePrint = 0;
      using RawReader = o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS>;
      auto& rawErrorReader = reinterpret_cast<RawReader&>(*mRawReader);
      //----------------Clusters
      mROFRef.clear();
      auto& currROFIR = mROFRef.getBCData();
      auto& currROFEntry = mROFRef.getEntry();
      Int_t mClustersCount = 0; //fullClus->size();
                                //------------------------------------
      int curlinkid = -1;

      mFileDone = 1; //AID

      while ((mChipData = mRawReader->getNextChipData(mChips))) {
        auto chipID = mChipData->getChipID();
        if (NChip < NChipMax)
          break;
        //---------------------------------------------------------------------------
        if (doCluster) {
          //------------------------getting ROF:
          if (!(mChipData->getInteractionRecord() == currROFIR)) {              // new ROF starts
            mROFRef.setNEntries(mClustersCount - currROFEntry.getFirstEntry()); // number of entries in this ROF
            if (!currROFIR.isDummy()) {
              LOG(INFO) << "[AID] ITS: clusterizing new ROFrame, Orbit :" << mChipData->getInteractionRecord().orbit << " BC: " << mChipData->getInteractionRecord().bc;

              currROFEntry.setFirstEntry(mClustersCount);
              currROFIR = mChipData->getInteractionRecord();
              mROFRef.setROFrame(mChipData->getROFrame()); // TODO: outphase this
            }
          }
          //--------------------mask pixels
          if (mMaxBCSeparationToMask > 0) { // mask pixels fired from the previous ROF
            if (mChipsOld.size() < mChips.size()) {
              mChipsOld.resize(mChips.size()); // expand buffer of previous ROF data
            }
            const auto& chipInPrevROF = mChipsOld[chipID];
            if (std::abs(currROFIR.differenceInBC(chipInPrevROF.getInteractionRecord())) < mMaxBCSeparationToMask) {
              mChipData->maskFiredInSample(mChipsOld[chipID]);
            }
          }
          //-----------------------------------
          auto validPixID = mChipData->getFirstUnmasked();
          if (validPixID < mChipData->getData().size()) { // chip data may have all of its pixels masked!
            initChip(validPixID++);
            for (; validPixID < mChipData->getData().size(); validPixID++) {
              if (!mChipData->getData()[validPixID].isMasked()) {
                updateChip(validPixID);
              }
            }

            finishChip(fullClusPtr); //we are not going to use MC true
          }
          if (mMaxBCSeparationToMask > 0) { // current chip data will be used in the next ROF to mask overflow pixels
            mChipsOld[chipID].swap(*mChipData);
          }
        }
        //---------------------------------------------------------------------------------

        // cout << "linkid0 status " << mRawReader->mRUDecodeVec[0].links[0]  << endl;
        // cout << "linkid1 status " << mRawReader->mRUDecodeVec[0].links[1]  << endl;
        // cout << "linkid2 status " << mRawReader->mRUDecodeVec[0].links[2]  << endl;

        //	cout << "Pass Chip" << endl;
        //
        //	[AID]  step1: we are getting all readers: pixels and errors arrays
        //
        // LOG(INFO) <<"LINK0 " << mRawReader->mRUDecodeVec[0].links[0] << " LINK1 " << mRawReader->mRUDecodeVec[0].links[1] << " LINK2 " << mRawReader->mRUDecodeVec[0].links[2];
        if (mRawReader->mRUDecodeVec[0].links[0] == 0) {
          curlinkid = 0;
        } else if (mRawReader->mRUDecodeVec[0].links[1] == 0) {

          curlinkid = 1;
        } else if (mRawReader->mRUDecodeVec[0].links[2] == 0) {
          curlinkid = 2;
        }
        if (curlinkid == -1) {
          cout << "Incorrect link id" << endl;
          continue;
        }
        //        const auto* ruInfo = mRawReader.getCurrRUDecodeData()->ruInfo;
        const auto* ruInfo = rawErrorReader.getCurrRUDecodeData()->ruInfo;

        const auto& statRU = rawErrorReader.getRUDecodingStatSW(ruInfo->idSW, curlinkid);
        // cout << "linkid isss" << curlinkid << endl;

        const auto& pixels = mChipData->getData();
        ChipWasReaded = 1;
        int pixelSize = mChipData->getData().size();
        //cout << "pixel size isss" << pixelSize << endl;
        NEvent = statRU->nPackets;
        //  cout << "NEvent is " << NEvent <<endl;
        //  mTotalPixelSize = mTotalPixelSize + pixelSize; [AID comment]
        /*
 //AID comment
 //
 //Here we are geeting size of samples, compare it with config file and create ANOTHER array with only description of samples. WHY?!
        if (NEvent > (mEventRegistered + 1) * mEventPerPush) {

          if (mTotalPixelSize < mMaxPixelSize || mTotalPixelSize == mMaxPixelSize) {
            cout << "Digit OK for 1 Push" << endl;
            mNDigits.push_back(mTotalPixelSize);
            mEventRegistered = mEventRegistered + 1;
            mErrorsVec.push_back(mErrors);
            cout << "TotalPixelSize = " << mTotalPixelSize << "  Pushed" << endl;
          }

          if (mTotalPixelSize > mMaxPixelSize) {
            cout << "Digit Split into 2 Pushes" << endl;
            mNDigits.push_back(mTotalPixelSize / 2);
            mNDigits.push_back(mTotalPixelSize / 2);
            mErrorsVec.push_back(mErrors);
            mErrorsVec.push_back(mErrors);
            mEventRegistered = mEventRegistered + 1;
            cout << "TotalPixelSize = " << mTotalPixelSize << "  Pushed" << endl;
          }

          mTotalPixelSize = 0;
        }
*/

        //G(INFO) << "[AID] 8";
        //
        //Printing the Event Number Processed

        //[AID] step2: PrintingOUT

        if (NEvent % 1000000 == 0 && TimePrint == 0) {
          cout << "Event Number = " << NEvent << endl;
          TimePrint = 1;
        }

        if (NEvent % 100000 != 0)
          TimePrint = 0;

        for (int i = 0; i < o2::itsmft::GBTLinkDecodingStat::NErrorsDefined; i++) {
          if ((mErrors[i] + (unsigned int)statRU->errorCounts[i]) < 4294967295)
            mErrors[i] = mErrors[i] + (unsigned int)statRU->errorCounts[i];
        }

        //[AID] step 3: Proceeding of Errors
        if (mTrackError == 1) {
          if (NEventPre != NEvent) {
            ErrorDetcted = 0;

            for (int i = 0; i < o2::itsmft::GBTLinkDecodingStat::NErrorsDefined; i++) {
              if ((int)statRU->errorCounts[i] > 0)
                ErrorDetcted = 1;
            }

            if (ErrorDetcted == 1) {
              fout << "Event Number = " << NEvent << endl;
              fout << " ------------------------------------------------" << endl;
              for (int i = 0; i < o2::itsmft::GBTLinkDecodingStat::NErrorsDefined; i++) {
                if (statRU->errorCounts[i]) {
                  fout << "Error ID " << i << " " << o2::itsmft::GBTLinkDecodingStat::ErrNames[i] << endl;
                }
              }
              fout << " ------------------------------------------------" << endl;
            }
          }
        }
        //[AID] step 4: time to fill vector of DIGITS from pixelData:
        // auto ChipID = mChipData->getChipID();
        //  cout << "ChiID is " << ChipID << endl;
        //  if(ChipID == 0) cout << "ChipID is :" << ChipID << endl;
        for (auto& pixel : pixels) {
          if (Index < IndexMax)
            break;
          int col = pixel.getCol();
          int row = pixel.getRow();
          mDigits.emplace_back((UShort_t)chipID, row, col, 0);
          //	    mDigits.emplace_back( QcDigit(ChipID, row, col));  //qcDigits

          //  cout << "ChipID is in:" << (UShort_t)ChipID << " " << NEvent << " " << row << " " << col << endl;
          Index = Index + 1;
        }

        NChip = NChip + 1;
        NEventPre = NEvent;
        //[AID] step 5. (NEW) time to push to QC!
        if (NEvent > (mEventRegistered + 1) * mEventPerPush) {
          mFileDone = 0;
          mFileInfo = mFileDone + mFileRemain * 10;
          mEventRegistered = mEventRegistered + 1;

          //if (NEvent > 40000) break;
          pc.outputs().snapshot(Output{"ITS", "Run", 0, Lifetime::Timeframe}, mRunNumber);
          pc.outputs().snapshot(Output{"ITS", "File", 0, Lifetime::Timeframe}, mFileID);
          pc.outputs().snapshot(Output{"ITS", "Error", 0, Lifetime::Timeframe}, mErrors);
          //pc.outputs().snapshot(Output{"ITS", "EP", 0, Lifetime::Timeframe}, mEpNumber); //What its?
          pc.outputs().snapshot(Output{"ITS", "EP", 0, Lifetime::Timeframe}, mEpNumber); //qcDigit
                                                                                         //if (doCluster)
                                                                                         //pc.outputs().snapshot(Output{ "ITS", "CLUSTERS", 0, Lifetime::Timeframe }, fullClus);
          LOG(INFO) << "Sampled mDigits.size() " << mDigits.size();
          LOG(INFO) << "FileDone = " << mFileDone;
          LOG(INFO) << "FileRemain = " << mFileRemain;

          pc.outputs().snapshot(Output{"ITS", "Finish", 0, Lifetime::Timeframe}, mFileInfo); //?????
          //DigitsArrayInp = gsl::span<const o2::itsmft::QcDigit>(mDigits->data(), (*mDigits).size());
          //	DigitsArrayInp = gsl::make_span(&(*mDigits)[0], (*mDigits).size());
          //	pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, DigitsArrayInp);
          pc.outputs().snapshot(Output{"ITS", "DIGITS", 0, Lifetime::Timeframe}, mDigits);
          mDigits.clear();
          //mFileDone=0;
          LOG(INFO) << " NEvent= " << NEvent << " mEventPerPush " << mEventPerPush << " " << mEventRegistered << " Cluster Vector: " << fullClus.size();
          break; //AID really important
        }
      }

      LOG(INFO) << "Run " << mNowFolderNames[i] << " File " << mInputName << "    Integrated Raw Pixel Pushed " << mDigits.size();

      if (mFileDone == 1) {
        LOG(INFO) << "END OF THE FILE" << endl;
        mFileInfo = mFileDone + mFileRemain * 10;
        if (mFolderNames.size() < mNowFolderNames.size())
          mFileNames.push_back(NewNextFold);
        mFileNames[i].push_back(mInputName);

        pc.outputs().snapshot(Output{"ITS", "Run", 0, Lifetime::Timeframe}, mRunNumber);
        pc.outputs().snapshot(Output{"ITS", "File", 0, Lifetime::Timeframe}, mFileID);
        pc.outputs().snapshot(Output{"ITS", "Error", 0, Lifetime::Timeframe}, mErrors);
        pc.outputs().snapshot(Output{"ITS", "EP", 0, Lifetime::Timeframe}, mEpNumber);     //What its?
        pc.outputs().snapshot(Output{"ITS", "Finish", 0, Lifetime::Timeframe}, mFileInfo); //?????
                                                                                           // if (doCluster)
        pc.outputs().snapshot(Output{"ITS", "CLUSTERS", 0, Lifetime::Timeframe}, fullClus);

        //DigitsArrayInp = gsl::span<const o2::itsmft::QcDigit>(mDigits->data(), (*mDigits).size());
        // DigitsArrayInp = gsl::make_span(&(*mDigits)[0], (*mDigits).size());
        //        pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, DigitsArrayInp);
        pc.outputs().snapshot(Output{"ITS", "DIGITS", 0, Lifetime::Timeframe}, mDigits);

        for (int i = 0; i < o2::itsmft::GBTLinkDecodingStat::NErrorsDefined; i++)
          mErrors[i] = 0;
        mDigits.clear();
      }
    } //[AID] What do we have those loops over directories?
  }

  mFolderNames.clear();
  NewNextFold.clear();
  mFolderNames = mNowFolderNames;

  mNowFolderNames.clear();
  mNowFileNames.clear();
  mDiffFileNames.clear();
  mDiffFolderName.clear();

  LOG(DEBUG) << "Pushing Reset Histogram Decision";

  cout << "Start Sleeping" << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << " " << endl;

  //std::this_thread::sleep_for(std::chrono::milliseconds(3000));

  if (ChipWasReaded == 0)
    std::this_thread::sleep_for(std::chrono::milliseconds(5000)); //AID WHY ITS HERE!??
}

std::vector<string> TestDataReader::GetFName(std::string folder)
{

  DIR* dirp;

  char cstr[folder.size() + 1];
  strcpy(cstr, folder.c_str());
  dirp = opendir(cstr);
  std::vector<string> names;
  //string search_path = folder + "/*";
  if (dirp) {
    struct dirent* directory;
    while ((directory = readdir(dirp)) != nullptr) {

      //printf("%s\n", directory->d_name);

      if (!(!strcmp(directory->d_name, ".") || !strcmp(directory->d_name, "..")))
        names.push_back(folder + "/" + directory->d_name);
    }

    closedir(dirp);
  }

  cout << "names size = " << names.size() << endl;
  return (names);
}

DataProcessorSpec getTestDataReaderSpec()
{
  return DataProcessorSpec{
    "Raw-Pixel-Reader",
    Inputs{},
    Outputs{
      OutputSpec{"ITS", "DIGITS", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "CLUSTERS", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "TEST", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "Error", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "Run", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "File", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "Finish", 0, Lifetime::Timeframe},
      OutputSpec{"ITS", "EP", 0, Lifetime::Timeframe},
    },
    AlgorithmSpec{adaptFromTask<TestDataReader>()},
  };
}

void TestDataReader::updateChip(UInt_t ip)
{
  const auto pix = mChipData->getData()[ip];
  UShort_t row = pix.getRowDirect(); // can use getRowDirect since the pixel is not masked
  if (mCol != pix.getCol()) {        // switch the buffers
    swapColumnBuffers();
    resetColumn(mCurr);
    mNoLeftColumn = false;
    if (pix.getCol() > mCol + 1) {
      mCol = pix.getCol();
      addNewPrecluster(ip, row);
      mNoLeftColumn = true;
      return;
    }
    mCol = pix.getCol();
  }

  Bool_t orphan = true;

  if (mNoLeftColumn) { // check only the row above
    if (mCurr[row - 1] >= 0) {
      expandPreCluster(ip, row, mCurr[row - 1]); // attach to the precluster of the previous row
      return;
    }
  } else {
    int neighbours[]{mCurr[row - 1], mPrev[row], mPrev[row + 1], mPrev[row - 1]};
    for (auto pci : neighbours) {
      if (pci < 0) {
        continue;
      }
      if (orphan) {
        expandPreCluster(ip, row, pci); // attach to the adjascent precluster
        orphan = false;
        continue;
      }
      // reassign precluster index to smallest one
      if (mPreClusterIndices[pci] < mPreClusterIndices[mCurr[row]]) {
        mPreClusterIndices[mCurr[row]] = mPreClusterIndices[pci];
      } else {
        mPreClusterIndices[pci] = mPreClusterIndices[mCurr[row]];
      }
    }
  }
  if (orphan) {
    addNewPrecluster(ip, row); // start new precluster
  }
}

void TestDataReader::initChip(UInt_t first)
{
  // init chip with the 1st unmasked pixel (entry "from" in the mChipData)
  mPrev = mColumn1 + 1;
  mCurr = mColumn2 + 1;
  resetColumn(mCurr);

  mPixels.clear();
  mPreClusterHeads.clear();
  mPreClusterIndices.clear();
  auto pix = mChipData->getData()[first];
  mCol = pix.getCol();

  //addNewPrecluster(first, pix.getRowDirect()); // save on .size() calls ?
  mCurr[pix.getRowDirect()] = 0; // can use getRowDirect since the pixel is not masked
  // start the first pre-cluster
  mPreClusterHeads.push_back(0);
  mPreClusterIndices.push_back(0);
  mPixels.emplace_back(-1, first); // id of current pixel
  mNoLeftColumn = true;            // flag that there is no column on the left to check yet
}

void TestDataReader::finishChip(std::vector<Cluster>* fullClus)
{

  constexpr Float_t SigmaX2 = Segmentation::PitchRow * Segmentation::PitchRow / 12.; // FIXME
  constexpr Float_t SigmaY2 = Segmentation::PitchCol * Segmentation::PitchCol / 12.; // FIXME

  const auto& pixData = mChipData->getData();
  //cout<<"Finish Chip 1"<<endl;
  for (int i1 = 0; i1 < mPreClusterHeads.size(); ++i1) {
    const auto ci = mPreClusterIndices[i1];
    if (ci < 0) {
      continue;
    }
    UShort_t rowMax = 0, rowMin = 65535;
    UShort_t colMax = 0, colMin = 65535;
    int nlab = 0, npix = 0;
    int next = mPreClusterHeads[i1];
    //cout<<"Finish Chip 2"<<endl;
    while (next >= 0) {
      const auto& pixEntry = mPixels[next];
      const auto pix = pixData[pixEntry.second];
      if (npix < mPixArrBuff.size()) {
        mPixArrBuff[npix++] = pix; // needed for cluster topology
        adjustBoundingBox(pix, rowMin, rowMax, colMin, colMax);
        next = pixEntry.first;
      } else {
        cout << "[ERROR] Cluster size " << npix + 1 << " exceeds the buffer size";
      }
    }
    // cout<<"Finish Chip 3"<<endl;

    mPreClusterIndices[i1] = -1;
    for (int i2 = i1 + 1; i2 < mPreClusterHeads.size(); ++i2) {
      if (mPreClusterIndices[i2] != ci) {
        continue;
      }
      next = mPreClusterHeads[i2];
      while (next >= 0) {
        const auto& pixEntry = mPixels[next];
        const auto pix = pixData[pixEntry.second]; // PixelData
        if (npix < mPixArrBuff.size()) {
          mPixArrBuff[npix++] = pix; // needed for cluster topology
          adjustBoundingBox(pix, rowMin, rowMax, colMin, colMax);
          next = pixEntry.first;
        } else {
          cout << "[ERROR] Cluster size " << npix + 1 << " exceeds the buffer size";
        }
      }
      mPreClusterIndices[i2] = -1;
    }
    //cout<<"Finish Chip 4"<<endl;
    UShort_t rowSpan = rowMax - rowMin + 1, colSpan = colMax - colMin + 1;
    Cluster clus;
    //clus.setROFrame(mChipData->getROFrame()); //AID removed in new O2 version
    clus.setSensorID(mChipData->getChipID());
    clus.setNxNzN(rowSpan, colSpan, npix);
#ifdef _ClusterTopology_
    UShort_t colSpanW = colSpan, rowSpanW = rowSpan;
    if (colSpan * rowSpan > Cluster::kMaxPatternBits) { // need to store partial info

      if (colSpan > rowSpan) {
        if ((colSpanW = Cluster::kMaxPatternBits / rowSpan) == 0) {
          colSpanW = 1;
          rowSpanW = Cluster::kMaxPatternBits;
        }
      } else {
        if ((rowSpanW = Cluster::kMaxPatternBits / colSpan) == 0) {
          rowSpanW = 1;
          colSpanW = Cluster::kMaxPatternBits;
        }
      }
    }
    clus.setPatternRowSpan(rowSpanW, rowSpanW < rowSpan);
    clus.setPatternColSpan(colSpanW, colSpanW < colSpan);
    clus.setPatternRowMin(rowMin);
    clus.setPatternColMin(colMin);
    for (int i = 0; i < npix; i++) {
      const auto pix = mPixArrBuff[i];
      unsigned short ir = pix.getRowDirect() - rowMin, ic = pix.getCol() - colMin;
      if (ir < rowSpanW && ic < colSpanW) {
        clus.setPixel(ir, ic);
      }
    }
#endif //_ClusterTopology_

    //cout<<"Finish Chip 5"<<endl;
    if (fullClus) { // do we need conventional clusters with full topology and coordinates?
                    //cout<<"Finish Chip 6"<<endl;
      fullClus->push_back(clus);
      //cout<<"After Push"<<endl;
      //    if (clus.getNPix()>1) cout<<"[AID] added new cluster! with size: " << clus.getNPix()<<endl;
      Cluster& c = fullClus->back();
      Float_t x = 0., z = 0.;
      for (int i = npix; i--;) {
        x += mPixArrBuff[i].getRowDirect();
        z += mPixArrBuff[i].getCol();
      }
      Point3D<float> xyzLoc;
      Segmentation::detectorToLocalUnchecked(x / npix, z / npix, xyzLoc);

      auto xyzTra = mGeometry->getMatrixT2L(mChipData->getChipID()) ^ (xyzLoc); // inverse transform from Local to Tracking frame
      c.setPos(xyzTra);
      c.setErrors(SigmaX2, SigmaY2, 0.f);
    }
    mClustersCount++;
  }
}

} // namespace its
} // namespace o2
