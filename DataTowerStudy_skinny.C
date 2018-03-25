// 'DataTowerStudy_skinny.C'
// Derek Anderson
// 09.10.2017
//
// Use this to extract the tower distributions
// from the data 'femtoDst' tree.


#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t nDist(2);
static const UInt_t nEffPar(5);
static const UInt_t nNumBins(4);
static const UInt_t nSpecies(2);
//static const UInt_t nHotTwrs(302);
static const UInt_t nHotTwrs(41);
static const UInt_t nTrkMax(1000);
static const UInt_t nTwrMax(10000);
static const UInt_t nMatMax(100);
static const UInt_t fTrgMode(1);
// io parameters
static const TString sTree("Gfmtodst");
static const TString sInDefault("/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.epCalc2.root");
static const TString sOutDefault("pp200r9.et9vz55nHot41pi0.root");
static const TString sTrkDirs[nDist + 1] = {"AllTracks", "RecoilTracks", "Track2D"};
static const TString sTwrDirs[nDist + 1] = {"AllTowers", "RecoilTowers", "Tower2D"};
static const TString sEpDir("TowerEp");



void DataTowerStudy_skinny(const Bool_t isInBatchMode=false, const TString sInput=sInDefault, const TString sOutput=sOutDefault) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning tower script..." << endl;

  // constants
  const Bool_t   isInTrgMode0  = (fTrgMode == 0);
  const Bool_t   isInTrgMode1  = (fTrgMode == 1);
  const Bool_t   isInTrgMode2  = (fTrgMode == 2);
  const Bool_t   isInTrgMode3  = (fTrgMode == 3);
  const Bool_t   removeHotTwrs = true;
  const Double_t fRecoilMax    = TMath::PiOver2();
  const Double_t mPion         = 0.140;

  // event / trigger constants
  const Int_t    adcMax   = 6004;
  const Double_t zVtxMax  = 55.;
  const Double_t rVtxMax  = 2.;
  const Double_t eHmin    = 0.5;
  const Double_t eFmin    = 0.5;
  const Double_t pProjMax = 3.;
  const Double_t hTrgMax  = 0.9;
  const Double_t eTtrgMin = 9.;
  const Double_t eTtrgMax = 20.;
  const Double_t tspMin   = 0.;
  const Double_t tspMax   = 0.08;

  // track constants
  const Double_t nFitMin  = 15;
  const Double_t rFitMin  = 0.52;
  const Double_t dcaMax   = 1.0;
  const Double_t hTrkMax  = 1.0;
  const Double_t pTtrkMin = 0.2;

  // tower constants
  const Double_t hTwrMax  = 0.9;
  const Double_t eTwrMin  = 0.2;
  const Double_t eCorrMin = 0.2;
  const Double_t eTwrHard = 10.;
  const Double_t pidCut   = 3.;


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fOutput || !fInput) {
    cerr << "PANIC: couldn't open file!" << endl;
    return;
  }

  // get TTree
  TTree *tInput;
  fInput -> GetObject(sTree.Data(), tInput);
  if (!tInput) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  cout << "    Files opened, trees grabbed." << endl;


  // declre leaf types
  UInt_t   fUniqueID;
  UInt_t   fBits;
  Long64_t runNumber;
  Long64_t eventNumber;
  Int_t    trigID;
  Int_t    nGlobalTracks;
  Int_t    nPrimaryTracks;
  Int_t    refMult;
  Double_t vpdVz;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;
  Double_t bbcZVertex;
  Double_t zdcCoincidenceRate;
  Double_t bbcCoincidenceRate;
  Double_t backgroundRate;
  Double_t bbcBlueBackgroundRate;
  Double_t bbcYellowBackgroundRate;
  Double_t refMultPos;
  Double_t refMultNeg;
  Double_t bTOFTrayMultiplicity;
  Int_t    nVerticies;
  Double_t MagF;
  Double_t VrtxRank;
  Int_t    FlagEvent_TrgTrkMisMtch;
  Float_t  Etsp;
  Int_t    ETwrdidT;
  Int_t    ETwradc11;
  Float_t  ETwreneT0;
  Float_t  ETwreT;
  Float_t  ETwrENET0;
  Float_t  ETwrphT;
  Float_t  ETwrPTower;
  Float_t  ETwrpidTower;
  Int_t    ETwrmoduleT;
  Float_t  EClustEneT0;
  Float_t  EClustetav1;
  Float_t  EClustphiv1;
  Float_t  EEstrpen01;
  Float_t  EEstrpen02;
  Float_t  EEstrpen03;
  Float_t  EEstrpen0;
  Float_t  EEstrpen1;
  Float_t  EEstrpen2;
  Float_t  EEstrpen3;
  Float_t  EEstrpen4;
  Float_t  EEstrpen5;
  Float_t  EEstrpen6;
  Float_t  EEstrpen7;
  Float_t  EEstrpen8;
  Float_t  EEstrpen9;
  Float_t  EEstrpen10;
  Float_t  EEstrpen11;
  Float_t  EEstrpen12;
  Float_t  EEstrpen13;
  Float_t  EEstrpen14;
  Float_t  EEstrpen15;
  Int_t    ETwrdidE;
  Float_t  EPstripenp01;
  Float_t  EPstripenp02;
  Float_t  EPstripenp03;
  Float_t  EPstripenp0;
  Float_t  EPstripenp1;
  Float_t  EPstripenp2;
  Float_t  EPstripenp3;
  Float_t  EPstripenp4;
  Float_t  EPstripenp5;
  Float_t  EPstripenp6;
  Float_t  EPstripenp7;
  Float_t  EPstripenp8;
  Float_t  EPstripenp9;
  Float_t  EPstripenp10;
  Float_t  EPstripenp11;
  Float_t  EPstripenp12;
  Float_t  EPstripenp13;
  Float_t  EPstripenp14;
  Float_t  EPstripenp15;
  Float_t  EclustEnnq1;
  Float_t  EclustEnnq20;
  Float_t  EclustEnnq19;
  Float_t  EclustEnpq1;
  Float_t  EclustEnpq20;
  Float_t  EclustEnpq19;
  Float_t  EclustEnpq21;
  Int_t    PrimaryTrackArray_;
  UInt_t   PrimaryTrackArray_fUniqueID[nTrkMax];
  UInt_t   PrimaryTrackArray_fBits[nTrkMax];
  Double_t PrimaryTrackArray_nHitsFit[nTrkMax];
  Double_t PrimaryTrackArray_nHitsPoss[nTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[nTrkMax];
  Double_t PrimaryTrackArray_pZ[nTrkMax];
  Double_t PrimaryTrackArray_pX[nTrkMax];
  Double_t PrimaryTrackArray_pY[nTrkMax];
  Double_t PrimaryTrackArray_pT[nTrkMax];
  Double_t PrimaryTrackArray_dEdx[nTrkMax];
  Double_t PrimaryTrackArray_charge[nTrkMax];
  Double_t PrimaryTrackArray_tofBeta[nTrkMax];
  Double_t PrimaryTrackArray_eta[nTrkMax];
  Double_t PrimaryTrackArray_phi[nTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[nTrkMax];
  Double_t PrimaryTrackArray_nSigPion[nTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[nTrkMax];
  Double_t PrimaryTrackArray_nSigProton[nTrkMax];
  Double_t PrimaryTrackArray_dcag[nTrkMax];
  Double_t PrimaryTrackArray_nHits[nTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[nTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[nTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[nTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[nTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[nTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[nTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[nTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[nTrkMax];
  Double_t PrimaryTrackArray_pathLength[nTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[nTrkMax];
  Int_t    TowerArray_;
  UInt_t   TowerArray_fUniqueID[nTwrMax];
  UInt_t   TowerArray_fBits[nTwrMax];
  Int_t    TowerArray_TwrId[nTwrMax];
  Float_t  TowerArray_TwrEng[nTwrMax];
  Float_t  TowerArray_TwrEta[nTwrMax];
  Float_t  TowerArray_TwrPhi[nTwrMax];
  Float_t  TowerArray_TwrADC[nTwrMax];
  Float_t  TowerArray_TwrPed[nTwrMax];
  Float_t  TowerArray_TwrRMS[nTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[nTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[nTwrMax];
  Float_t  TowerArray_TwrMatchP[nTwrMax];
  Float_t  TowerArray_TwrPx[nTwrMax];
  Float_t  TowerArray_TwrPy[nTwrMax];
  Float_t  TowerArray_TwrPz[nTwrMax];
  Int_t    TowerArray_fNAssocTracks[nTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_P[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigPi[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigK[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigP[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigE[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_dcag[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_eta[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_pT[nTwrMax][nMatMax];
  Int_t    TowerArray_fMatchedTracksArray_nFit[nTwrMax][nMatMax];
  Int_t    TowerArray_fMatchedTracksArray_nPos[nTwrMax][nMatMax];

  // set branch addresses
  tInput -> SetMakeClass(1);
  tInput -> SetBranchAddress("fUniqueID", &fUniqueID);
  tInput -> SetBranchAddress("fBits", &fBits);
  tInput -> SetBranchAddress("runNumber", &runNumber);
  tInput -> SetBranchAddress("eventNumber", &eventNumber);
  tInput -> SetBranchAddress("trigID", &trigID);
  tInput -> SetBranchAddress("nGlobalTracks", &nGlobalTracks);
  tInput -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks);
  tInput -> SetBranchAddress("refMult", &refMult);
  tInput -> SetBranchAddress("vpdVz", &vpdVz);
  tInput -> SetBranchAddress("xVertex", &xVertex);
  tInput -> SetBranchAddress("yVertex", &yVertex);
  tInput -> SetBranchAddress("zVertex", &zVertex);
  tInput -> SetBranchAddress("bbcZVertex", &bbcZVertex);
  tInput -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate);
  tInput -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate);
  tInput -> SetBranchAddress("backgroundRate", &backgroundRate);
  tInput -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate);
  tInput -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate);
  tInput -> SetBranchAddress("refMultPos", &refMultPos);
  tInput -> SetBranchAddress("refMultNeg", &refMultNeg);
  tInput -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity);
  tInput -> SetBranchAddress("nVerticies", &nVerticies);
  tInput -> SetBranchAddress("MagF", &MagF);
  tInput -> SetBranchAddress("VrtxRank", &VrtxRank);
  tInput -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch);
  tInput -> SetBranchAddress("Etsp", &Etsp);
  tInput -> SetBranchAddress("ETwrdidT", &ETwrdidT);
  tInput -> SetBranchAddress("ETwradc11", &ETwradc11);
  tInput -> SetBranchAddress("ETwreneT0", &ETwreneT0);
  tInput -> SetBranchAddress("ETwreT", &ETwreT);
  tInput -> SetBranchAddress("ETwrENET0", &ETwrENET0);
  tInput -> SetBranchAddress("ETwrphT", &ETwrphT);
  tInput -> SetBranchAddress("ETwrPTower", &ETwrPTower);
  tInput -> SetBranchAddress("ETwrpidTower", &ETwrpidTower);
  tInput -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT);
  tInput -> SetBranchAddress("EClustEneT0", &EClustEneT0);
  tInput -> SetBranchAddress("EClustetav1", &EClustetav1);
  tInput -> SetBranchAddress("EClustphiv1", &EClustphiv1);
  tInput -> SetBranchAddress("EEstrpen01", &EEstrpen01);
  tInput -> SetBranchAddress("EEstrpen02", &EEstrpen02);
  tInput -> SetBranchAddress("EEstrpen03", &EEstrpen03);
  tInput -> SetBranchAddress("EEstrpen0", &EEstrpen0);
  tInput -> SetBranchAddress("EEstrpen1", &EEstrpen1);
  tInput -> SetBranchAddress("EEstrpen2", &EEstrpen2);
  tInput -> SetBranchAddress("EEstrpen3", &EEstrpen3);
  tInput -> SetBranchAddress("EEstrpen4", &EEstrpen4);
  tInput -> SetBranchAddress("EEstrpen5", &EEstrpen5);
  tInput -> SetBranchAddress("EEstrpen6", &EEstrpen6);
  tInput -> SetBranchAddress("EEstrpen7", &EEstrpen7);
  tInput -> SetBranchAddress("EEstrpen8", &EEstrpen8);
  tInput -> SetBranchAddress("EEstrpen9", &EEstrpen9);
  tInput -> SetBranchAddress("EEstrpen10", &EEstrpen10);
  tInput -> SetBranchAddress("EEstrpen11", &EEstrpen11);
  tInput -> SetBranchAddress("EEstrpen12", &EEstrpen12);
  tInput -> SetBranchAddress("EEstrpen13", &EEstrpen13);
  tInput -> SetBranchAddress("EEstrpen14", &EEstrpen14);
  tInput -> SetBranchAddress("EEstrpen15", &EEstrpen15);
  tInput -> SetBranchAddress("ETwrdidE", &ETwrdidE);
  tInput -> SetBranchAddress("EPstripenp01", &EPstripenp01);
  tInput -> SetBranchAddress("EPstripenp02", &EPstripenp02);
  tInput -> SetBranchAddress("EPstripenp03", &EPstripenp03);
  tInput -> SetBranchAddress("EPstripenp0", &EPstripenp0);
  tInput -> SetBranchAddress("EPstripenp1", &EPstripenp1);
  tInput -> SetBranchAddress("EPstripenp2", &EPstripenp2);
  tInput -> SetBranchAddress("EPstripenp3", &EPstripenp3);
  tInput -> SetBranchAddress("EPstripenp4", &EPstripenp4);
  tInput -> SetBranchAddress("EPstripenp5", &EPstripenp5);
  tInput -> SetBranchAddress("EPstripenp6", &EPstripenp6);
  tInput -> SetBranchAddress("EPstripenp7", &EPstripenp7);
  tInput -> SetBranchAddress("EPstripenp8", &EPstripenp8);
  tInput -> SetBranchAddress("EPstripenp9", &EPstripenp9);
  tInput -> SetBranchAddress("EPstripenp10", &EPstripenp10);
  tInput -> SetBranchAddress("EPstripenp11", &EPstripenp11);
  tInput -> SetBranchAddress("EPstripenp12", &EPstripenp12);
  tInput -> SetBranchAddress("EPstripenp13", &EPstripenp13);
  tInput -> SetBranchAddress("EPstripenp14", &EPstripenp14);
  tInput -> SetBranchAddress("EPstripenp15", &EPstripenp15);
  tInput -> SetBranchAddress("EclustEnnq1", &EclustEnnq1);
  tInput -> SetBranchAddress("EclustEnnq20", &EclustEnnq20);
  tInput -> SetBranchAddress("EclustEnnq19", &EclustEnnq19);
  tInput -> SetBranchAddress("EclustEnpq1", &EclustEnpq1);
  tInput -> SetBranchAddress("EclustEnpq20", &EclustEnpq20);
  tInput -> SetBranchAddress("EclustEnpq19", &EclustEnpq19);
  tInput -> SetBranchAddress("EclustEnpq21", &EclustEnpq21);
  tInput -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_);
  tInput -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID);
  tInput -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss);
  tInput -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag);
  tInput -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ);
  tInput -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX);
  tInput -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY);
  tInput -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx);
  tInput -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta);
  tInput -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta);
  tInput -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight);
  tInput -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength);
  tInput -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex);
  tInput -> SetBranchAddress("TowerArray", &TowerArray_);
  tInput -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID);
  tInput -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits);
  tInput -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId);
  tInput -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng);
  tInput -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta);
  tInput -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi);
  tInput -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC);
  tInput -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed);
  tInput -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS);
  tInput -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex);
  tInput -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk);
  tInput -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP);
  tInput -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx);
  tInput -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy);
  tInput -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz);
  tInput -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigPi[10]", TowerArray_fMatchedTracksArray_nSigPi);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigK[10]", TowerArray_fMatchedTracksArray_nSigK);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigP[10]", TowerArray_fMatchedTracksArray_nSigP);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigE[10]", TowerArray_fMatchedTracksArray_nSigE);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_dcag[10]", TowerArray_fMatchedTracksArray_dcag);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_eta[10]", TowerArray_fMatchedTracksArray_eta);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_pT[10]", TowerArray_fMatchedTracksArray_pT);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nFit[10]", TowerArray_fMatchedTracksArray_nFit);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_nPos[10]", TowerArray_fMatchedTracksArray_nPos);
  cout << "    Branches set." << endl;


  // load hot towers
  //UInt_t hotTwrs[nHotTwrs] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};
  UInt_t hotTwrs[nHotTwrs] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  cout << "    Hot towers loaded:\n"
       << "      " << nHotTwrs << " hot towers."
       << endl;


  fOutput -> cd();
  // trigger/event histograms
  TH1D *hNumTrg;
  TH1D *hNumTrk;
  TH1D *hTrgEt;
  TH1D *hTrgEne;
  TH1D *hTrgEta;
  TH1D *hTrgPhi;
  TH1D *hClustEta;
  TH1D *hClustPhi;
  TH1D *hTrgTsp;
  TH1D *hEffPt[nEffPar];
  TH2D *hTrgVsClustEta;
  TH2D *hTrgPhiVsEtaTwr;
  TH2D *hTrgPhiVsEtaClust;
  TH2D *hHotTrgPhiVsEta;
  TH2D *hHotTwrPhiVsEta;
  // track (triggered) histograms
  TH1D *hTrkEne[nDist];
  TH1D *hTrkPt[nDist];
  TH1D *hTrkEta[nDist];
  TH1D *hTrkPhi[nDist];
  TH1D *hTrkDeta[nDist];
  TH1D *hTrkDphi[nDist];
  TH2D *hTrkPhiVsEta;
  TH2D *hTrkPtVsEta;
  TH2D *hTrkPtVsPhi;
  TH2D *hTrkPtVsDeta;
  TH2D *hTrkPtVsDphi;
  // tower (triggered) histograms
  TH1D *hTwrEne[nDist];
  TH1D *hTwrEneH09[nDist];
  TH1D *hTwrCorr[nDist];
  TH1D *hTwrPt[nDist];
  TH1D *hTwrPtH09[nDist];
  TH1D *hTwrEta[nDist];
  TH1D *hTwrEtaH09[nDist];
  TH1D *hTwrPhi[nDist];
  TH1D *hTwrSumP[nDist];
  TH1D *hTwrDeta[nDist];
  TH1D *hTwrDphi[nDist];
  TH2D *hTwrPhiVsEta;
  TH2D *hTwrPhiVsEtaG10;
  TH2D *hTwrPtVsEta;
  TH2D *hTwrPtVsPhi;
  TH2D *hTwrPtVsDeta;
  TH2D *hTwrPtVsDphi;
  // tower E/p (triggered) histograms
  TH1D *hTwrEp[nSpecies][nNumBins];
  TH1D *hTwrEpG2[nSpecies][nNumBins];
  TH1D *hTwrTrkP[nSpecies];
  TH1D *hTwrPe[nNumBins];
  TH1D *hTwrNtrk;

  const UInt_t   nN  = 200;
  const UInt_t   nS  = 1000;
  const UInt_t   nE  = 200;
  const UInt_t   nH  = 80;
  const UInt_t   nDh = 160;
  const UInt_t   nF  = 120;
  const UInt_t   nDf = 240;
  const UInt_t   nP  = 100;
  const Double_t n[2]  = {0., 200.};
  const Double_t s[2]  = {0., 1.};
  const Double_t e[2]  = {0., 100.};
  const Double_t h[2]  = {-2., 2.};
  const Double_t dH[2] = {-4., 4.};
  const Double_t f[2]  = {-1. * TMath::Pi(), TMath::Pi()};
  const Double_t dF[2] = {-1. * TMath::TwoPi(), TMath::TwoPi()};
  const Double_t p[2]  = {0., 10.};
  // trigger histograms
  hNumTrg           = new TH1D("hNumTrg", "No. of triggers", nN, n[0], n[1]);
  hNumTrk           = new TH1D("hNumTrk", "No. of tracks", nN, n[0], n[1]);
  hTrgEt            = new TH1D("hTrgEt", "Trigger E_{T}", nE, e[0], e[1]);
  hTrgEne           = new TH1D("hTrgEne", "Trigger energy", nE, e[0], e[1]);
  hTrgEta           = new TH1D("hTrgEta", "Trigger (tower) #eta", nH, h[0], h[1]);
  hTrgPhi           = new TH1D("hTrgPhi", "Trigger (tower) #varphi", nF, f[0], f[1]);
  hClustEta         = new TH1D("hClustEta", "Trigger (cluster) #eta", nH, h[0], h[1]);
  hClustPhi         = new TH1D("hClustPhi", "Trigger (cluster) #varphi", nF, f[0], f[1]);
  hTrgTsp           = new TH1D("hTrgTsp", "Trigger TSP", nS, s[0], s[1]);
  hEffPt[0]         = new TH1D("hPionPt", "#pi^{#pm} p_{T}", nE, e[0], e[1]);
  hEffPt[1]         = new TH1D("hKaonPt", "K^{#pm} p_{T}", nE, e[0], e[1]);
  hEffPt[2]         = new TH1D("hProtonPt", "p^{#pm} p_{T}", nE, e[0], e[1]);
  hEffPt[3]         = new TH1D("hElectronPt", "e^{#pm} p_{T}", nE, e[0], e[1]);
  hEffPt[4]         = new TH1D("hTotalPt", "All particles p_{T}", nE, e[0], e[1]);
  hTrgVsClustEta    = new TH2D("hTrgVsClustEta", "Trigger tower vs. cluster eta;#eta^{twr};#eta^{clust}", nH, h[0], h[1], nH, h[0], h[1]);
  hTrgPhiVsEtaTwr   = new TH2D("hTrgPhiVsEtaTwr", "Trigger #varphi vs #eta_{twr}", nH ,h[0], h[1], nF, f[0], f[1]);
  hTrgPhiVsEtaClust = new TH2D("hTrgPhiVsEtaClust", "Trigger #varphi vs #eta_{clust}", nH, h[0], h[1], nF, f[0], f[1]);
  hHotTrgPhiVsEta   = new TH2D("hHotTrgPhiVsEta", "Hot trigger #varphi vs. #eta", nH, h[0], h[1], nF, f[0], f[1]);
  hHotTwrPhiVsEta   = new TH2D("hHotTwrPhiVsEta", "Hot tower #varphi vs. #eta", nH, h[0], h[1], nF, f[0], f[1]);
  // triggered track histograms
  hTrkEne[0]        = new TH1D("hTrkEneAll", "Track energy, all", nE, e[0], e[1]);
  hTrkEne[1]        = new TH1D("hTrkEneRec", "Track energy, recoil", nE, e[0], e[1]);
  hTrkPt[0]         = new TH1D("hTrkPtAll", "Track p_{T}, all", nE, e[0], e[1]);
  hTrkPt[1]         = new TH1D("hTrkPtRec", "Track p_{T}, recoil", nE, e[0], e[1]);
  hTrkEta[0]        = new TH1D("hTrkEtaAll", "Track #eta, all", nH, h[0], h[1]);
  hTrkEta[1]        = new TH1D("hTrkEtaRec", "Track #eta, recoil", nH, h[0], h[1]);
  hTrkPhi[0]        = new TH1D("hTrkPhiAll", "Track #varphi, all", nF, f[0], f[1]);
  hTrkPhi[1]        = new TH1D("hTrkPhiRec", "Track #varphi, recoil", nF, f[0], f[1]);
  hTrkDeta[0]       = new TH1D("hTrkDetaAll", "Track #Delta#eta, all", nDh, dH[0], dH[1]);
  hTrkDeta[1]       = new TH1D("hTrkDetaRec", "Track #Delta#eta, recoil", nDh, dH[0],dH[1]);
  hTrkDphi[0]       = new TH1D("hTrkDphiAll", "Track #Delta#varphi, all", nDf, dF[0], dF[1]);
  hTrkDphi[1]       = new TH1D("hTrkDphiRec", "Track #Delta#varphi, recoil", nDf, dF[0], dF[1]);
  hTrkPhiVsEta      = new TH2D("hTrkPhiVsEta", "Track #varphi vs. #eta", nH, h[0], h[1], nF, f[0], f[1]);
  hTrkPtVsEta       = new TH2D("hTrkPtVsEta", "Track p_{T} vs. #eta", nH, h[0], h[1], nE, e[0], e[1]);
  hTrkPtVsPhi       = new TH2D("hTrkPtVsPhi", "Track p_{T} vs. #varphi", nF, f[0], f[1], nE, e[0], e[1]);
  hTrkPtVsDeta      = new TH2D("hTrkPtVsDeta", "Track p_{T} vs. #Delta#eta", nDh, dH[0], dH[1], nE, e[0], e[1]);
  hTrkPtVsDphi      = new TH2D("hTrkPtVsDphi", "Track p_{T} vs. #Delta#varphi", nDf, dF[0], dF[1], nE, e[0], e[1]);
  // triggered tower histograms
  hTwrEne[0]        = new TH1D("hTwrEneAll", "Tower energy (raw), all", nE, e[0], e[1]);
  hTwrEne[1]        = new TH1D("hTwrEneRec", "Tower energy (raw), recoil", nE, e[0], e[1]);
  hTwrEneH09[0]     = new TH1D("hTwrEneH09all", "Tower energy (raw), all (|#eta^{clust}| < 0.9)", nE, e[0], e[1]);
  hTwrEneH09[1]     = new TH1D("hTwrEneH09rec", "Tower energy (raw), recoil (|#eta^{clust}| < 0.9)", nE, e[0], e[1]);
  hTwrCorr[0]       = new TH1D("hTwrCorrAll", "Tower energy (corr.), all", nE, e[0], e[1]);
  hTwrCorr[1]       = new TH1D("hTwrCorrRec", "Tower energy (corr.), recoil", nE, e[0], e[1]);
  hTwrPt[0]         = new TH1D("hTwrPtAll", "Tower p_{T}, all", nE, e[0], e[1]);
  hTwrPt[1]         = new TH1D("hTwrPtRec", "Tower p_{T}, recoil", nE, e[0], e[1]);
  hTwrPtH09[0]      = new TH1D("hTwrPtH09all", "Tower p_{T}, all (|#eta^{clust}| < 0.9)", nE, e[0], e[1]);
  hTwrPtH09[1]      = new TH1D("hTwrPtH09rec", "Tower p_{T}, recoil (|#eta^{clust}| < 0.9)", nE, e[0], e[1]);
  hTwrSumP[0]       = new TH1D("hTwrSumPall", "Tower #sump_{match}, all", nE, e[0], e[1]);
  hTwrSumP[1]       = new TH1D("hTwrSumPrec", "Tower #sump_{match}, recoil", nE, e[0], e[1]);
  hTwrEta[0]        = new TH1D("hTwrEtaAll", "Tower #eta, all", nH, h[0], h[1]);
  hTwrEta[1]        = new TH1D("hTwrEtaRec", "Tower #eta, recoil", nH, h[0], h[1]);
  hTwrEtaH09[0]     = new TH1D("hTwrEtaH09all", "Tower #eta, all (|#eta^{clust}| < 0.9)", nH, h[0], h[1]);
  hTwrEtaH09[1]     = new TH1D("hTwrEtaH09rec", "Tower #eta, recoil (|#eta^{clust}| < 0.9)", nH, h[0], h[1]);
  hTwrPhi[0]        = new TH1D("hTwrPhiAll", "Tower #varphi, all", nF, f[0], f[1]);
  hTwrPhi[1]        = new TH1D("hTwrPhiRec", "Tower #varphi, recoil", nF, f[0], f[1]);
  hTwrDeta[0]       = new TH1D("hTwrDetaAll", "Tower #Delta#eta, all", nDh, dH[0], dH[1]);
  hTwrDeta[1]       = new TH1D("hTwrDetaRec", "Tower #Delta#eta, recoil", nDh, dH[0], dH[1]);
  hTwrDphi[0]       = new TH1D("hTwrDphiAll", "Tower #Delta#varphi, all", nDf, dF[0], dF[1]);
  hTwrDphi[1]       = new TH1D("hTwrDphiRec", "Tower #Delta#varphi, recoil", nDf, dF[0], dF[1]);
  hTwrPhiVsEta      = new TH2D("hTwrPhiVsEta", "Tower #varphi vs. #eta", nH, h[0], h[1], nF, f[0], f[1]);
  hTwrPhiVsEtaG10   = new TH2D("hTwrPhiVsEtaG10", "Tower #varphi vs. #eta, E_{T}^{twr} > 10 GeV", nH, h[0], h[1], nF, f[0], f[1]);
  hTwrPtVsEta       = new TH2D("hTwrPtVsEta", "Tower p_{T} vs. #eta", nH, h[0], h[1], nE, e[0], e[1]);
  hTwrPtVsPhi       = new TH2D("hTwrPtVsPhi", "Tower p_{T} vs. #varphi", nF, f[0], f[1], nE, e[0], e[1]);
  hTwrPtVsDeta      = new TH2D("hTwrPtVsDeta", "Tower p_{T} vs. #Delta#eta", nDh, dH[0], dH[1], nE, e[0], e[1]);
  hTwrPtVsDphi      = new TH2D("hTwrPtVsDphi", "Tower p_{T} vs. #Delta#varphi", nDf, dF[0], dF[1], nE, e[0], e[1]);
  // tower E/p (triggered) histograms
  hTwrEp[0][0]      = new TH1D("hTwrEp_elecN1", "Tower E/p (e^{#pm}), N_{match} = 1", nP, p[0], p[1]);
  hTwrEp[0][1]      = new TH1D("hTwrEp_elecN2", "Tower E/p (e^{#pm}), N_{match} = 2", nP, p[0], p[1]);
  hTwrEp[0][2]      = new TH1D("hTwrEp_elecN3", "Tower E/p (e^{#pm}), N_{match} #geq 3", nP, p[0], p[1]);
  hTwrEp[0][3]      = new TH1D("hTwrEp_elecNA", "Tower E/p (e^{#pm}), N_{match} > 0", nP, p[0], p[1]);
  hTwrEp[1][0]      = new TH1D("hTwrEp_hadrN1", "Tower E/p (h^{#pm}), N_{match} = 1", nP, p[0], p[1]);
  hTwrEp[1][1]      = new TH1D("hTwrEp_hadrN2", "Tower E/p (h^{#pm}), N_{match} = 2", nP, p[0], p[1]);
  hTwrEp[1][2]      = new TH1D("hTwrEp_hadrN3", "Tower E/p (h^{#pm}), N_{match} #geq 3", nP, p[0], p[1]);
  hTwrEp[1][3]      = new TH1D("hTwrEp_hadrNA", "Tower E/p (h^{#pm}), N_{match} > 0", nP, p[0], p[1]);
  hTwrEpG2[0][0]    = new TH1D("hTwrEpG2_elecN1", "Tower E/p (e^{#pm}), N_{match} = 1 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[0][1]    = new TH1D("hTwrEpG2_elecN2", "Tower E/p (e^{#pm}), N_{match} = 2 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[0][2]    = new TH1D("hTwrEpG2_elecN3", "Tower E/p (e^{#pm}), N_{match} #geq 3 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[0][3]    = new TH1D("hTwrEpG2_elecNA", "Tower E/p (e^{#pm}), N_{match} > 0 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[1][0]    = new TH1D("hTwrEpG2_hadrN1", "Tower E/p (h^{#pm}), N_{match} = 1 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[1][1]    = new TH1D("hTwrEpG2_hadrN2", "Tower E/p (h^{#pm}), N_{match} = 2 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[1][2]    = new TH1D("hTwrEpG2_hadrN3", "Tower E/p (h^{#pm}), N_{match} #geq 3 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrEpG2[1][3]    = new TH1D("hTwrEpG2_hadrNA", "Tower E/p (h^{#pm}), N_{match} > 0 and p_{T}^{match} > 2 GeV/c", nP, p[0], p[1]);
  hTwrTrkP[0]       = new TH1D("hTwrTrkP_elec", "Matched track p (e^{#pm})", nE, e[0], e[1]);
  hTwrTrkP[1]       = new TH1D("hTwrTrkP_hadr", "Matched track p (h^{#pm})", nE, e[0], e[1]);
  hTwrPe[0]         = new TH1D("hTwrPeN1", "Tower #sump/E, N_{match} = 1", nP, p[0], p[1]);
  hTwrPe[1]         = new TH1D("hTwrPeN2", "Tower #sump/E, N_{match} = 2", nP, p[0], p[1]);
  hTwrPe[2]         = new TH1D("hTwrPeN3", "Tower #sump/E, N_{match} #geq 3", nP, p[0], p[1]);
  hTwrPe[3]         = new TH1D("hTwrPeNA", "Tower #sump/E, N_{match} > 0", nP, p[0], p[1]);
  hTwrNtrk          = new TH1D("hTwrNtrk", "No. of matched tracks", nN, n[0], n[1]);
  // error
  hNumTrk           -> Sumw2();
  hTrgEt            -> Sumw2();
  hTrgEne           -> Sumw2();
  hTrgEta           -> Sumw2();
  hTrgPhi           -> Sumw2();
  hClustEta         -> Sumw2();
  hClustPhi         -> Sumw2();
  hTrgTsp           -> Sumw2();
  hTrgVsClustEta    -> Sumw2();
  hTrgPhiVsEtaTwr   -> Sumw2();
  hTrgPhiVsEtaClust -> Sumw2();
  hHotTrgPhiVsEta   -> Sumw2();
  hHotTwrPhiVsEta   -> Sumw2();
  for (Long64_t iEffPar = 0; iEffPar < nEffPar; iEffPar++) {
    hEffPt[iEffPar] -> Sumw2();
  }
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    hTrkEne[iDist]    -> Sumw2();
    hTrkPt[iDist]     -> Sumw2();
    hTrkEta[iDist]    -> Sumw2();
    hTrkPhi[iDist]    -> Sumw2();
    hTrkDeta[iDist]   -> Sumw2();
    hTrkDphi[iDist]   -> Sumw2();
    hTwrEne[iDist]    -> Sumw2();
    hTwrEneH09[iDist] -> Sumw2();
    hTwrCorr[iDist]   -> Sumw2();
    hTwrPt[iDist]     -> Sumw2();
    hTwrPtH09[iDist]  -> Sumw2();
    hTwrSumP[iDist]   -> Sumw2();
    hTwrEta[iDist]    -> Sumw2();
    hTwrEtaH09[iDist] -> Sumw2();
    hTwrPhi[iDist]    -> Sumw2();
    hTwrDeta[iDist]   -> Sumw2();
    hTwrDphi[iDist]   -> Sumw2();
  }
  for (Long64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (Long64_t iNumBin = 0; iNumBin < nNumBins; iNumBin++) {
      hTwrEp[iSpecies][iNumBin]   -> Sumw2();
      hTwrEpG2[iSpecies][iNumBin] -> Sumw2();
    }
    hTwrTrkP[iSpecies] -> Sumw2();
  }
  for (Long64_t iNumBin = 0; iNumBin < nNumBins; iNumBin++) {
    hTwrPe[iNumBin] -> Sumw2();
  }
  hTwrNtrk        -> Sumw2();
  hTrkPhiVsEta    -> Sumw2();
  hTrkPtVsEta     -> Sumw2();
  hTrkPtVsPhi     -> Sumw2();
  hTrkPtVsDeta    -> Sumw2();
  hTrkPtVsDphi    -> Sumw2();
  hTwrPhiVsEta    -> Sumw2();
  hTwrPhiVsEtaG10 -> Sumw2();
  hTwrPtVsEta     -> Sumw2();
  hTwrPtVsPhi     -> Sumw2();
  hTwrPtVsDeta    -> Sumw2();
  hTwrPtVsDphi    -> Sumw2();


  Long64_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " evts. to process." << endl;

  // event loop
  Int_t    iTrig  = 0;
  Long64_t nTrg   = 0;
  Long64_t nBytes = 0;
  Double_t fTrig  = 0;
  Double_t hTrig  = 0;
  for (Long64_t iEvt = 0; iEvt < nEvts; iEvt++) {

    const Long64_t bytes = tInput -> GetEntry(iEvt);
    if (bytes < 0) {
      cerr << "WARNING: Something weird at event " << iEvt + 1 << "!" << endl;
      break;
    }
    nBytes += bytes;

    if (!isInBatchMode) {
      cout << "      processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
      if ((iEvt + 1) == nEvts) cout << endl;
    } 
    else {
      cout << "      processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
    }


    // reset trigger values
    iTrig = 0;
    fTrig = 0;
    hTrig = 0;

    // event info
    const Int_t    nTrks = nPrimaryTracks;
    const Int_t    nTwrs = TowerArray_;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isInZCut = (TMath::Abs(zVtx) < zVtxMax);
    const Bool_t isInRCut = (TMath::Abs(rVtx) < rVtxMax);
    if (!isInRCut || !isInZCut) continue;


    // trigger info
    const Int_t    trgID  = ETwrdidT;
    const Int_t    adc    = ETwradc11;
    const Float_t  tsp    = Etsp;
    const Double_t eHstrp = EEstrpen4;
    const Double_t eFstrp = EPstripenp4;
    const Double_t pProj  = ETwrPTower;
    const Double_t hTrg   = ETwreT;
    const Double_t fTrg   = ETwrphT;
    const Double_t hClust = EClustetav1;
    const Double_t fClust = EClustphiv1;
    const Double_t eTrg   = EClustEneT0;
    const Double_t eTtrg  = eTrg * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hClust)));

    // trigger cuts
    if (isInTrgMode1) {
      // hot tower check
      Bool_t isHotTrg = false;
      if (removeHotTwrs) {
        for (Long64_t iHot = 0; iHot < nHotTwrs; iHot++) {
          if (trgID == hotTwrs[iHot]) {
            isHotTrg = true;
            hHotTrgPhiVsEta -> Fill(hTrg, fTrg);
            break;
          }
        }
      }  // end hot tower check

      // TEST [02.23.2018]
      const Bool_t isInBadPhiZone = ((fTrg > 1.25) && (fTrg < 1.58));
      if (isInBadPhiZone) continue;

      const Bool_t isInAdcCut      = (adc <= adcMax);
      const Bool_t isInStrpCut     = ((eHstrp >= eHmin) && (eFstrp >= eFmin));
      const Bool_t isInProjCut     = (pProj < pProjMax);
      const Bool_t isInEtaTrgCut   = (TMath::Abs(hTrg) < hTrgMax);
      const Bool_t isInEtaClustCut = (TMath::Abs(hClust) < hTrgMax);
      const Bool_t isInEneTrgCut   = ((eTtrg > eTtrgMin) && (eTtrg < eTtrgMax));
      const Bool_t isInTspTrgCut   = ((tsp > tspMin) && (tsp < tspMax));
      const Bool_t foundTrg1       = (isInAdcCut && isInStrpCut && isInProjCut && isInEtaTrgCut && isInEneTrgCut && isInTspTrgCut && !isHotTrg);
      if (!foundTrg1) {
        continue;
      }
      else {
        hNumTrk           -> Fill(nTrks);
        hTrgEt            -> Fill(eTtrg);
        hTrgEne           -> Fill(eTrg);
        hTrgEta           -> Fill(hTrg);
        hTrgPhi           -> Fill(fTrg);
        hClustEta         -> Fill(hClust);
        hClustPhi         -> Fill(fClust);
        hTrgTsp           -> Fill(tsp);
        hTrgVsClustEta    -> Fill(hTrg, hClust);
        hTrgPhiVsEtaTwr   -> Fill(hTrg, fTrg);
        hTrgPhiVsEtaClust -> Fill(hClust, fClust);
        iTrig = trgID;
        fTrig = fClust;
        hTrig = hClust;
        nTrg++;
      }
    }  // end if (isInTrgMode1)


    // tower loop 1
    UInt_t numTrg = 0;
    if (isInTrgMode3) {
      Bool_t foundTrg3 = false;
      for (Long64_t iTwr = 0; iTwr <  nTwrs; iTwr++) {

        // reset flag
        foundTrg3 = false;

        // tower info
        const Int_t    twrID  = TowerArray_TwrId[iTwr];
        const Int_t    twrADC = TowerArray_TwrADC[iTwr];
        const UInt_t   fMatch = TowerArray_TwrMatchIdnex[iTwr];
        const UInt_t   nMatch = TowerArray_NoOfmatchedTrk[iTwr];
        const Double_t pProj  = TowerArray_TwrMatchP[iTwr];
        const Double_t fTwr   = TowerArray_TwrPhi[iTwr];
        const Double_t hTwr   = TowerArray_TwrEta[iTwr];
        const Double_t eTwr   = TowerArray_TwrEng[iTwr];
        const Double_t eTtwr  = eTwr * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTwr)));
        const Double_t pXtwr  = TowerArray_TwrPx[iTwr];
        const Double_t pYtwr  = TowerArray_TwrPy[iTwr];
        const Double_t pTtwr  = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr));
        const Double_t pTwr   = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr) + (pZtwr * pZtwr));

        // calculate projectd phi, eta
        Double_t hClust = TMath::Log(1. / (TMath::Tan(0.5 * TMath::ASin(pTtwr / pTwr))));
        Double_t fClust = TMath::ATan(pYtwr / pXtwr);
        if (pZtwr < 0.) hClust *= -1.;

        // hot tower check
        if (removeHotTwrs) {
        Bool_t isHot = false;
          for (Long64_t iHot = 0; iHot < nHotTwrs; iHot++) {
            if (twrID == hotTwrs[iHot]) {
              isHot = true;
              break;
            }
          }
          if (isHot) continue;
        }

        // tower cuts
        const Bool_t isInAdcTrgCut  = (twrADC <= adcMax);
        const Bool_t isInProjTrgCut = (pProj < pProjMax);
        const Bool_t isTriggerTwr   = (twrID == trgID);
        const Bool_t isInEtaTwrCut  = (TMath::Abs(hTwr) < hTwrMax);
        const Bool_t isInEneTwrCut  = (eTwr > eTwrMin);
        const Bool_t isInEtTrgCut   = ((eTtwr > eTtrgMin) && (eTtwr < eTtrgMax));
        if (!isInAdcTrgCut || !isInProjTrgCut || !isInEneTrgCut || !isInEtaTrgCut || !isInEtTrgCut) continue;

        // check if trigger
        if (isInEtTrgCut) {
          iTrig     = twrID;
          fTrig     = fTwr;
          hTrig     = hTwr;
          foundTrg3 = true;
          numTrg++;
          hNumTrk           -> Fill(nTrks);
          hTrgEt            -> Fill(eTtwr);
          hTrgEne           -> Fill(eTwr);
          hTrgEta           -> Fill(hTwr);
          hTrgPhi           -> Fill(fTwr);
          hClustEta         -> Fill(hClust);
          hClustPhi         -> Fill(fClust);
          hTrgVsClustEta    -> Fill(hTwr, hClust);
          hTrgPhiVsEtaTwr   -> Fill(hTwr, fTwr);
          hTrgPhiVsEtaClust -> Fill(hClust, fClust);
          break;
        }

      }  // end tower loop 1
      if (!foundTrg3) {
        continue;
      }
      else {
        nTrg++;
      }
    }  // end if (isInTrgMode3)
    hNumTrg -> Fill(numTrg);


    // track loop
    for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    nFitTrk = PrimaryTrackArray_nHitsFit[iTrk];
      const Int_t    nPosTrk = PrimaryTrackArray_nHitsPoss[iTrk];
      const Float_t  nPiTrk  = TMath::Abs(PrimaryTrackArray_nSigPion[iTrk]);
      const Float_t  nKtrk   = TMath::Abs(PrimaryTrackArray_nSigKaon[iTrk]);
      const Float_t  nPtrk   = TMath::Abs(PrimaryTrackArray_nSigProton[iTrk]);
      const Float_t  nEtrk   = TMath::Abs(PrimaryTrackArray_nSigElectron[iTrk]);
      const Double_t rFitTrk = (Double_t) nFitTrk / (Double_t) nPosTrk;
      const Double_t dca     = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk    = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk    = PrimaryTrackArray_phi[iTrk];
      const Double_t pZtrk   = PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk   = PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk    = TMath::Sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (mPion * mPion));

      Double_t dHtrk = hTrk - hTrig;
      Double_t dFtrk = fTrk - fTrig;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();


      // track cuts
      const Bool_t isInNFitTrkCut = (nFitTrk > nFitMin);
      const Bool_t isInRFitTrkCut = (rFitTrk > rFitMin);
      const Bool_t isInDcaTrkCut  = (dca < dcaMax);
      const Bool_t isInEtaTrkCut  = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtTrkCut   = (pTtrk > pTtrkMin);
      if (!isInNFitTrkCut || !isInRFitTrkCut || !isInDcaTrkCut || !isInEtaTrkCut || !isInPtTrkCut) continue;

      // fill histograms
      const Bool_t isPion      = ((nPiTrk < pidCut) && (nKtrk > pidCut) && (nPtrk > pidCut) && (nEtrk > pidCut));
      const Bool_t isKaon      = ((nPiTrk > pidCut) && (nKtrk < pidCut) && (nPtrk > pidCut) && (nEtrk > pidCut));
      const Bool_t isProton    = ((nPiTrk > pidCut) && (nKtrk > pidCut) && (nPtrk < pidCut) && (nEtrk > pidCut));
      const Bool_t isElectron  = ((nPiTrk > pidCut) && (nKtrk > pidCut) && (nPtrk > pidCut) && (nEtrk < pidCut));
      const Bool_t isRecoilTrk = (TMath::Abs(dFtrk - TMath::Pi()) < fRecoilMax);
      hTrkEne[0]  -> Fill(eTrk);
      hTrkPt[0]   -> Fill(pTtrk);
      hTrkEta[0]  -> Fill(hTrk);
      hTrkPhi[0]  -> Fill(fTrk);
      hTrkDeta[0] -> Fill(dHtrk);
      hTrkDphi[0] -> Fill(dFtrk);
      if (isRecoilTrk) {
        hTrkEne[1]  -> Fill(eTrk);
        hTrkPt[1]   -> Fill(pTtrk);
        hTrkEta[1]  -> Fill(hTrk);
        hTrkPhi[1]  -> Fill(fTrk);
        hTrkDeta[1] -> Fill(dHtrk);
        hTrkDphi[1] -> Fill(dFtrk);
      }
      hTrkPhiVsEta -> Fill(hTrk, fTrk);
      hTrkPtVsEta  -> Fill(hTrk, pTtrk);
      hTrkPtVsPhi  -> Fill(fTrk, pTtrk);
      hTrkPtVsDeta -> Fill(dHtrk, pTtrk);
      hTrkPtVsDphi -> Fill(dFtrk, pTtrk);

      // for efficiency
      if (isPion)     hEffPt[0] -> Fill(pTtrk);
      if (isKaon)     hEffPt[1] -> Fill(pTtrk);
      if (isProton)   hEffPt[2] -> Fill(pTtrk);
      if (isElectron) hEffPt[3] -> Fill(pTtrk);
      hEffPt[4] -> Fill(pTtrk);

    }  // end track loop


    // tower loop 2
    for (Long64_t iTwr = 0; iTwr <  nTwrs; iTwr++) {

      // tower info
      const Int_t    twrID    = TowerArray_TwrId[iTwr];
      const UInt_t   fMatch   = TowerArray_TwrMatchIdnex[iTwr];
      const UInt_t   nMatch   = TowerArray_NoOfmatchedTrk[iTwr];
      const Double_t pProjTwr = TowerArray_TwrMatchP[iTwr];
      const Double_t fTwr     = TowerArray_TwrPhi[iTwr];
      const Double_t hTwr     = TowerArray_TwrEta[iTwr];
      const Double_t eTwr     = TowerArray_TwrEng[iTwr];
      const Double_t eTtwr    = eTwr * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * eTwr)));
      const Double_t pXtwr    = TowerArray_TwrPx[iTwr];
      const Double_t pYtwr    = TowerArray_TwrPy[iTwr];
      const Double_t pTtwr    = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr));

      Double_t dHtwr = hTwr - hTrig;
      Double_t dFtwr = fTwr - fTrig;
      if (dFtwr < (-1. * TMath::PiOver2())) dFtwr += TMath::TwoPi();
      if (dFtwr > (3. * TMath::PiOver2()))  dFtwr -= TMath::TwoPi();

      // hot tower check
      if (removeHotTwrs) {
        Bool_t isHot = false;
        for (Long64_t iHot = 0; iHot < nHotTwrs; iHot++) {
          if (twrID == hotTwrs[iHot]) {
            isHot = true;
            hHotTwrPhiVsEta -> Fill(hTwr, fTwr);
            break;
          }
        }
        if (isHot) continue;
      }  // end hot tower check

      // check if trigger
      const Bool_t isTrigger    = (twrID == iTrig);
      const Bool_t isTriggerTwr = (twrID == trgID);
      if (isTrigger || isTriggerTwr) continue;

      // tower cuts
      const Bool_t isInEtaTwrCut  = (TMath::Abs(hTwr) < hTwrMax);
      const Bool_t isInEneTwrCut  = (eTwr > eTwrMin);
      const Bool_t isInEneHardCut = (eTwr > eTwrHard);
      if (!isInEtaTwrCut || !isInEneTwrCut) continue;

      // fill histograms
      const Bool_t isRecoilTwr = (TMath::Abs(dFtwr - TMath::Pi()) < fRecoilMax);
      hTwrEne[0]  -> Fill(eTwr);
      hTwrPt[0]   -> Fill(pTtwr);
      hTwrEta[0]  -> Fill(hTwr);
      hTwrPhi[0]  -> Fill(fTwr);
      hTwrSumP[0] -> Fill(pProjTwr);
      hTwrDeta[0] -> Fill(dHtwr);
      hTwrDphi[0] -> Fill(dFtwr);
      if (isInEtaClustCut) {
        hTwrEneH09[0] -> Fill(eTwr);
        hTwrPtH09[0]  -> Fill(pTtwr);
        hTwrEtaH09[0] -> Fill(hTwr);
      }
      if (isRecoilTwr) {
        hTwrEne[1]  -> Fill(eTwr);
        hTwrPt[1]   -> Fill(pTtwr);
        hTwrEta[1]  -> Fill(hTwr);
        hTwrPhi[1]  -> Fill(fTwr);
        hTwrSumP[1] -> Fill(pProjTwr);
        hTwrDeta[1] -> Fill(dHtwr);
        hTwrDphi[1] -> Fill(dFtwr);
        if (isInEtaClustCut) {
          hTwrEneH09[1] -> Fill(eTwr);
          hTwrPtH09[1]  -> Fill(pTtwr);
          hTwrEtaH09[1] -> Fill(hTwr);
        }
      }
      hTwrNtrk     -> Fill(nMatch);
      hTwrPhiVsEta -> Fill(hTwr, fTwr);
      hTwrPtVsEta  -> Fill(hTwr, pTtwr);
      hTwrPtVsPhi  -> Fill(fTwr, pTtwr);
      hTwrPtVsDeta -> Fill(dHtwr, pTtwr);
      hTwrPtVsDphi -> Fill(dFtwr, pTtwr);
      if (isInEneHardCut) hTwrPhiVsEtaG10 -> Fill(hTwr, fTwr);

      // E/p calculation
      const Bool_t isMatched         = (fMatch == 1);
      const Bool_t hasMatches        = (nMatch >= 1);
      const Bool_t hasOneMatch       = (nMatch == 1);
      const Bool_t hasTwoMatches     = (nMatch == 2);
      const Bool_t hasSeveralMatches = (nMatch >= 3);
      if (isMatched && hasMatches) {
        for (Long64_t iMatch = 0; iMatch < nMatch; iMatch++) {
          const Int_t    nFitMat = TowerArray_fMatchedTracksArray_nFit[iTwr][iMatch];
          const Int_t    nPosMat = TowerArray_fMatchedTracksArray_nPos[iTwr][iMatch];
          const Float_t  nPiMat  = TMath::Abs(TowerArray_fMatchedTracksArray_nSigPi[iTwr][iMatch]);
          const Float_t  nKmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigK[iTwr][iMatch]);
          const Float_t  nPmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigP[iTwr][iMatch]);
          const Float_t  nEmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigE[iTwr][iMatch]);
          const Float_t  dcaMat  = TowerArray_fMatchedTracksArray_dcag[iTwr][iMatch];
          const Float_t  hMat    = TowerArray_fMatchedTracksArray_eta[iTwr][iMatch];
          const Float_t  pTmat   = TowerArray_fMatchedTracksArray_pT[iTwr][iMatch];
          const Float_t  pMat    = TowerArray_fMatchedTracksArray_P[iTwr][iMatch];
          const Double_t rFitMat = (Double_t) nFitMat / (Double_t) nPosMat;

          // match cuts
          const Bool_t isInElecCut = ((nPiMat > pidCut) && (nKmat > pidCut) && (nPmat > pidCut) && (nEmat < pidCut));
          const Bool_t isInHadrCut = ((nPiMat < pidCut) && (nKmat < pidCut) && (nPmat < pidCut) && (nEmat > pidCut));
          const Bool_t isInNFitCut = (nFitMat > nFitMin);
          const Bool_t isInRFitCut = (rFitMat > rFitMin);
          const Bool_t isInDcaCut  = (dcaMat < dcaMax);
          const Bool_t isInEtaCut  = (TMath::Abs(hMat) < hTrkMax);
          const Bool_t isInPtCut   = (pTmat > pTtrkMin);
          const Bool_t isInG2Cut   = (pTmat > 2.);
          if (!isInNFitCut || !isInRFitCut || !isInDcaCut || !isInEtaCut || !isInEtaCut || !isInPtCut) continue;

          // electrons
          if (isInElecCut) {
            if (hasOneMatch) {
              hTwrEp[0][0] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[0][0] -> Fill(eTwr / pMat);
            }
            if (hasTwoMatches) {
              hTwrEp[0][1] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[0][1] -> Fill(eTwr / pMat);
            }
            if (hasSeveralMatches) {
              hTwrEp[0][2] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[0][2] -> Fill(eTwr / pMat);
            }
            hTwrEp[0][3] -> Fill(eTwr / pMat);
            hTwrTrkP[0]  -> Fill(pMat);
            if (isInG2Cut) hTwrEpG2[0][3] -> Fill(eTwr / pMat);
          }
          // hadrons
          if (isInHadrCut) {
            if (hasOneMatch) {
              hTwrEp[1][0] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[1][0] -> Fill(eTwr / pMat);
            }
            if (hasTwoMatches) {
              hTwrEp[1][1] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[1][1] -> Fill(eTwr / pMat);
            }
            if (hasSeveralMatches) {
              hTwrEp[1][2] -> Fill(eTwr / pMat);
              if (isInG2Cut) hTwrEpG2[1][2] -> Fill(eTwr / pMat);
            }
            hTwrEp[1][3] -> Fill(eTwr / pMat);
            hTwrTrkP[1]  -> Fill(pMat);
            if (isInG2Cut) hTwrEpG2[1][3] -> Fill(eTwr / pMat);
          }
        }  // end match loop
        if (hasOneMatch)
          hTwrPe[0] -> Fill(pProjTwr / eTwr);
        if (hasTwoMatches)
          hTwrPe[1] -> Fill(pProjTwr / eTwr);
        if (hasSeveralMatches)
          hTwrPe[2] -> Fill(pProjTwr / eTwr);
        if (hasMatches)
          hTwrPe[3] -> Fill(pProjTwr / eTwr);
      }  // end E/p calculation

    }  // end tower loop 2

  }  // end event loop

  cout << "    Event loop finished: " << nTrg << " triggers found." << endl;


  // make directories
  TDirectory *dTrks[nDist + 1];
  TDirectory *dTwrs[nDist + 1];
  TDirectory *dEp;
  for (Long64_t iDist = 0; iDist < nDist + 1; iDist++) {
    dTrks[iDist] = (TDirectory*) fOutput -> mkdir(sTrkDirs[iDist].Data());
    dTwrs[iDist] = (TDirectory*) fOutput -> mkdir(sTwrDirs[iDist].Data());
  }
  dEp = (TDirectory*) fOutput -> mkdir(sEpDir.Data());
  cout << "    Directories made." << endl;


  fOutput           -> cd();
  hNumTrg           -> Write();
  hNumTrk           -> Write();
  hTrgEt            -> Write();
  hTrgEne           -> Write();
  hTrgEta           -> Write();
  hTrgPhi           -> Write();
  hClustEta         -> Write();
  hClustPhi         -> Write();
  hTrgTsp           -> Write();
  hTrgVsClustEta    -> Write();
  hTrgPhiVsEtaTwr   -> Write();
  hTrgPhiVsEtaClust -> Write();
  hHotTrgPhiVsEta   -> Write();
  hHotTwrPhiVsEta   -> Write();
  for (Long64_t iEffPar = 0; iEffPar < nEffPar; iEffPar++) {
    hEffPt[iEffPar] -> Write();
  }
  for (Long64_t iDist = 0; iDist < nDist + 1; iDist++) {
    dTrks[iDist] -> cd();
    if (iDist == nDist) {
      hTrkPhiVsEta -> Write();
      hTrkPtVsEta  -> Write();
      hTrkPtVsPhi  -> Write();
      hTrkPtVsDeta -> Write();
      hTrkPtVsDphi -> Write();
    }
    else {
      hTrkEne[iDist]  -> Write();
      hTrkPt[iDist]   -> Write();
      hTrkEta[iDist]  -> Write();
      hTrkPhi[iDist]  -> Write();
      hTrkDeta[iDist] -> Write();
      hTrkDphi[iDist] -> Write();
    }
    dTwrs[iDist] -> cd();
    if (iDist == nDist) {
      hTwrPhiVsEta    -> Write();
      hTwrPhiVsEtaG10 -> Write();
      hTwrPtVsEta     -> Write();
      hTwrPtVsPhi     -> Write();
      hTwrPtVsDeta    -> Write();
      hTwrPtVsDphi    -> Write();
    }
    else {
      hTwrEne[iDist]    -> Write();
      hTwrEneH09[iDist] -> Write();
      hTwrCorr[iDist]   -> Write();
      hTwrPt[iDist]     -> Write();
      hTwrPtH09[iDist]  -> Write();
      hTwrSumP[iDist]   -> Write();
      hTwrEta[iDist]    -> Write();
      hTwrEtaH09[iDist] -> Write();
      hTwrPhi[iDist]    -> Write();
      hTwrDeta[iDist]   -> Write();
      hTwrDphi[iDist]   -> Write();
    }
  }
  dEp -> cd();
  for (Long64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (Long64_t iNumBin = 0; iNumBin < nNumBins; iNumBin++) {
      hTwrEp[iSpecies][iNumBin]   -> Write();
      hTwrEpG2[iSpecies][iNumBin] -> Write();
    }
    hTwrTrkP[iSpecies] -> Write();
  }
  for (Long64_t iNumBin = 0; iNumBin < nNumBins; iNumBin++) {
    hTwrPe[iNumBin] -> Write();
  }
  hTwrNtrk -> Write();
  fOutput  -> Close();
  fInput   -> cd();
  fInput   -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
