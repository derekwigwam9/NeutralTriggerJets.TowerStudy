// 'SimpleSpectrumPi0.C'
// Derek Anderson
// 03.28.2018
//
// Use this to (simply) get a pi0-
// triggered spectrum from the
// output of the 'StJetTreeThird-
// Maker' class.


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"

using namespace std;


// global constants
static const UInt_t   NTrkMax(1000);
static const UInt_t   NTwrMax(10000);
static const UInt_t   NMatMax(100);
static const UInt_t   NHotTwr(41);
// io parameters
static const TString sTree("Gfmtodst");
static const TString sInDefault("/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.epCalc2.root");
static const TString sOutDefault("test.root");



void SimpleSpectrumPi0(const Bool_t isInBatchMode=false, const TString sInput=sInDefault, const TString sOutput=sOutDefault) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning (simple) spectrum script..." << endl;

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
  const Double_t pTtwrMin = 0.2;
  // delta-phi cut
  const Double_t dFrecoil(TMath::PiOver2());
  const Double_t dFcut1(TMath::Pi() - dFrecoil);
  const Double_t dFcut2(TMath::Pi() + dFrecoil);


  // hot towers
  const UInt_t hotTwrs[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  cout << "    Hot towers loaded:\n"
       << "      " << NHotTwr << " hot towers."
       << endl;


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
  UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
  UInt_t   PrimaryTrackArray_fBits[NTrkMax];
  Double_t PrimaryTrackArray_nHitsFit[NTrkMax];
  Double_t PrimaryTrackArray_nHitsPoss[NTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
  Double_t PrimaryTrackArray_pZ[NTrkMax];
  Double_t PrimaryTrackArray_pX[NTrkMax];
  Double_t PrimaryTrackArray_pY[NTrkMax];
  Double_t PrimaryTrackArray_pT[NTrkMax];
  Double_t PrimaryTrackArray_dEdx[NTrkMax];
  Double_t PrimaryTrackArray_charge[NTrkMax];
  Double_t PrimaryTrackArray_tofBeta[NTrkMax];
  Double_t PrimaryTrackArray_eta[NTrkMax];
  Double_t PrimaryTrackArray_phi[NTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_nSigPion[NTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_nSigProton[NTrkMax];
  Double_t PrimaryTrackArray_dcag[NTrkMax];
  Double_t PrimaryTrackArray_nHits[NTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[NTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[NTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[NTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[NTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[NTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[NTrkMax];
  Double_t PrimaryTrackArray_pathLength[NTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[NTrkMax];
  Int_t    TowerArray_;
  UInt_t   TowerArray_fUniqueID[NTwrMax];
  UInt_t   TowerArray_fBits[NTwrMax];
  Int_t    TowerArray_TwrId[NTwrMax];
  Float_t  TowerArray_TwrEng[NTwrMax];
  Float_t  TowerArray_TwrEta[NTwrMax];
  Float_t  TowerArray_TwrPhi[NTwrMax];
  Float_t  TowerArray_TwrADC[NTwrMax];
  Float_t  TowerArray_TwrPed[NTwrMax];
  Float_t  TowerArray_TwrRMS[NTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[NTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[NTwrMax];
  Float_t  TowerArray_TwrMatchP[NTwrMax];
  Float_t  TowerArray_TwrPx[NTwrMax];
  Float_t  TowerArray_TwrPy[NTwrMax];
  Float_t  TowerArray_TwrPz[NTwrMax];
  Int_t    TowerArray_fNAssocTracks[NTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_P[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigPi[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigK[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigP[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigE[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_dcag[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_eta[NTwrMax][NMatMax];
  Float_t  TowerArray_fMatchedTracksArray_pT[NTwrMax][NMatMax];
  Int_t    TowerArray_fMatchedTracksArray_nFit[NTwrMax][NMatMax];
  Int_t    TowerArray_fMatchedTracksArray_nPos[NTwrMax][NMatMax];

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


  // declare histograms
  const UInt_t  nNum(200);
  const UInt_t  nEta(80);
  const UInt_t  nPhi(120);
  const UInt_t  nPt(200);
  const UInt_t  nDh(160);
  const UInt_t  nDf(120);
  const Float_t num[2] = {0., 200.};
  const Float_t eta[2] = {-2., 2.};
  const Float_t phi[2] = {-3.15, 3.15};
  const Float_t pT[2]  = {0., 100.};
  const Float_t dH[2]  = {-4., 4.};
  const Float_t dF[2]  = {-1.5, 4.8};
  // trigger / event histograms
  TH1D *hTrkNum = new TH1D("hTrkNum", "", nNum, num[0], num[1]);
  TH1D *hTwrNum = new TH1D("hTwrNum", "", nNum, num[0], num[1]);
  TH1D *hTrgEta = new TH1D("hTrgEta", "", nEta, eta[0], eta[1]);
  TH1D *hTrgPhi = new TH1D("hTrgPhi", "", nPhi, phi[0], phi[1]);
  TH1D *hTrgEt  = new TH1D("hTrgEt", "", nPt, pT[0], pT[1]);
  // track histograms
  TH1D *hTrkEta = new TH1D("hTrkEta", "", nEta, eta[0], eta[1]);
  TH1D *hTrkPhi = new TH1D("hTrkPhi", "", nPhi, phi[0], phi[1]);
  TH1D *hTrkPt  = new TH1D("hTrkPt", "", nPt, pT[0], pT[1]);
  TH1D *hTrkDh  = new TH1D("hTrkDh", "", nDh, dH[0], dH[1]);
  TH1D *hTrkDf  = new TH1D("hTrkDf", "", nDf, dF[0], dF[1]);
  // tower histograms
  TH1D *hTwrEta = new TH1D("hTwrEta", "", nEta, eta[0], eta[1]);
  TH1D *hTwrPhi = new TH1D("hTwrPhi", "", nPhi, phi[0], phi[1]);
  TH1D *hTwrPt  = new TH1D("hTwrPt", "", nPt, pT[0], pT[1]);
  TH1D *hTwrDh  = new TH1D("hTwrDh", "", nDh, dH[0], dH[1]);
  TH1D *hTwrDf  = new TH1D("hTwrDf", "", nDf, dF[0], dF[1]);
  // errors
  hTrkNum -> Sumw2();
  hTwrNum -> Sumw2();
  hTrgEta -> Sumw2();
  hTrgPhi -> Sumw2();
  hTrgEt  -> Sumw2();
  hTrkEta -> Sumw2();
  hTrkPhi -> Sumw2();
  hTrkPt  -> Sumw2();
  hTrkDh  -> Sumw2();
  hTrkDf  -> Sumw2();
  hTwrEta -> Sumw2();
  hTwrPhi -> Sumw2();
  hTwrPt  -> Sumw2();
  hTwrDh  -> Sumw2();
  hTwrDf  -> Sumw2();
  cout << "    Declared histograms." << endl;


  UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " evts. to process." << endl;

  // event loop
  Int_t  bytes(0);
  UInt_t nTrgs(0);
  UInt_t nBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    bytes = tInput -> GetEntry(iEvt);
    if (bytes < 0) {
      cerr << "WARNING: something weird at event " << iEvt + 1 << "!" << endl;
      break;
    }
    else {
      if (!isInBatchMode) {
        cout << "      processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
        if ((iEvt + 1) == nEvts) cout << endl;
      }
      else {
        cout << "      processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      }
    }


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
    const Double_t hClust = EClustetav1;
    const Double_t fTrg   = EClustphiv1;
    const Double_t eTrg   = EClustEneT0;
    const Double_t eTtrg  = eTrg * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hClust)));

    Bool_t isHotTrg(false);
    for (UInt_t iHot = 0; iHot < NHotTwr; iHot++) {
      if (trgID == hotTwrs[iHot]) {
        isHotTrg = true;
        break;
      }
    }

    // trigger cuts
    const Bool_t isInAdcCut    = (adc <= adcMax);
    const Bool_t isInStrpCut   = ((eHstrp >= eHmin) && (eFstrp >= eFmin));
    const Bool_t isInProjCut   = (pProj < pProjMax);
    const Bool_t isInTspTrgCut = ((tsp > tspMin) && (tsp < tspMax));
    const Bool_t isInEtaTrgCut = (TMath::Abs(hTrg) < hTrgMax);
    const Bool_t isInEtTrgCut  = ((eTtrg > eTtrgMin) && (eTtrg < eTtrgMax));
    if (isHotTrg || !isInAdcCut || !isInStrpCut || !isInProjCut || !isInTspTrgCut || !isInEtaTrgCut || !isInEtTrgCut) continue;


    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      const Double_t nFit  = PrimaryTrackArray_nHitsFit[iTrk];
      const Double_t nPos  = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFit  = nFit / nPos;
      const Double_t dca   = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk  = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk  = PrimaryTrackArray_phi[iTrk];
      const Double_t pTtrk = PrimaryTrackArray_pT[iTrk];

      Double_t dHtrk = hTrk - hClust;
      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();


      // track cuts
      const Bool_t isInFitCut    = (nFit > nFitMin);
      const Bool_t isInRatioCut  = (rFit > rFitMin);
      const Bool_t isInDcaCut    = (dca < dcaMax);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtTrkCut  = (pTtrk > pTtrkMin);
      if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut || !isInPtTrkCut) continue;

      // recoil condition
      const Bool_t isRecoil = ((dFtrk > dFcut1) && (dFtrk < dFcut2));
      if (!isRecoil) continue;

      // fill track histograms
      hTrkEta -> Fill(hTrk);
      hTrkPhi -> Fill(fTrk);
      hTrkPt  -> Fill(pTtrk);
      hTrkDh  -> Fill(dHtrk);
      hTrkDf  -> Fill(dFtrk);

    }  // end track loop


    // tower loop
    for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

      const Int_t    twrID = TowerArray_TwrId[iTwr];
      const Double_t hTwr  = TowerArray_TwrEta[iTwr];
      const Double_t fTwr  = TowerArray_TwrPhi[iTwr];
      const Double_t eTwr  = TowerArray_TwrEng[iTwr];
      const Double_t pXtwr = TowerArray_TwrPx[iTwr];
      const Double_t pYtwr = TowerArray_TwrPy[iTwr];
      const Double_t pTtwr = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr));

      Double_t dHtwr = hTwr - hClust;
      Double_t dFtwr = fTwr - fTrg;
      if (dFtwr < (-1. * TMath::PiOver2())) dFtwr += TMath::TwoPi();
      if (dFtwr > (3. * TMath::PiOver2()))  dFtwr -= TMath::TwoPi();

      Bool_t isHotTwr(false);
      for (UInt_t iHot = 0; iHot < NHotTwr; iHot++) {
        if (twrID == hotTwrs[iHot]) {
          isHotTwr = true;
          break;
        }
      }


      // tower cuts
      const Bool_t isNotTrg      = (twrID != trgID);
      const Bool_t isInEneCut    = (eTwr > eTwrMin);
      const Bool_t isInPtTwrCut  = (pTtwr > pTtwrMin);
      const Bool_t isInEtaTwrCut = (TMath::Abs(hTwr) < hTwrMax);
      if (isHotTwr || !isNotTrg || !isInEneCut || !isInEtaTwrCut) continue;

      // recoil condition
      const Bool_t isRecoil = ((dFtwr > dFcut1) && (dFtwr < dFcut2));
      if (!isRecoil) continue;

      // fill tower histograms
      hTwrEta -> Fill(hTwr);
      hTwrPhi -> Fill(fTwr);
      hTwrPt  -> Fill(pTtwr);
      hTwrDh  -> Fill(dHtwr);
      hTwrDf  -> Fill(dFtwr);

    }  // end tower loop


    // fill trigger histograms
    hTrkNum -> Fill(nTrks);
    hTwrNum -> Fill(nTwrs);
    hTrgEta -> Fill(hClust);
    hTrgPhi -> Fill(fTrg);
    hTrgEt  -> Fill(eTtrg);
    nTrgs++;

  }  // end event loop

  cout << "    Event loop finished: found " << nTrgs << " triggers." << endl;


  // save histograms
  fOutput -> cd();
  hTrkNum -> Write();
  hTwrNum -> Write();
  hTrgEta -> Write();
  hTrgPhi -> Write();
  hTrgEt  -> Write();
  hTrkEta -> Write();
  hTrkPhi -> Write();
  hTrkPt  -> Write();
  hTrkDh  -> Write();
  hTrkDf  -> Write();
  hTwrEta -> Write();
  hTwrPhi -> Write();
  hTwrPt  -> Write();
  hTwrDh  -> Write();
  hTwrDf  -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
