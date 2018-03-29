// 'EmbeddingTowerStudy.C'
// Derek Anderson
// 07.11.2017
//
// Use this to extract the tower distributions
// from the embedding 'femtoDst' tree.


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
static const UInt_t nEbins(5);
static const UInt_t nPtBins(5);
static const UInt_t nHotTwrs(3);
static const UInt_t nTrkMax(1000);
static const UInt_t nTwrMax(10000);
static const UInt_t nMatMax(100);
static const UInt_t fTrgMode(3);
// io parameters
static const TString sTree("GfmtoDst_mu");
static const TString sInDefault("../../Embedding/Run12pp/MuDstMatching/output/merged/pt5_-1.ePcalc2.root");
static const TString sOutDefault("pp200r12.highTwr.d3m9y2017.root");
static const TString sHotTwrList("text/HotTowerCheck.pp200r12.d15m8y2017.list");
static const TString sTrkDirs[nDist] = {"AllTracks", "RecoilTracks"};
static const TString sTwrDirs[nDist] = {"AllTowers", "RecoilTowers"};



void EmbeddingTowerStudy(const Bool_t isInBatchMode=false, const TString sInput=sInDefault, const TString sOutput=sOutDefault) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning tower script..." << endl;

  // constants
  const Bool_t   isInTrgMode0  = (fTrgMode == 0);
  const Bool_t   isInTrgMode1  = (fTrgMode == 1);
  const Bool_t   isInTrgMode2  = (fTrgMode == 2);
  const Bool_t   isInTrgMode3  = (fTrgMode == 3);
  const Bool_t   removeHotTwrs = true;
  const Double_t fRecoilMax    = TMath::PiOver4();
  const Double_t mPion         = 0.140;
  const Double_t dFmax         = 0.025;
  const Double_t dHmax         = 0.025;
  const Double_t hTrgMax       = 0.9;
  const Double_t eTtrgMin      = 9.0;
  const Double_t eTtrgMax      = 20.;
  const Double_t nFitMin       = 15;
  const Double_t rFitMin       = 0.52;
  const Double_t dcaMax        = 1.0;
  const Double_t hTrkMax       = 1.0;
  const Double_t pTtrkMin      = 0.2;
  const Double_t pTtrgMin      = 9.0;
  const Double_t hTwrMax       = 1.0;
  const Double_t eTwrMin       = 0.2;
  const Double_t eCorrMin      = 0.2;
  const Double_t pidCut        = 3.;
  const Double_t pidMax        = 999.;
  // hot-tower / E/p constants
  const Double_t eCut[nEbins]      = {0., 1., 5., 10., 25.};
  const Double_t pTmatMin[nPtBins] = {0., 0., 2., 5., 10.};
  const Double_t pTmatMax[nPtBins] = {100., 2., 5., 10., 100.};


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


  // declarare leaf types
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
  Int_t    PrimaryTrackArray_nHitsFit[nTrkMax];
  Int_t    PrimaryTrackArray_nHitsPoss[nTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[nTrkMax];
  Int_t    PrimaryTrackArray_pdgId[nTrkMax];
  Int_t    PrimaryTrackArray_geantId[nTrkMax];
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
  tInput -> SetBranchAddress("PrimaryTrackArray.pdgId", PrimaryTrackArray_pdgId);
  tInput -> SetBranchAddress("PrimaryTrackArray.geantId", PrimaryTrackArray_geantId);
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
  ifstream hotTwrList(sHotTwrList.Data());
  if (!hotTwrList && removeHotTwrs) {
    cerr << "PANIC: couldn't open hot tower list!" << endl;
    return;
  }

  UInt_t iHotTwr(0);
  UInt_t hotTwrID(0);
  UInt_t hotTwrs[nHotTwrs];
  if (removeHotTwrs) {
    while (hotTwrList) {
      hotTwrList >> hotTwrID;
      hotTwrs[iHotTwr] = hotTwrID;
      iHotTwr++;
      if (iHotTwr == nHotTwrs) break;
    }
    hotTwrList.close();
    cout << "    Hot towers loaded:\n"
         << "      " << iHotTwr << " towers."
         << endl;
  }


  fOutput -> cd();
  // trigger histograms
  TH1D *hNumTrg;
  TH1D *hTrgEt;
  TH1D *hTrgEne;
  TH1D *hTrgEta;
  // E/p histograms
  TH1D *hEpTwrEne;
  TH1D *hEpTwrEneL5;
  TH1D *hEpTrkPtE[nPtBins];
  TH1D *hEpTrkPtH[nPtBins];
  TH1D *hEpElec[nPtBins];
  TH1D *hEpHadr[nPtBins];
  TH1D *hEpTrkPtEL5[nPtBins];
  TH1D *hEpTrkPtHL5[nPtBins];
  TH1D *hEpElecL5[nPtBins];
  TH1D *hEpHadrL5[nPtBins];
  // track (triggered) histograms
  TH1D *hTrkEne[nDist];
  TH1D *hTrkPt[nDist];
  TH1D *hTrkEta[nDist];
  TH1D *hTrkPhi[nDist];
  TH1D *hTrkDeta[nDist];
  TH1D *hTrkDphi[nDist];
  // tower (triggered) histograms
  TH1D *hTwrEne[nDist];
  TH1D *hTwrCorr[nDist];
  TH1D *hTwrPt[nDist];
  TH1D *hTwrEta[nDist];
  TH1D *hTwrPhi[nDist];
  TH1D *hTwrDeta[nDist];
  TH1D *hTwrDphi[nDist];
  // hot tower histograms
  TH1D *hHotTwrID[nEbins];
  TH1D *hHotTwrEne[nEbins];
  TH1D *hHotTwrEta[nEbins];
  TH2D *hHotTwrIDvsEne[nEbins];
  TH2D *hHotTwrIDvsEta[nEbins];

  const UInt_t   nN = 200;
  const UInt_t   nE = 1050;
  const UInt_t   nH = 40;
  const UInt_t   nF = 720;
  const UInt_t   nP = 200;
  const UInt_t   nT = 4800;
  const Double_t n1 = 0.;
  const Double_t n2 = 200.;
  const Double_t e1 = -5.;
  const Double_t e2 = 100.;
  const Double_t h1 = -2.;
  const Double_t h2 = 2.;
  const Double_t f1 = -1. * TMath::TwoPi();
  const Double_t f2 = TMath::TwoPi();
  const Double_t p1 = 0.;
  const Double_t p2 = 20.;
  const Double_t t1 = 0.;
  const Double_t t2 = 4800.;
  // trigger histograms
  hNumTrg           = new TH1D("hNumTrg", "No. of triggers", nN, n1, n2);
  hTrgEt            = new TH1D("hTrgEt", "Trigger E_{T}", nE, e1, e2);
  hTrgEne           = new TH1D("hTrgEne", "Trigger energy", nE, e1, e2);
  hTrgEta           = new TH1D("hTrgEta", "Trigger #eta", nH, h1, h2);
  // E/p histograms
  hEpTwrEne         = new TH1D("hEpTwrEne", "Tower energy, trigger excluded and has one match (for E/p)", nE, e1, e2);
  hEpTwrEneL5       = new TH1D("hEpTwrEneL5", "Tower energy for E_{twr} < 5, trigger excluded and has one match (for E/p)", nE, e1, e2);
  hEpTrkPtE[0]      = new TH1D("hEpTrkPtE0", "Matched e^{#pm} p_{T}, (0, 100) GeV/c", nE, e1, e2);
  hEpTrkPtE[1]      = new TH1D("hEpTrkPtE1", "Matched e^{#pm} p_{T}, (0, 2) GeV/c", nE, e1, e2);
  hEpTrkPtE[2]      = new TH1D("hEpTrkPtE2", "Matched e^{#pm} p_{T}, (2, 5) GeV/c", nE, e1, e2);
  hEpTrkPtE[3]      = new TH1D("hEpTrkPtE3", "Matched e^{#pm} p_{T}, (5, 10) GeV/c", nE, e1, e2);
  hEpTrkPtE[4]      = new TH1D("hEpTrkPtE4", "Matched e^{#pm} p_{T}, (10, 100) GeV/c", nE, e1, e2);
  hEpTrkPtH[0]      = new TH1D("hEpTrkPtH0", "Matched h^{#pm} p_{T}, (0, 100) GeV/c", nE, e1, e2);
  hEpTrkPtH[1]      = new TH1D("hEpTrkPtH1", "Matched h^{#pm} p_{T}, (0, 2) GeV/c", nE, e1, e2);
  hEpTrkPtH[2]      = new TH1D("hEpTrkPtH2", "Matched h^{#pm} p_{T}, (2, 5) GeV/c", nE, e1, e2);
  hEpTrkPtH[3]      = new TH1D("hEpTrkPtH3", "Matched h^{#pm} p_{T}, (5, 10) GeV/c", nE, e1, e2);
  hEpTrkPtH[4]      = new TH1D("hEpTrkPtH4", "Matched h^{#pm} p_{T}, (10, 100) GeV/c", nE, e1, e2);
  hEpElec[0]        = new TH1D("hEpElec0", "E/p e^{#pm}, p_{T} #in (0, 100) GeV/c", nP, p1, p2);
  hEpElec[1]        = new TH1D("hEpElec1", "E/p e^{#pm}, p_{T} #in (0, 2) GeV/c", nP, p1, p2);
  hEpElec[2]        = new TH1D("hEpElec2", "E/p e^{#pm}, p_{T} #in (2, 5) GeV/c", nP, p1, p2);
  hEpElec[3]        = new TH1D("hEpElec3", "E/p e^{#pm}, p_{T} #in (5, 10) GeV/c", nP, p1, p2);
  hEpElec[4]        = new TH1D("hEpElec4", "E/p e^{#pm}, p_{T} #in (10, 100) GeV/c", nP, p1, p2);
  hEpHadr[0]        = new TH1D("hEpHadr0", "E/p h^{#pm}, p_{T} #in (0, 100) GeV/c", nP, p1, p2);
  hEpHadr[1]        = new TH1D("hEpHadr1", "E/p h^{#pm}, p_{T} #in (0, 2) GeV/c", nP, p1, p2);
  hEpHadr[2]        = new TH1D("hEpHadr2", "E/p h^{#pm}, p_{T} #in (2, 5) GeV/c", nP, p1, p2);
  hEpHadr[3]        = new TH1D("hEpHadr3", "E/p h^{#pm}, p_{T} #in (5, 10) GeV/c", nP, p1, p2);
  hEpHadr[4]        = new TH1D("hEpHadr4", "E/p h^{#pm}, p_{T} #in (10, 100) GeV/c", nP, p1, p2);
  hEpTrkPtEL5[0]    = new TH1D("hEpTrkPtE0L5", "Matched e^{#pm} p_{T}, (0, 100) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtEL5[1]    = new TH1D("hEpTrkPtE1L5", "Matched e^{#pm} p_{T}, (0, 2) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtEL5[2]    = new TH1D("hEpTrkPtE2L5", "Matched e^{#pm} p_{T}, (2, 5) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtEL5[3]    = new TH1D("hEpTrkPtE3L5", "Matched e^{#pm} p_{T}, (5, 10) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtEL5[4]    = new TH1D("hEpTrkPtE4L5", "Matched e^{#pm} p_{T}, (10, 100) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtHL5[0]    = new TH1D("hEpTrkPtH0L5", "Matched h^{#pm} p_{T}, (0, 100) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtHL5[1]    = new TH1D("hEpTrkPtH1L5", "Matched h^{#pm} p_{T}, (0, 2) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtHL5[2]    = new TH1D("hEpTrkPtH2L5", "Matched h^{#pm} p_{T}, (2, 5) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtHL5[3]    = new TH1D("hEpTrkPtH3L5", "Matched h^{#pm} p_{T}, (5, 10) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpTrkPtHL5[4]    = new TH1D("hEpTrkPtH4L5", "Matched h^{#pm} p_{T}, (10, 100) GeV/c and E_{twr} < 5 GeV", nE, e1, e2);
  hEpElecL5[0]      = new TH1D("hEpElec0L5", "E/p e^{#pm}, p_{T} #in (0, 100) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpElecL5[1]      = new TH1D("hEpElec1L5", "E/p e^{#pm}, p_{T} #in (0, 2) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpElecL5[2]      = new TH1D("hEpElec2L5", "E/p e^{#pm}, p_{T} #in (2, 5) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpElecL5[3]      = new TH1D("hEpElec3L5", "E/p e^{#pm}, p_{T} #in (5, 10) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpElecL5[4]      = new TH1D("hEpElec4L5", "E/p e^{#pm}, p_{T} #in (10, 100) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpHadrL5[0]      = new TH1D("hEpHadr0L5", "E/p h^{#pm}, p_{T} #in (0, 100) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpHadrL5[1]      = new TH1D("hEpHadr1L5", "E/p h^{#pm}, p_{T} #in (0, 2) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpHadrL5[2]      = new TH1D("hEpHadr2L5", "E/p h^{#pm}, p_{T} #in (2, 5) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpHadrL5[3]      = new TH1D("hEpHadr3L5", "E/p h^{#pm}, p_{T} #in (5, 10) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  hEpHadrL5[4]      = new TH1D("hEpHadr4L5", "E/p h^{#pm}, p_{T} #in (10, 100) GeV/c and E_{twr} < 5 GeV", nP, p1, p2);
  // triggered track histograms
  hTrkEne[0]        = new TH1D("hTrkEneAll", "Track energy, all", nE, e1, e2);
  hTrkEne[1]        = new TH1D("hTrkEneRec", "Track energy, recoil", nE, e1, e2);
  hTrkPt[0]         = new TH1D("hTrkPtAll", "Track p_{T}, all", nE, e1, e2);
  hTrkPt[1]         = new TH1D("hTrkPtRec", "Track p_{T}, recoil", nE, e1, e2);
  hTrkEta[0]        = new TH1D("hTrkEtaAll", "Track #eta, all", nH, h1, h2);
  hTrkEta[1]        = new TH1D("hTrkEtaRec", "Track #eta, recoil", nH, h1, h2);
  hTrkPhi[0]        = new TH1D("hTrkPhiAll", "Track #varphi, all", nF, f1, f2);
  hTrkPhi[1]        = new TH1D("hTrkPhiRec", "Track #varphi, recoil", nF, f1, f2);
  hTrkDeta[0]       = new TH1D("hTrkDetaAll", "Track #Delta#eta, all", nH, h1, h2);
  hTrkDeta[1]       = new TH1D("hTrkDetaRec", "Track #Delta#eta, recoil", nH, h1, h2);
  hTrkDphi[0]       = new TH1D("hTrkDphiAll", "Track #Delta#varphi, all", nF, f1, f2);
  hTrkDphi[1]       = new TH1D("hTrkDphiRec", "Track #Delta#varphi, recoil", nF, f1, f2);
  // triggered tower histograms
  hTwrEne[0]        = new TH1D("hTwrEneAll", "Tower energy (raw), all", nE, e1, e2);
  hTwrEne[1]        = new TH1D("hTwrEneRec", "Tower energy (raw), recoil", nE, e1, e2);
  hTwrCorr[0]       = new TH1D("hTwrCorrAll", "Tower energy (corr.), all", nE, e1, e2);
  hTwrCorr[1]       = new TH1D("hTwrCorrRec", "Tower energy (corr.), recoil", nE, e1, e2);
  hTwrPt[0]         = new TH1D("hTwrPtAll", "Tower p_{T}, all", nE, e1, e2);
  hTwrPt[1]         = new TH1D("hTwrPtRec", "Tower p_{T}, recoil", nE, e1, e2);
  hTwrEta[0]        = new TH1D("hTwrEtaAll", "Tower #eta, all", nH, h1, h2);
  hTwrEta[1]        = new TH1D("hTwrEtaRec", "Tower #eta, recoil", nH, h1, h2);
  hTwrPhi[0]        = new TH1D("hTwrPhiAll", "Tower #varphi, all", nF, f1, f2);
  hTwrPhi[1]        = new TH1D("hTwrPhiRec", "Tower #varphi, recoil", nF, f1, f2);
  hTwrDeta[0]       = new TH1D("hTwrDetaAll", "Tower #Delta#eta, all", nH, h1, h2);
  hTwrDeta[1]       = new TH1D("hTwrDetaRec", "Tower #Delta#eta, recoil", nH, h1, h2);
  hTwrDphi[0]       = new TH1D("hTwrDphiAll", "Tower #Delta#varphi, all", nF, f1, f2);
  hTwrDphi[1]       = new TH1D("hTwrDphiRec", "Tower #Delta#varphi, recoil", nF, f1, f2);
  // hot tower histograms
  hHotTwrID[0]      = new TH1D("hHotTwrID_e0", "Tower ID (for hot towers), all", nT, t1, t2);
  hHotTwrID[1]      = new TH1D("hHotTwrID_e1", "Tower ID (for hot towers), E > 1", nT, t1, t2);
  hHotTwrID[2]      = new TH1D("hHotTwrID_e5", "Tower ID (for hot towers), E > 5", nT ,t1, t2);
  hHotTwrID[3]      = new TH1D("hHotTwrID_e10", "Tower ID (for hot towers), E > 10", nT, t1, t2);
  hHotTwrID[4]      = new TH1D("hHotTwrID_e25", "Tower ID (for hot towers), E > 25", nT, t1, t2);
  hHotTwrEne[0]     = new TH1D("hHotTwrEne_e0", "Tower energy (for hot towers), all", nE, e1, e2);
  hHotTwrEne[1]     = new TH1D("hHotTwrEne_e1", "Tower energy (for hot towers), E > 1", nE, e1, e2);
  hHotTwrEne[2]     = new TH1D("hHotTwrEne_e5", "Tower energy (for hot towers), E > 5", nE, e1, e2);
  hHotTwrEne[3]     = new TH1D("hHotTwrEne_e10", "Tower energy (for hot towers), E > 10", nE, e1, e2);
  hHotTwrEne[4]     = new TH1D("hHotTwrEne_e25", "TowerEnergy (for hot towers), E > 25", nE, e1, e2);
  hHotTwrEta[0]     = new TH1D("hHotTwrEta_e0", "Tower #eta (for hot towers), all", nH, h1, h2);
  hHotTwrEta[1]     = new TH1D("hHotTwrEta_e1", "Tower #eta (for hot towers), E > 1", nH, h1, h2);
  hHotTwrEta[2]     = new TH1D("hHotTwrEta_e5", "Tower #eta (for hot towers), E > 5", nH, h1, h2);
  hHotTwrEta[3]     = new TH1D("hHotTwrEta_e10", "Tower #eta (for hot towers), E > 10", nH, h1, h2);
  hHotTwrEta[4]     = new TH1D("hHotTwrEta_e25", "TowerEnergy (for hot towers), E > 25", nH, h1, h2);
  hHotTwrIDvsEne[0] = new TH2D("hTwrIDvsEne_e0", "Tower ID vs. energy (for hot towers), all", nE, e1, e2, nT, t1, t2);
  hHotTwrIDvsEne[1] = new TH2D("hTwrIDvsEne_e1", "Tower ID vs. energy (for hot towers), E > 1", nE, e1, e2, nT, t1, t2);
  hHotTwrIDvsEne[2] = new TH2D("hTwrIDvsEne_e5", "Tower ID vs. energy (for hot towers), E > 5", nE ,e1, e2, nT, t1, t2);
  hHotTwrIDvsEne[3] = new TH2D("hTwrIDvsEne_e10", "Tower ID vs. energy (for hot towers), E > 10", nE, e1, e2, nT, t1, t2);
  hHotTwrIDvsEne[4] = new TH2D("hTwrIDvsEne_e25", "Tower ID vs. energy (for hot towers), E > 25", nE, e1, e2, nT, t1, t2);
  hHotTwrIDvsEta[0] = new TH2D("hTwrIDvsEta_e0", "Tower ID vs. #eta (for hot towers), all", nH, h1, h2, nT, t1, t2);
  hHotTwrIDvsEta[1] = new TH2D("hTwrIDvsEta_e1", "Tower ID vs. #eta (for hot towers), E > 1", nH, h1, h2, nT, t1, t2);
  hHotTwrIDvsEta[2] = new TH2D("hTwrIDvsEta_e5", "Tower ID vs. #eta (for hot towers), E > 5", nH ,h1, h2, nT, t1, t2);
  hHotTwrIDvsEta[3] = new TH2D("hTwrIDvsEta_e10", "Tower ID vs. #eta (for hot towers), E > 10", nH, h1, h2, nT, t1, t2);
  hHotTwrIDvsEta[4] = new TH2D("hTwrIDvsEta_e25", "Tower ID vs. #eta (for hot towers), E > 25", nH, h1, h2, nT, t1, t2);
  // errors
  hTrgEt      -> Sumw2();
  hTrgEne     -> Sumw2();
  hTrgEta     -> Sumw2();
  hEpTwrEne   -> Sumw2();
  hEpTwrEneL5 -> Sumw2();
  for (Long64_t iPtBin = 0; iPtBin < nPtBins; iPtBin++) {
    hEpTrkPtE[iPtBin]   -> Sumw2();
    hEpTrkPtH[iPtBin]   -> Sumw2();
    hEpElec[iPtBin]     -> Sumw2();
    hEpHadr[iPtBin]     -> Sumw2();
    hEpTrkPtEL5[iPtBin] -> Sumw2();
    hEpTrkPtHL5[iPtBin] -> Sumw2();
    hEpElecL5[iPtBin]   -> Sumw2();
    hEpHadrL5[iPtBin]   -> Sumw2();
  }
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    hTrkEne[iDist]  -> Sumw2();
    hTrkPt[iDist]   -> Sumw2();
    hTrkEta[iDist]  -> Sumw2();
    hTrkPhi[iDist]  -> Sumw2();
    hTrkDeta[iDist] -> Sumw2();
    hTrkDphi[iDist] -> Sumw2();
    hTwrEne[iDist]  -> Sumw2();
    hTwrCorr[iDist] -> Sumw2();
    hTwrPt[iDist]   -> Sumw2();
    hTwrEta[iDist]  -> Sumw2();
    hTwrPhi[iDist]  -> Sumw2();
    hTwrDeta[iDist] -> Sumw2();
    hTwrDphi[iDist] -> Sumw2();
  }
  for (Long64_t iEbin = 0; iEbin < nEbins; iEbin++) {
    hHotTwrEne[iEbin]     -> Sumw2();
    hHotTwrEta[iEbin]     -> Sumw2();
    hHotTwrIDvsEne[iEbin] -> Sumw2();
    hHotTwrIDvsEta[iEbin] -> Sumw2();
  }


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
    else
      cout << "      processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;


    if (isInTrgMode0) {
      iTrig = -1;
      fTrig = 0.;
      hTrig = 0.;
      nTrg++;
    }


    // trigger info
    const Int_t    trgID  = ETwrdidT;
    const Int_t    nTrks  = nPrimaryTracks;
    const Int_t    nTwrs  = TowerArray_;
    const Double_t hTrg   = ETwreT;
    const Double_t fClust = EClustphiv1;
    const Double_t hClust = EClustetav1;
    const Double_t eTrg   = EClustEneT0;
    const Double_t eTtrg  = eTrg * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hClust)));

    // trigger cuts
    if (isInTrgMode1) {
      const Bool_t isInEtaTrgCut = (TMath::Abs(hTrg) < hTrgMax);
      const Bool_t isInEneTrgCut = (eTtrg > eTtrgMin);
      const Bool_t foundTrg1     = (isInEtaTrgCut && isInEneTrgCut);
      if (!foundTrg1)
        continue;
      else {
        iTrig = trgID;
        fTrig = fClust;
        hTrig = hClust;
        nTrg++;
      }
    }


    // track loop 1
    if (isInTrgMode2) {
      Bool_t foundTrg2 = false;
      for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

        // track info
        const Int_t    nFit  = PrimaryTrackArray_nHitsFit[iTrk];
        const Int_t    nPos  = PrimaryTrackArray_nHitsPoss[iTrk];
        const Double_t rFit  = (Double_t) nFit / (Double_t) nPos;
        const Double_t dca   = PrimaryTrackArray_dcag[iTrk];
        const Double_t hTrk  = PrimaryTrackArray_eta[iTrk];
        const Double_t fTrk  = PrimaryTrackArray_phi[iTrk];
        const Double_t pTtrk = PrimaryTrackArray_pT[iTrk];

        // track cuts
        const Bool_t isInNFitTrkCut = (nFit > nFitMin);
        const Bool_t isInRFitTrkCut = (rFit > rFitMin);
        const Bool_t isInDcaTrkCut  = (dca < dcaMax);
        const Bool_t isInEtaTrkCut  = (TMath::Abs(hTrk) < hTrkMax);
        const Bool_t isInPtTrkCut   = (pTtrk > pTtrkMin);
        if (!isInNFitTrkCut || !isInRFitTrkCut || !isInDcaTrkCut || !isInEtaTrkCut || !isInPtTrkCut) continue;

        // check if trigger and match to tower
        const Bool_t isInPtTrgCut = (pTtrk > pTtrgMin);
        if (isInPtTrgCut) {
          Bool_t isMatched = false;
          for (Long64_t iTwr = 0; iTwr < nTwrs; iTwr++) {
            const Double_t fTwr    = TowerArray_TwrPhi[iTwr];
            const Double_t hTwr    = TowerArray_TwrEta[iTwr];
            const Double_t dPhi    = fTrk - fTwr;
            const Double_t dEta    = hTrk - hTwr;
            const Bool_t   isInPhi = (dPhi < dFmax);
            const Bool_t   isInEta = (dEta < dHmax);
            if (isInPhi && isInEta) {
              iTrig     = TowerArray_TwrId[iTwr];
              fTrig     = fTwr;
              hTrig     = hTwr;
              isMatched = true;
              break;
            }
          }
          if (isMatched) {
            foundTrg2 = true;
            break;
          }
        }

      }  // end track loop 1
      if (!foundTrg2)
        continue;
      else
        nTrg++;
    }


    // tower loop 1
    UInt_t numTrg = 0;
    if (isInTrgMode3) {
      Bool_t foundTrg3 = false;
      for (Long64_t iTwr = 0; iTwr <  nTwrs; iTwr++) {

        // tower info
        const Int_t    twrID  = TowerArray_TwrId[iTwr];
        const UInt_t   fMatch = TowerArray_TwrMatchIdnex[iTwr];
        const UInt_t   nMatch = TowerArray_NoOfmatchedTrk[iTwr];
        const Double_t fTwr   = TowerArray_TwrPhi[iTwr];
        const Double_t hTwr   = TowerArray_TwrEta[iTwr];
        const Double_t eTwr   = TowerArray_TwrEng[iTwr];
        const Double_t eTtwr  = eTwr * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTwr)));
        const Double_t pXtwr  = TowerArray_TwrPx[iTwr];
        const Double_t pYtwr  = TowerArray_TwrPy[iTwr];
        const Double_t pTtwr  = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr));
        // for E/p
        const Bool_t isInEtaTwrCut = (TMath::Abs(hTwr) < hTwrMax);
        const Bool_t isInEneTwrCut = (eTwr > eTwrMin);
        const Bool_t isLessThan5   = (eTwr < 5);
        const Bool_t hasOneMatch   = (nMatch == 1);


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

        for (Long64_t iEbin = 0; iEbin < nEbins; iEbin++) {
          if (eTwr > eCut[iEbin]) {
            hHotTwrID[iEbin]      -> Fill(twrID);
            hHotTwrEne[iEbin]     -> Fill(eTwr);
            hHotTwrEta[iEbin]     -> Fill(hTwr);
            hHotTwrIDvsEne[iEbin] -> Fill(eTwr, twrID);
            hHotTwrIDvsEta[iEbin] -> Fill(hTwr, twrID);
          }
        }


        // calculate corrected energy and E/p
        if (isInEtaTwrCut && isInEneTwrCut && hasOneMatch) {
          hEpTwrEne -> Fill(eTwr);
          if (isLessThan5) hEpTwrEneL5 -> Fill(eTwr);
        }

        Int_t    nFitMat = 0;
        Int_t    nPosMat = 0;
        Float_t  pMat    = 0.;
        Float_t  nPiMat  = 0.;
        Float_t  nKmat   = 0.;
        Float_t  nPmat   = 0.;
        Float_t  nEmat   = 0.;
        Float_t  dcaMat  = 0.;
        Float_t  hMat    = 0.;
        Float_t  pTmat   = 0.;
        Double_t rFitMat = 0.;
        Double_t eSum    = 0.;
        Double_t eCorr   = 0.;

        const Bool_t isMatched    = (fMatch == 1);
        const Bool_t isNotMatched = (fMatch == 0);
        const Bool_t hasMatches   = (nMatch >= 1);
        const Bool_t hasNoMatches = (nMatch == 0);
        if (isMatched && hasMatches) {
          for (Long64_t iMatch = 0; iMatch < nMatch; iMatch++) {
            nPiMat  = TMath::Abs(TowerArray_fMatchedTracksArray_nSigPi[iTwr][iMatch]);
            nKmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigK[iTwr][iMatch]);
            nPmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigP[iTwr][iMatch]);
            nEmat   = TMath::Abs(TowerArray_fMatchedTracksArray_nSigE[iTwr][iMatch]);
            nFitMat = TowerArray_fMatchedTracksArray_nFit[iTwr][iMatch];
            nPosMat = TowerArray_fMatchedTracksArray_nPos[iTwr][iMatch];
            rFitMat = (Double_t) nFitMat / (Double_t) nPosMat;
            dcaMat  = TowerArray_fMatchedTracksArray_dcag[iTwr][iMatch];
            hMat    = TowerArray_fMatchedTracksArray_eta[iTwr][iMatch];
            pTmat   = TowerArray_fMatchedTracksArray_pT[iTwr][iMatch];
            pMat    = TowerArray_fMatchedTracksArray_P[iTwr][iMatch];
            eSum   += pMat;

            // E/p calculation
            if (isInEtaTwrCut && isInEneTwrCut && hasOneMatch) { 
              const Bool_t isInElecCut = ((nPiMat > pidCut) && (nKmat > pidCut) && (nPmat > pidCut) && (nEmat < pidCut));
              const Bool_t isInHadrCut = ((nPiMat < pidCut) && (nKmat < pidCut) && (nPmat < pidCut) && (nEmat > pidCut));
              const Bool_t isInNFitCut = (nFitMat > nFitMin);
              const Bool_t isInRFitCut = (rFitMat > rFitMin);
              const Bool_t isInDcaCut  = (dcaMat < dcaMax);
              const Bool_t isInEtaCut  = (TMath::Abs(hMat) < hTrkMax);
              const Bool_t isInPtCut   = (pTmat > pTtrkMin);
              if (!isInNFitCut || !isInRFitCut || !isInDcaCut || !isInEtaCut || !isInPtCut) continue;

              Bool_t isInPtBin;
              for (Long64_t iPtBin = 0; iPtBin < nPtBins; iPtBin++) {
                if (iPtBin == 0)
                  isInPtBin = true;
                else
                  isInPtBin = ((pTmat > pTmatMin[iPtBin]) && (pTmat < pTmatMax[iPtBin]));
                // electrons
                if (isInElecCut && isInPtBin) {
                  hEpTrkPtE[iPtBin] -> Fill(pTmat);
                  hEpElec[iPtBin]   -> Fill(eTwr / pMat);
                  if (isLessThan5) {
                    hEpTrkPtEL5[iPtBin] -> Fill(pTmat);
                    hEpElecL5[iPtBin]   -> Fill(eTwr / pMat);
                  }
                }
                // hadrons
                if (isInHadrCut && isInPtBin) {
                  hEpTrkPtH[iPtBin] -> Fill(pTmat);
                  hEpHadr[iPtBin]   -> Fill(eTwr / pMat);
                  if (isLessThan5) {
                    hEpTrkPtHL5[iPtBin] -> Fill(pTmat);
                    hEpHadrL5[iPtBin]   -> Fill(eTwr / pMat);
                  }
                }
              }  // end bin loop
            }  // tower cuts
          }  // end track loop
          eCorr = eTwr - eSum;
          if (eCorr < 0.) eCorr = 0.;
        }  // end match condition
        else if (isNotMatched && hasNoMatches) {
          eCorr = eTwr;
        }


        // tower cuts
        const Bool_t isInCorrTwrCut = (eCorr >= eCorrMin);
        if (!isInEtaTwrCut || !isInEneTwrCut/* || !isInCorrTwrCut*/) continue;

        // check if trigger
        const Bool_t isInEtTrgCut = ((eTtwr > eTtrgMin) && (eTtwr < eTtrgMax));
        if (isInEtTrgCut) {
          iTrig     = twrID;
          fTrig     = fTwr;
          hTrig     = hTwr;
          foundTrg3 = true;
          numTrg++;
          hTrgEt  -> Fill(eTtwr);
          hTrgEne -> Fill(eTwr);
          hTrgEta -> Fill(hTwr);
          break;
        }

      }  // end tower loop 1
      if (!foundTrg3)
        continue;
      else
        nTrg++;
    }
    hNumTrg -> Fill(numTrg);


    // track loop 2
    for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    nFitTrk = PrimaryTrackArray_nHitsFit[iTrk];
      const Int_t    nPosTrk = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFitTrk = (Double_t) nFitTrk / (Double_t) nPosTrk;
      const Double_t dca     = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk    = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk    = PrimaryTrackArray_phi[iTrk];
      const Double_t pZtrk   = PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk   = PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk    = TMath::Sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (mPion * mPion));

      Double_t dEta = hTrk - hTrig;
      Double_t dPhi = fTrk - fTrig;
      if (dPhi < (-1. * TMath::PiOver2())) dPhi += TMath::TwoPi();
      if (dPhi > (3. * TMath::PiOver2()))  dPhi -= TMath::TwoPi();


      // track cuts
      const Bool_t isInNFitTrkCut = (nFitTrk > nFitMin);
      const Bool_t isInRFitTrkCut = (rFitTrk > rFitMin);
      const Bool_t isInDcaTrkCut  = (dca < dcaMax);
      const Bool_t isInEtaTrkCut  = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtTrkCut   = (pTtrk > pTtrkMin);
      if (!isInNFitTrkCut || !isInRFitTrkCut || !isInDcaTrkCut || !isInEtaTrkCut || !isInPtTrkCut) continue;

      // fill histograms
      const Bool_t isRecoilTrk = (TMath::Abs(dPhi - TMath::Pi()) < fRecoilMax);
      hTrkEne[0]  -> Fill(eTrk);
      hTrkPt[0]   -> Fill(pTtrk);
      hTrkEta[0]  -> Fill(hTrk);
      hTrkPhi[0]  -> Fill(fTrk);
      hTrkDeta[0] -> Fill(dEta);
      hTrkDphi[0] -> Fill(dPhi);
      if (isRecoilTrk) {
        hTrkEne[1]  -> Fill(eTrk);
        hTrkPt[1]   -> Fill(pTtrk);
        hTrkEta[1]  -> Fill(hTrk);
        hTrkPhi[1]  -> Fill(fTrk);
        hTrkDeta[1] -> Fill(dEta);
        hTrkDphi[1] -> Fill(dPhi);
      }

    }  // end track loop 2


    // tower loop 2
    for (Long64_t iTwr = 0; iTwr <  nTwrs; iTwr++) {

      // tower info
      const Int_t    twrID  = TowerArray_TwrId[iTwr];
      const UInt_t   fMatch = TowerArray_TwrMatchIdnex[iTwr];
      const UInt_t   nMatch = TowerArray_NoOfmatchedTrk[iTwr];
      const Double_t fTwr   = TowerArray_TwrPhi[iTwr];
      const Double_t hTwr   = TowerArray_TwrEta[iTwr];
      const Double_t eTwr   = TowerArray_TwrEng[iTwr];
      const Double_t eTtwr  = eTwr * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * eTwr)));
      const Double_t pXtwr  = TowerArray_TwrPx[iTwr];
      const Double_t pYtwr  = TowerArray_TwrPy[iTwr];
      const Double_t pTtwr  = TMath::Sqrt((pXtwr * pXtwr) + (pYtwr * pYtwr));

      Double_t dEta = hTwr - hTrig;
      Double_t dPhi = fTwr - fTrig;
      if (dPhi < (-1. * TMath::PiOver2())) dPhi += TMath::TwoPi();
      if (dPhi > (3. * TMath::PiOver2()))  dPhi -= TMath::TwoPi();


      // calculate corrected energy
      Double_t eSum  = 0.;
      Double_t eCorr = 0.;

      const Bool_t isMatched    = (fMatch == 1);
      const Bool_t isNotMatched = (fMatch == 0);
      const Bool_t hasMatches   = (nMatch >= 1);
      const Bool_t hasNoMatches = (nMatch == 0);
      if (isMatched && hasMatches) {
        for (Long64_t iMatch = 0; iMatch < nMatch; iMatch++) {
          eSum += TowerArray_fMatchedTracksArray_P[iTwr][iMatch];
        }
        eCorr = eTwr - eSum;
        if (eCorr < 0.) eCorr = 0.;
      }
      else if (isNotMatched && hasNoMatches) {
        eCorr = eTwr;
      }


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
      const Bool_t isTrigger      = (twrID == iTrig);
      const Bool_t isInEtaTwrCut  = (TMath::Abs(hTwr) < hTwrMax);
      const Bool_t isInEneTwrCut  = (eTwr > eTwrMin);
      const Bool_t isInCorrTwrCut = (eCorr >= eCorrMin);
      if (isTrigger)
        continue;
      if (!isInEtaTwrCut || !isInEneTwrCut/* || !isInCorrTwrCut*/)
        continue;

      // fill histograms
      const Bool_t isRecoilTwr = (TMath::Abs(dPhi - TMath::Pi()) < fRecoilMax);
      hTwrEne[0]  -> Fill(eTwr);
      hTwrCorr[0] -> Fill(eCorr);
      hTwrPt[0]   -> Fill(pTtwr);
      hTwrEta[0]  -> Fill(hTwr);
      hTwrPhi[0]  -> Fill(fTwr);
      hTwrDeta[0] -> Fill(dEta);
      hTwrDphi[0] -> Fill(dPhi);
      if (isRecoilTwr) {
        hTwrEne[1]  -> Fill(eTwr);
        hTwrCorr[1] -> Fill(eCorr);
        hTwrPt[1]   -> Fill(pTtwr);
        hTwrEta[1]  -> Fill(hTwr);
        hTwrPhi[1]  -> Fill(fTwr);
        hTwrDeta[1] -> Fill(dEta);
        hTwrDphi[1] -> Fill(dPhi);
      }

    }  // end tower loop 2

  }  // end event loop

  cout << "    Event loop finished: " << nTrg << " triggers found." << endl;


/*
  // normalize histograms
  const Double_t eBin  = (e2 - e1) / (Double_t) nE;
  const Double_t hBin  = (h2 - h1) / (Double_t) nH;
  const Double_t fBin  = (f2 - f1) / (Double_t) nF;
  const Double_t eNorm = eBin * nTrg;
  const Double_t hNorm = hBin * nTrg;
  const Double_t fNorm = fBin * nTrg;
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    hTrkEne[iDist]  -> Scale(1. / eNorm);
    hTrkPt[iDist]   -> Scale(1. / eNorm);
    hTrkEta[iDist]  -> Scale(1. / hNorm);
    hTrkPhi[iDist]  -> Scale(1. / fNorm);
    hTrkDeta[iDist] -> Scale(1. / hNorm);
    hTrkDphi[iDist] -> Scale(1. / fNorm);
    hTwrEne[iDist]  -> Scale(1. / eNorm);
    hTwrCorr[iDist] -> Scale(1. / eNorm);
    hTwrPt[iDist]   -> Scale(1. / eNorm);
    hTwrEta[iDist]  -> Scale(1. / hNorm);
    hTwrPhi[iDist]  -> Scale(1. / fNorm);
    hTwrDeta[iDist] -> Scale(1. / hNorm);
    hTwrDphi[iDist] -> Scale(1. / fNorm);
  }
  cout << "    Distributions scaled." << endl;
*/


  // make directories
  TDirectory *dEps;
  TDirectory *dHots;
  TDirectory *dTrks[nDist];
  TDirectory *dTwrs[nDist];
  dEps  = (TDirectory*) fOutput -> mkdir("EoverP");
  dHots = (TDirectory*) fOutput -> mkdir("HotTowers");
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    dTrks[iDist] = (TDirectory*) fOutput -> mkdir(sTrkDirs[iDist].Data());
    dTwrs[iDist] = (TDirectory*) fOutput -> mkdir(sTwrDirs[iDist].Data());
  }
  cout << "    Directories made." << endl;


  fOutput     -> cd();
  hNumTrg     -> Write();
  hTrgEt      -> Write();
  hTrgEne     -> Write();
  hTrgEta     -> Write();
  dEps        -> cd();
  hEpTwrEne   -> Write();
  hEpTwrEneL5 -> Write();
  for (Long64_t iPtBin = 0; iPtBin < nPtBins; iPtBin++) {
    hEpTrkPtE[iPtBin]   -> Write();
    hEpTrkPtH[iPtBin]   -> Write();
    hEpElec[iPtBin]     -> Write();
    hEpHadr[iPtBin]     -> Write();
    hEpTrkPtEL5[iPtBin] -> Write();
    hEpTrkPtHL5[iPtBin] -> Write();
    hEpElecL5[iPtBin]   -> Write();
    hEpHadrL5[iPtBin]   -> Write();
  }
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    dTrks[iDist]    -> cd();
    hTrkEne[iDist]  -> Write();
    hTrkPt[iDist]   -> Write();
    hTrkEta[iDist]  -> Write();
    hTrkPhi[iDist]  -> Write();
    hTrkDeta[iDist] -> Write();
    hTrkDphi[iDist] -> Write();
    dTwrs[iDist]    -> cd();
    hTwrEne[iDist]  -> Write();
    hTwrCorr[iDist] -> Write();
    hTwrPt[iDist]   -> Write();
    hTwrEta[iDist]  -> Write();
    hTwrPhi[iDist]  -> Write();
    hTwrDeta[iDist] -> Write();
    hTwrDphi[iDist] -> Write();
  }
  dHots   -> cd();
  for (Long64_t iEbin = 0; iEbin < nEbins; iEbin++) {
    hHotTwrID[iEbin]      -> Write();
    hHotTwrEne[iEbin]     -> Write();
    hHotTwrEta[iEbin]     -> Write();
    hHotTwrIDvsEne[iEbin] -> Write();
    hHotTwrIDvsEta[iEbin] -> Write();
  }
  fOutput     -> Close();
  fInput      -> cd();
  fInput      -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
