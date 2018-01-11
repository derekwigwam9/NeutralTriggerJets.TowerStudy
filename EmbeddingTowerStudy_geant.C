// 'EmbeddingTowerStudy_geant.C'
// Derek Anderson
// 09.12.2017
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
static const UInt_t nEffPar(5);
static const UInt_t nTrkMax(1000);
static const UInt_t nTwrMax(10000);
static const UInt_t nMatMax(100);
// io parameters
static const TString sTree("GfmtoDst_gnt");
static const TString sInDefault("/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/embedding/pt35_-1.ePcalc2.root");
static const TString sOutDefault("pp200r12pt35.test.root");
static const TString sTrkDirs[nDist + 1] = {"AllTracks", "RecoilTracks", "Track2D"};
static const TString sTwrDirs[nDist + 1] = {"AllTowers", "RecoilTowers", "Tower2D"};



void EmbeddingTowerStudy_geant(const Bool_t isInBatchMode=false, const TString sInput=sInDefault, const TString sOutput=sOutDefault) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning tower script..." << endl;

  // constants
  const Double_t fRecoilMax = TMath::PiOver4();
  const Double_t mPion      = 0.140;
  const Double_t dFmax      = 0.025;
  const Double_t dHmax      = 0.025;
  const Double_t zVtxMax    = 55.;
  const Double_t rVtxMax    = 2.;
  const Double_t hTrgMax    = 0.9;
  const Double_t eTtrgMin   = 9.0;
  const Double_t eTtrgMax   = 20.;
  const Double_t qTrg       = 0.;
  const Double_t hTrkMax    = 1.0;
  const Double_t pTtrkMin   = 0.2;
  const Int_t    gIDtrg     = 7;
  const Int_t    gIDpip     = 8;
  const Int_t    gIDpim     = 9;
  const Int_t    gIDkap     = 11;
  const Int_t    gIDkam     = 12;
  const Int_t    gIDpop     = 14;
  const Int_t    gIDpom     = 15;
  const Int_t    gIDelp     = 2;
  const Int_t    gIDelm     = 2;


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


  fOutput -> cd();
  // trigger histograms
  TH1D *hNumTrg;
  TH1D *hTrgEt;
  TH1D *hTrgEne;
  TH1D *hTrgEta;
  TH1D *hEffPt[nEffPar];
  // charged track (triggered) histograms
  TH1D *hTrkEne[nDist];
  TH1D *hTrkPt[nDist];
  TH1D *hTrkPtH1[nDist];
  TH1D *hTrkPtH05[nDist];
  TH1D *hTrkEta[nDist];
  TH1D *hTrkPhi[nDist];
  TH1D *hTrkChrg[nDist];
  TH1D *hTrkDeta[nDist];
  TH1D *hTrkDphi[nDist];
  TH2D *hTrkPtVsEta;
  TH2D *hTrkPtVsPhi;
  TH2D *hTrkPtVsDeta;
  TH2D *hTrkPtVsDphi;
  // neutral track (triggered) histograms
  TH1D *hTwrEne[nDist];
  TH1D *hTwrPt[nDist];
  TH1D *hTwrEta[nDist];
  TH1D *hTwrPhi[nDist];
  TH1D *hTwrChrg[nDist];
  TH1D *hTwrDeta[nDist];
  TH1D *hTwrDphi[nDist];
  TH2D *hTwrPtVsEta;
  TH2D *hTwrPtVsPhi;
  TH2D *hTwrPtVsDeta;
  TH2D *hTwrPtVsDphi;

  const UInt_t   nN  = 200;
  const UInt_t   nE  = 1000;
  const UInt_t   nH  = 40;
  const UInt_t   nDh = 80;
  const UInt_t   nF  = 720;
  const UInt_t   nDf = 720;
  const UInt_t   nQ = 6;
  const Double_t n1  = 0.;
  const Double_t n2  = 200.;
  const Double_t e1  = 0.;
  const Double_t e2  = 100.;
  const Double_t h1  = -2.;
  const Double_t h2  = 2.;
  const Double_t dH1 = -4.;
  const Double_t dH2 = 4.;
  const Double_t f1  = -1. * TMath::TwoPi();
  const Double_t f2  = TMath::TwoPi();
  const Double_t dF1 = -1. * TMath::TwoPi();
  const Double_t dF2 = TMath::TwoPi();
  const Double_t q1 = -3.;
  const Double_t q2 = 3.;
  // trigger histograms
  hNumTrg      = new TH1D("hNumTrg", "No. of triggers", nN, n1, n2);
  hTrgEt       = new TH1D("hTrgEt", "Trigger E_{T}", nE, e1, e2);
  hTrgEne      = new TH1D("hTrgEne", "Trigger energy", nE, e1, e2);
  hTrgEta      = new TH1D("hTrgEta", "Trigger #eta", nH, h1, h2);
  hTrgPhi      = new TH1D("hTrgPhi", "Trigger #varphi", nF, f1, f2);
  hTrgChrg     = new TH1D("hTrgChrg", "Trigger charge", nQ, q1, q2);
  hEffPt[0]    = new TH1D("hPionPt", "#pi^{#pm} p_{T}", nE, e1, e2);
  hEffPt[1]    = new TH1D("hKaonPt", "K^{#pm} p_{T}", nE, e1, e2);
  hEffPt[2]    = new TH1D("hProtonPt", "p^{#pm} p_{T}", nE, e1, e2);
  hEffPt[3]    = new TH1D("hElectronPt", "e^{#pm} p_{T}", nE, e1, e2);
  hEffPt[4]    = new TH1D("hTotalPt", "All particles p_{T}", nE, e1, e2);
  // triggered track histograms
  hTrkEne[0]   = new TH1D("hTrkEneAll", "Charged Track energy, all", nE, e1, e2);
  hTrkEne[1]   = new TH1D("hTrkEneRec", "Charged Track energy, recoil", nE, e1, e2);
  hTrkPt[0]    = new TH1D("hTrkPtAll", "Charged Track p_{T}, all", nE, e1, e2);
  hTrkPt[1]    = new TH1D("hTrkPtRec", "Charged Track p_{T}, recoil", nE, e1, e2);
  hTrkPtH1[0]  = new TH1D("hTrkPtH1all", "Charged track p_{T}, all w/ |#eta| < 1", nE, e1, e2);
  hTrkPtH1[1]  = new TH1D("hTrkPtH1rec", "Charged track p_{T}, recoil w/ |#eta| < 1", nE, e1, e2);
  hTrkPtH05[0] = new TH1D("hTrkPtH05all", "Charged track p_{T}, all w/ |#eta| < 0.5", nE, e1, e2);
  hTrkPtH05[1] = new TH1D("hTrkPtH05rec", "Charged track p_{T}, recoil w/ |#eta| < 0.5", nE, e1, e2);
  hTrkEta[0]   = new TH1D("hTrkEtaAll", "Charged Track #eta, all", nH, h1, h2);
  hTrkEta[1]   = new TH1D("hTrkEtaRec", "Charged Track #eta, recoil", nH, h1, h2);
  hTrkPhi[0]   = new TH1D("hTrkPhiAll", "Charged Track #varphi, all", nF, f1, f2);
  hTrkPhi[1]   = new TH1D("hTrkPhiRec", "Charged Track #varphi, recoil", nF, f1, f2);
  hTrkChrg[0]  = new TH1D("hTrkChrgAll", "Charged track charge, all", nQ, q1, q2);
  hTrkChrg[1]  = new TH1D("hTrkChrgRec", "Charged track charge, recoil", nQ, q1, q2);
  hTrkDeta[0]  = new TH1D("hTrkDetaAll", "Charged Track #Delta#eta, all", nDh, dH1, dH2);
  hTrkDeta[1]  = new TH1D("hTrkDetaRec", "Charged Track #Delta#eta, recoil", nDh, dH1, dH2);
  hTrkDphi[0]  = new TH1D("hTrkDphiAll", "Charged Track #Delta#varphi, all", nDf, dF1, dF2);
  hTrkDphi[1]  = new TH1D("hTrkDphiRec", "Charged Track #Delta#varphi, recoil", nDf, dF1, dF2);
  hTrkPtVsEta  = new TH2D("hTrkPtVsEta", "Track p_{T} vs. #eta", nH, h1, h2, nE, e1, e2);
  hTrkPtVsPhi  = new TH2D("hTrkPtVsPhi", "Track p_{T} vs. #varphi", nF, f1, f2, nE, e1, e2);
  hTrkPtVsDeta = new TH2D("hTrkPtVsDeta", "Track p_{T} vs. #Delta#eta", nDh, dH1, dH2, nE, e1, e2);
  hTrkPtVsDphi = new TH2D("hTrkPtVsDphi", "Track p_{T} vs. #Delta#varphi", nDf, dF1, dF2, nE, e1, e2);
  // triggered tower histograms
  hTwrEne[0]   = new TH1D("hTwrEneAll", "Neutral Track energy, all", nE, e1, e2);
  hTwrEne[1]   = new TH1D("hTwrEneRec", "Neutral Track energy, recoil", nE, e1, e2);
  hTwrPt[0]    = new TH1D("hTwrPtAll", "Neutral Track p_{T}, all", nE, e1, e2);
  hTwrPt[1]    = new TH1D("hTwrPtRec", "Neutral Track p_{T}, recoil", nE, e1, e2);
  hTwrEta[0]   = new TH1D("hTwrEtaAll", "Neutral Track #eta, all", nH, h1, h2);
  hTwrEta[1]   = new TH1D("hTwrEtaRec", "Neutral Track #eta, recoil", nH, h1, h2);
  hTwrPhi[0]   = new TH1D("hTwrPhiAll", "Neutral Track #varphi, all", nF, f1, f2);
  hTwrPhi[1]   = new TH1D("hTwrPhiRec", "Neutral Track #varphi, recoil", nF, f1, f2);
  hTwrChrg[0]  = new TH1D("hTwrChrgAll", "Neutral track charge, all", nQ, q1, q2);
  hTwrChrg[1]  = new TH1D("hTwrChrgRec", "Neutral track charge, recoil", nQ, q1, q2);
  hTwrDeta[0]  = new TH1D("hTwrDetaAll", "Neutral Track #Delta#eta, all", nDh, dH1, dH2);
  hTwrDeta[1]  = new TH1D("hTwrDetaRec", "Neutral Track #Delta#eta, recoil", nDh, dH1, dH2);
  hTwrDphi[0]  = new TH1D("hTwrDphiAll", "Neutral Track #Delta#varphi, all", nDf, dF1, dF2);
  hTwrDphi[1]  = new TH1D("hTwrDphiRec", "Neutral Track #Delta#varphi, recoil", nDf, dF1, dF2);
  hTwrPtVsEta  = new TH2D("hTwrPtVsEta", "Tower p_{T} vs. #eta", nH, h1, h2, nE, e1, e2);
  hTwrPtVsPhi  = new TH2D("hTwrPtVsPhi", "Tower p_{T} vs. #varphi", nF, f1, f2, nE, e1, e2);
  hTwrPtVsDeta = new TH2D("hTwrPtVsDeta", "Tower p_{T} vs. #Delta#eta", nDh, dH1, dH2, nE, e1, e2);
  hTwrPtVsDphi = new TH2D("hTwrPtVsDphi", "Tower p_{T} vs. #Delta#varphi", nDf, dF1, dF2, nE, e1, e2);
  // errors
  hTrgEt  -> Sumw2();
  hTrgEne -> Sumw2();
  hTrgEta -> Sumw2();
  hTrgPhi -> Sumw2();
  for (Long64_t iEffPar = 0; iEffPar < nEffPar; iEffPar++) {
    hEffPt[iEffPar] -> Sumw2();
  }
  for (Long64_t iDist = 0; iDist < nDist; iDist++) {
    hTrkEne[iDist]   -> Sumw2();
    hTrkPt[iDist]    -> Sumw2();
    hTrkPtH1[iDist]  -> Sumw2();
    hTrkPtH05[iDist] -> Sumw2();
    hTrkEta[iDist]   -> Sumw2();
    hTrkPhi[iDist]   -> Sumw2();
    hTrkDeta[iDist]  -> Sumw2();
    hTrkDphi[iDist]  -> Sumw2();
    hTwrEne[iDist]   -> Sumw2();
    hTwrPt[iDist]    -> Sumw2();
    hTwrEta[iDist]   -> Sumw2();
    hTwrPhi[iDist]   -> Sumw2();
    hTwrDeta[iDist]  -> Sumw2();
    hTwrDphi[iDist]  -> Sumw2();
  }
  hTrkPtVsEta  -> Sumw2();
  hTrkPtVsPhi  -> Sumw2();
  hTrkPtVsDeta -> Sumw2();
  hTrkPtVsDphi -> Sumw2();
  hTwrPtVsEta  -> Sumw2();
  hTwrPtVsPhi  -> Sumw2();
  hTwrPtVsDeta -> Sumw2();
  hTwrPtVsDphi -> Sumw2();


  Long64_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " evts. to process." << endl;

  // event loop
  Int_t    iTrig    = 0;
  Bool_t   foundTrg = false;
  Long64_t nTrig    = 0;
  Long64_t nBytes   = 0;
  Double_t fTrig    = 0;
  Double_t hTrig    = 0;
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
    iTrig    = 0;
    fTrig    = 0;
    hTrig    = 0;
    foundTrg = false;

    // event info
    const UInt_t   nTrks = (UInt_t) nPrimaryTracks;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isInZCut = (abs(zVtx) < zVtxMax);
    const Bool_t isInRCut = (abs(rVtx) < rVtxMax);
    //if (!isInRCut || !isInZCut) continue;


    // track loop 1
    for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    gIDtrk = PrimaryTrackArray_geantId[iTrk];
      const Double_t hTrk   = PrimaryTrackArray_eta[iTrk];
      const Double_t qTrk   = PrimaryTrackArray_charge[iTrk];
      const Double_t pTtrk  = PrimaryTrackArray_pT[iTrk];
      const Double_t pXtrk  = PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk  = PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk  = PrimaryTrackArray_pZ[iTrk];
      const Double_t eTrk   = TMath::Sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (mPion * mPion));
      const Double_t eTtrk  = eTrk * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTrk)));

      // calculate phi
      Double_t fTrk = TMath::ATan(pYtrk / pXtrk);
      if ((pXtrk < 0.) && (pYtrk > 0.)) fTrk += TMath::Pi();
      if ((pXtrk < 0.) && (pYtrk < 0.)) fTrk += TMath::Pi();


      // trigger cuts
      const Bool_t isInEtaTrgCut  = (TMath::Abs(hTrk) < hTrgMax);
      const Bool_t isInPtTrgCut   = ((eTtrk > eTtrgMin) && (eTtrk < eTtrgMax));
      const Bool_t isInChrgTrgCut = (qTrk == qTrg);
      const Bool_t isInPidTrgCut  = (gIDtrk == gIDtrg);
      if (isInEtaTrgCut && isInPtTrgCut && isInChrgTrgCut && isInPidTrgCut) {
        hTrgEt   -> Fill(pTtrk);
        hTrgEne  -> Fill(eTrk);
        hTrgEta  -> Fill(hTrk);
        hTrgPhi  -> Fill(fTrk);
        hTrgChrg -> Fill(qTrk);
        iTrig    = iTrk;
        fTrig    = fTrk;
        hTrig    = hTrk;
        foundTrg = true;
        break;
      }
      else {
        continue;
      }

    }  // end track loop 1
    if (foundTrg) {
      nTrig++;
    }
    else {
      //continue;
    }


    // track loop 2 
    for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    gIDtrk = PrimaryTrackArray_geantId[iTrk];
      const Double_t hTrk   = PrimaryTrackArray_eta[iTrk];
      const Double_t qTrk   = PrimaryTrackArray_charge[iTrk];
      const Double_t pXtrk  = PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk  = PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk  = PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk  = PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk   = TMath::Sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (mPion * mPion));

      // calculate phi
      Double_t fTrk = TMath::ATan(pYtrk / pXtrk);
      if ((pXtrk < 0.) && (pYtrk > 0.)) fTrk += TMath::Pi();
      if ((pXtrk < 0.) && (pYtrk < 0.)) fTrk += TMath::Pi();

      // calculate delta-eta,phi
      Double_t dHtrk = hTrk - hTrig;
      Double_t dFtrk = fTrk - fTrig;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();


      // track cuts
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtTrkCut  = (pTtrk > pTtrkMin);
      const Bool_t isTrigger     = (iTrk == iTrig);
      if (!isInEtaTrkCut || !isInPtTrkCut || isTrigger) continue;

      // fill histograms
      const Bool_t isPion      = ((gIDtrk == gIDpip) || (gIDtrk == gIDpim));
      const Bool_t isKaon      = ((gIDtrk == gIDkap) || (gIDtrk == gIDkam));
      const Bool_t isProton    = ((gIDtrk == gIDpop) || (gIDtrk == gIDpom));
      const Bool_t isElectron  = ((gIDtrk == gIDelp) || (gIDtrk == gIDelm));
      const Bool_t isCharged   = (qTrk != 0.);
      const Bool_t isRecoilTrk = (TMath::Abs(dFtrk - TMath::Pi()) < fRecoilMax);
       // charged tracks
      if (isCharged) {
        hTrkEne[0]   -> Fill(eTrk);
        hTrkPt[0]    -> Fill(pTtrk);
        hTrkEta[0]   -> Fill(hTrk);
        hTrkPhi[0]   -> Fill(fTrk);
        hTrkChrg[0]  -> Fill(qTrk);
        hTrkDeta[0]  -> Fill(dHtrk);
        hTrkDphi[0]  -> Fill(dFtrk);
        hTrkPtVsEta  -> Fill(hTrk, pTtrk);
        hTrkPtVsPhi  -> Fill(fTrk, pTtrk);
        hTrkPtVsDeta -> Fill(dHtrk, pTtrk);
        hTrkPtVsDphi -> Fill(dFtrk, pTtrk);
        if (TMath::Abs(hTrk) < 1.)  hTrkPtH1[0]  -> Fill(pTtrk);
        if (TMath::Abs(hTrk) < 0.5) hTrkPtH05[0] -> Fill(pTtrk);

        // recoil tracks
        if (isRecoilTrk) {
          hTrkEne[1]  -> Fill(eTrk);
          hTrkPt[1]   -> Fill(pTtrk);
          hTrkEta[1]  -> Fill(hTrk);
          hTrkPhi[1]  -> Fill(fTrk);
          hTrkChrg[1] -> Fill(qTrk);
          hTrkDeta[1] -> Fill(dHtrk);
          hTrkDphi[1] -> Fill(dFtrk);
          if (TMath::Abs(hTrk) < 1.)  hTrkPtH1[1]  -> Fill(pTtrk);
          if (TMath::Abs(hTrk) < 0.5) hTrkPtH05[1] -> Fill(pTtrk);
        }

        // individual species
        if (isPion)     hEffPt[0] -> Fill(pTtrk);
        if (isKaon)     hEffPt[1] -> Fill(pTtrk);
        if (isProton)   hEffPt[2] -> Fill(pTtrk);
        if (isElectron) hEffPt[3] -> Fill(pTtrk);
        hEffPt[4] -> Fill(pTtrk);
      }
      // neutral tracks
      else {
        hTwrEne[0]   -> Fill(eTrk);
        hTwrPt[0]    -> Fill(pTtrk);
        hTwrEta[0]   -> Fill(hTrk);
        hTwrPhi[0]   -> Fill(fTrk);
        hTwrChrg[0]  -> Fill(qTrk);
        hTwrDeta[0]  -> Fill(dHtrk);
        hTwrDphi[0]  -> Fill(dFtrk);
        hTwrPtVsEta  -> Fill(hTrk, pTtrk);
        hTwrPtVsPhi  -> Fill(fTrk, pTtrk);
        hTwrPtVsDeta -> Fill(dHtrk, pTtrk);
        hTwrPtVsDphi -> Fill(dFtrk, pTtrk);
        if (isRecoilTrk) {
          hTwrEne[1]  -> Fill(eTrk);
          hTwrPt[1]   -> Fill(pTtrk);
          hTwrEta[1]  -> Fill(hTrk);
          hTwrPhi[1]  -> Fill(fTrk);
          hTwrChrg[1] -> Fill(qTrk);
          hTwrDeta[1] -> Fill(dHtrk);
          hTwrDphi[1] -> Fill(dFtrk);
        }
      }

    }  // end track loop 2

    hNumTrg -> Fill(nTrig);

  }  // end event loop

  cout << "    Event loop finished: " << nTrig << " triggers found." << endl;


  // make directories
  TDirectory *dTrks[nDist + 1];
  TDirectory *dTwrs[nDist + 1];
  for (Long64_t iDist = 0; iDist < nDist + 1; iDist++) {
    dTrks[iDist] = (TDirectory*) fOutput -> mkdir(sTrkDirs[iDist].Data());
    dTwrs[iDist] = (TDirectory*) fOutput -> mkdir(sTwrDirs[iDist].Data());
  }
  cout << "    Directories made." << endl;


  fOutput  -> cd();
  hNumTrg  -> Write();
  hTrgEt   -> Write();
  hTrgEne  -> Write();
  hTrgEta  -> Write();
  hTrgPhi  -> Write();
  hTrgChrg -> Write();
  for (Long64_t iEffPar = 0; iEffPar < nEffPar; iEffPar++) {
    hEffPt[iEffPar] -> Write();
  }
  for (Long64_t iDist = 0; iDist < nDist + 1; iDist++) {
    dTrks[iDist] -> cd();
    if (iDist == nDist) {
      hTrkPtVsEta  -> Write();
      hTrkPtVsPhi  -> Write();
      hTrkPtVsDeta -> Write();
      hTrkPtVsDphi -> Write();
    }
    else {
      hTrkEne[iDist]   -> Write();
      hTrkPt[iDist]    -> Write();
      hTrkPtH1[iDist]  -> Write();
      hTrkPtH05[iDist] -> Write();
      hTrkEta[iDist]   -> Write();
      hTrkPhi[iDist]   -> Write();
      hTrkChrg[iDist]  -> Write();
      hTrkDeta[iDist]  -> Write();
      hTrkDphi[iDist]  -> Write();
    }
    dTwrs[iDist] -> cd();
    if (iDist == nDist) {
      hTwrPtVsEta  -> Write();
      hTwrPtVsPhi  -> Write();
      hTwrPtVsDeta -> Write();
      hTwrPtVsDphi -> Write();
    }
    else {
      hTwrEne[iDist]  -> Write();
      hTwrPt[iDist]   -> Write();
      hTwrEta[iDist]  -> Write();
      hTwrPhi[iDist]  -> Write();
      hTwrChrg[iDist] -> Write();
      hTwrDeta[iDist] -> Write();
      hTwrDphi[iDist] -> Write();
    }
  }
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
