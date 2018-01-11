// 'EmbeddingTowerStudy_hot.C'
// Derek Anderson
// 09.29.2017
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
static const UInt_t nTrkHist(2);
static const UInt_t nEneBins(6);
static const UInt_t nTrkMax(1000);
static const UInt_t nTwrMax(10000);
static const UInt_t nMatMax(100);
// io parameters
static const TString sTree("GfmtoDst_mu");
static const TString sInDefault("/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/embedding/pt5_-1.ePcalc2.root");
static const TString sOutDefault("pp200r12pt5g.trgVsClustEta.root");



void EmbeddingTowerStudy_hot(const Bool_t isInBatchMode=false, const TString sInput=sInDefault, const TString sOutput=sOutDefault) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning tower script..." << endl;

  // constants
  const Float_t  eneCut[nEneBins]  = {0., 2., 5., 10., 20., 35.};
  const TString  sTrkEnd[nTrkHist] = {"_beforeQA", "_afterQA"};
  const TString  sEneEnd[nEneBins] = {"_e0", "_e2", "_e5", "_e10", "_e20", "_e35"};
  const Double_t dEdXunits         = 10e5;

  // cuts
  const Int_t    adcMax   = 6004;
  const UInt_t   nFitMin  = 15;
  const Float_t  rFitMin  = 0.52;
  const Float_t  dcaMax   = 1.;
  const Float_t  etaMax   = 1.;
  const Float_t  pTtrkMin = 0.2;
  const Float_t  pMatMin  = 2.;
  const Float_t  elecMin  = 3.4;
  const Float_t  elecMax  = 5.;
  const Double_t pProjMax = 3.;
  const Double_t hTrgMax  = 0.9;
  const Double_t eTwrMin  = 0.2;
  const Double_t eTtrgMin = 9.;
  const Double_t eTtrgMax = 20.;


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
  // event histograms
  TH1D *hVtxZ;
  TH1D *hVtxR;
  TH1D *hTrgEta;
  TH1D *hClustEta;
  TH1D *hTrgPhi;
  TH1D *hTrgEt;
  TH2D *hTrgPhiVsEta;
  TH2D *hClustVsTrgEta;
  // track histograms
  TH1D *hTrkP[nTrkHist];
  TH1D *hTrkPtH1[nTrkHist];
  TH1D *hTrkPtH05[nTrkHist];
  TH1D *hTrkEta[nTrkHist];
  TH1D *hTrkPhi[nTrkHist];
  TH1D *hTrkDeDx[nTrkHist];
  TH2D *hTrkPhiVsEta[nTrkHist];
  TH2D *hTrkDeDxVsP[nTrkHist];
  // tower histograms
  TH1D *hTwrId[nEneBins];
  TH1D *hTwrEne[nEneBins];
  TH1D *hTwrEta[nEneBins];
  TH1D *hTwrPhi[nEneBins];
  TH1D *hTwrPe[nEneBins];
  TH1D *hTwrEp[nEneBins][2];
  TH2D *hTwrEneVsId[nEneBins];
  TH2D *hTwrEtaVsId[nEneBins];
  TH2D *hTwrPhiVsId[nEneBins];
  TH2D *hTwrPhiVsEta[nEneBins];

  const UInt_t   nId = 4802;
  const UInt_t   nDx = 1000;
  const UInt_t   nV  = 2000;
  const UInt_t   nE  = 2000;
  const UInt_t   nP  = 1000;
  const UInt_t   nPx = 5000;
  const UInt_t   nH  = 80;
  const UInt_t   nF  = 120;
  const UInt_t   nPe = 500;
  const UInt_t   nEp = 100;
  const Double_t id[2] = {-1., 4801.};
  const Double_t dX[2] = {0., 0.0001};
  const Double_t v[2]  = {-100., 100.};
  const Double_t e[2]  = {0., 200.};
  const Double_t p[2]  = {0., 100.};
  const Double_t pX[2] = {0., 5.};
  const Double_t h[2]  = {-2., 2.};
  const Double_t f[2]  = {-1. * TMath::Pi(), TMath::Pi()};
  const Double_t pE[2] = {0., 50.};
  const Double_t eP[2] = {0., 10.};
  hVtxZ          = new TH1D("hVtxZ", "", nV, v[0], v[1]);
  hVtxR          = new TH1D("hVtxR", "", nV, v[0], v[1]);
  hTrgEta        = new TH1D("hTrgEta", "", nH, h[0], h[1]);
  hClustEta      = new TH1D("hClustEta", "", nH, h[0], h[1]);
  hTrgPhi        = new TH1D("hTrgPhi", "", nF, f[0], f[1]);
  hTrgEt         = new TH1D("hTrgEt", "", nP, p[0], p[1]);
  hTrgPhiVsEta   = new TH2D("hTrgPhiVsEta", "", nH, h[0], h[1], nF, f[0], f[1]);
  hClustVsTrgEta = new TH2D("hClustVsTrgEta", "Trigger tower vs. cluster #eta;#eta^{twr};#eta^{clust}", nH, h[0], h[1], nH, h[0], h[1]);
  hVtxZ          -> Sumw2();
  hVtxR          -> Sumw2();
  hTrgEta        -> Sumw2();
  hClustEta      -> Sumw2();
  hTrgPhi        -> Sumw2();
  hTrgEt         -> Sumw2();
  hTrgPhiVsEta   -> Sumw2();
  hClustVsTrgEta -> Sumw2();
  for (Long64_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {
    TString sIdName("hTwrId");
    TString sEneName("hTwrEne");
    TString sEtaName("hTwrEta");
    TString sPhiName("hTwrPhi");
    TString sPeName("hTwrPe");
    TString sEpEname("hTwrEpElec");
    TString sEpHname("hTwrEpHadr");
    TString sEneVsId("hTwrEneVsId");
    TString sEtaVsId("hTwrEtaVsId");
    TString sPhiVsId("hTwrPhiVsId");
    TString sPhiVsEta("hTwrPhiVsEta");
    sIdName.Append(sEneEnd[iEneBin].Data());
    sEneName.Append(sEneEnd[iEneBin].Data());
    sEtaName.Append(sEneEnd[iEneBin].Data());
    sPhiName.Append(sEneEnd[iEneBin].Data());
    sPeName.Append(sEneEnd[iEneBin].Data());
    sEpEname.Append(sEneEnd[iEneBin].Data());
    sEpHname.Append(sEneEnd[iEneBin].Data());
    sEneVsId.Append(sEneEnd[iEneBin].Data());
    sEtaVsId.Append(sEneEnd[iEneBin].Data());
    sPhiVsId.Append(sEneEnd[iEneBin].Data());
    sPhiVsEta.Append(sEneEnd[iEneBin].Data());
    hTwrId[iEneBin]       = new TH1D(sIdName.Data(), "", nId, id[0], id[1]);
    hTwrEne[iEneBin]      = new TH1D(sEneName.Data(), "", nE, e[0], e[1]);
    hTwrEta[iEneBin]      = new TH1D(sEtaName.Data(), "", nH, h[0], h[1]);
    hTwrPhi[iEneBin]      = new TH1D(sPhiName.Data(), "", nF, f[0], f[1]);
    hTwrPe[iEneBin]       = new TH1D(sPeName.Data(), "", nPe, pE[0], pE[1]);
    hTwrEp[iEneBin][0]    = new TH1D(sEpEname.Data(), "", nEp, eP[0], eP[1]);
    hTwrEp[iEneBin][1]    = new TH1D(sEpHname.Data(), "", nEp, eP[0], eP[1]);
    hTwrEneVsId[iEneBin]  = new TH2D(sEneVsId.Data(), "", nId, id[0], id[1], nE, e[0], e[1]);
    hTwrEtaVsId[iEneBin]  = new TH2D(sEtaVsId.Data(), "", nId, id[0], id[1], nH, h[0], h[1]);
    hTwrPhiVsId[iEneBin]  = new TH2D(sPhiVsId.Data(), "", nId, id[0], id[1], nF, f[0], f[1]);
    hTwrPhiVsEta[iEneBin] = new TH2D(sPhiVsEta.Data(), "", nH, h[0], h[1], nF, f[0], f[1]);
    hTwrId[iEneBin]       -> Sumw2();
    hTwrEne[iEneBin]      -> Sumw2();
    hTwrEta[iEneBin]      -> Sumw2();
    hTwrPhi[iEneBin]      -> Sumw2();
    hTwrPe[iEneBin]       -> Sumw2();
    hTwrEp[iEneBin][0]    -> Sumw2();
    hTwrEp[iEneBin][1]    -> Sumw2();
    hTwrEneVsId[iEneBin]  -> Sumw2();
    hTwrEtaVsId[iEneBin]  -> Sumw2();
    hTwrPhiVsId[iEneBin]  -> Sumw2();
    hTwrPhiVsEta[iEneBin] -> Sumw2();
  }
  for (Long64_t iTrkHist = 0; iTrkHist < nTrkHist; iTrkHist++) {
    TString sPname("hTrkP");
    TString sPtH1name("hTrkPtH1");
    TString sPtH05name("hTrkPtH05");
    TString sEtaName("hTrkEta");
    TString sPhiName("hTrkPhi");
    TString sDeDxName("hTrkDeDx");
    TString sPhiVsEta("hTrkPhiVsEta");
    TString sDeDxVsP("hTrkDeDxVsP");
    sPname.Append(sTrkEnd[iTrkHist].Data());
    sPtH1name.Append(sTrkEnd[iTrkHist].Data());
    sPtH05name.Append(sTrkEnd[iTrkHist].Data());
    sEtaName.Append(sTrkEnd[iTrkHist].Data());
    sPhiName.Append(sTrkEnd[iTrkHist].Data());
    sDeDxName.Append(sTrkEnd[iTrkHist].Data());
    sPhiVsEta.Append(sTrkEnd[iTrkHist].Data());
    sDeDxVsP.Append(sTrkEnd[iTrkHist].Data());
    hTrkP[iTrkHist]        = new TH1D(sPname.Data(), "", nP, p[0], p[1]);
    hTrkPtH1[iTrkHist]     = new TH1D(sPtH1name.Data(), "", nP, p[0], p[1]);
    hTrkPtH05[iTrkHist]    = new TH1D(sPtH05name.Data(), "", nP, p[0], p[1]);
    hTrkEta[iTrkHist]      = new TH1D(sEtaName.Data(), "", nH, h[0], h[1]);
    hTrkPhi[iTrkHist]      = new TH1D(sPhiName.Data(), "", nF, f[0], f[1]);
    hTrkDeDx[iTrkHist]     = new TH1D(sDeDxName.Data(), "", nDx, dX[0], dX[1]);
    hTrkPhiVsEta[iTrkHist] = new TH2D(sPhiVsEta.Data(), "", nH, h[0], h[1], nF, f[0], f[1]);
    hTrkDeDxVsP[iTrkHist]  = new TH2D(sDeDxVsP.Data(), "", nPx, pX[0], pX[1], nDx, dX[0], dX[1]);
    hTrkP[iTrkHist]        -> Sumw2();
    hTrkPtH1[iTrkHist]     -> Sumw2();
    hTrkPtH05[iTrkHist]    -> Sumw2();
    hTrkEta[iTrkHist]      -> Sumw2();
    hTrkPhi[iTrkHist]      -> Sumw2();
    hTrkDeDx[iTrkHist]     -> Sumw2();
    hTrkPhiVsEta[iTrkHist] -> Sumw2();
    hTrkDeDxVsP[iTrkHist]  -> Sumw2();
  }


  Long64_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " evts. to process." << endl;

  // event loop
  Long64_t nBytes = 0;
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

    // fill histograms
    const Double_t xVtx = xVertex;
    const Double_t yVtx = yVertex;
    const Double_t zVtx = zVertex;
    const Double_t rVtx = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));
    hVtxZ -> Fill(zVtx);
    hVtxR -> Fill(rVtx);


    // track loop
    const Int_t nTrks = PrimaryTrackArray_;
    for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFit  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPos  = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFit  = (Double_t) nFit / (Double_t) nPos;
      const Double_t dca   = PrimaryTrackArray_dcag[iTrk];
      const Double_t dEdX  = PrimaryTrackArray_dEdx[iTrk];
      const Double_t hTrk  = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk  = PrimaryTrackArray_phi[iTrk];
      const Double_t pXtrk = PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk = PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk = PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk = PrimaryTrackArray_pT[iTrk];
      const Double_t pTrk  = TMath::Sqrt((pXtrk * pXtrk) + (pYtrk * pYtrk) + (pZtrk * pZtrk));

      // fill histograms
      hTrkP[0]        -> Fill(pTrk);
      hTrkEta[0]      -> Fill(hTrk);
      hTrkPhi[0]      -> Fill(fTrk);
      hTrkDeDx[0]     -> Fill(dEdX);
      hTrkPhiVsEta[0] -> Fill(hTrk, fTrk);
      hTrkDeDxVsP[0]  -> Fill(pTrk, dEdX);
      if (TMath::Abs(hTrk) < 1.)  hTrkPtH1[0]  -> Fill(pTtrk);
      if (TMath::Abs(hTrk) < 0.5) hTrkPtH05[0] -> Fill(pTtrk);

      const Bool_t isInFitCut = (nFit >= nFitMin);
      const Bool_t isInRatCut = (rFit >= rFitMin);
      const Bool_t isInDcaCut = (dca < dcaMax);
      const Bool_t isInEtaCut = (TMath::Abs(hTrk) < etaMax);
      const Bool_t isInPtCut  = (pTtrk > pTtrkMin);
      if (isInFitCut && isInRatCut && isInDcaCut && isInEtaCut) {
        hTrkP[1]        -> Fill(pTrk);
        hTrkDeDx[1]     -> Fill(dEdX);
        hTrkDeDxVsP[1]  -> Fill(pTrk, dEdX);
        if (isInPtCut) {
          hTrkEta[1]      -> Fill(hTrk);
          hTrkPhi[1]      -> Fill(fTrk);
          hTrkPhiVsEta[1] -> Fill(hTrk, fTrk);
          if (TMath::Abs(hTrk) < 1.)  hTrkPtH1[1]  -> Fill(pTtrk);
          if (TMath::Abs(hTrk) < 0.5) hTrkPtH05[1] -> Fill(pTtrk);
        }
      }

    }  // end track loop


    // tower loop
    const Int_t nTwrs = TowerArray_;
    for (Long64_t iTwr = 0; iTwr < nTwrs; iTwr++) {

      // tower info
      const Int_t    twrID  = TowerArray_TwrId[iTwr];
      const Int_t    twrADC = TowerArray_TwrADC[iTwr];
      const UInt_t   nMatch = TowerArray_NoOfmatchedTrk[iTwr];
      const UInt_t   iMatch = TowerArray_fMatchedTracksArray_[iTwr][0];
      const Double_t pProj  = TowerArray_TwrMatchP[iTwr];
      const Double_t fTwr   = TowerArray_TwrPhi[iTwr];
      const Double_t hTwr   = TowerArray_TwrEta[iTwr];
      const Double_t eTwr   = TowerArray_TwrEng[iTwr];
      const Double_t eTtwr  = eTwr * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTwr)));

      // fill trigger histograms
      const Bool_t isInAdcTrgCut  = (twrADC <= adcMax);
      const Bool_t isInProjTrgCut = (pProj < pProjMax);
      const Bool_t isInEneTrgCut  = (eTwr > eTwrMin);
      const Bool_t isInEtaTrgCut  = (TMath::Abs(hTwr) < hTrgMax);
      const Bool_t isInEtTrgCut   = ((eTtwr > eTtrgMin) && (eTtwr < eTtrgMax));
      if (isInAdcTrgCut && isInProjTrgCut && isInEneTrgCut && isInEtaTrgCut && isInEtTrgCut) {
        hTrgEta        -> Fill(hTwr);
        hClustEta      -> Fill(hTwr);
        hTrgPhi        -> Fill(fTwr);
        hTrgEt         -> Fill(eTtwr);
        hTrgPhiVsEta   -> Fill(hTwr, fTwr);
        hClustVsTrgEta -> Fill(hTwr, hTwr);
      }

      // fill tower histograms
      for (Long64_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {
        if (eTwr > eneCut[iEneBin]) {
          hTwrId[iEneBin]       -> Fill(twrID);
          hTwrEne[iEneBin]      -> Fill(eTwr);
          hTwrEta[iEneBin]      -> Fill(hTwr);
          hTwrPhi[iEneBin]      -> Fill(fTwr);
          hTwrEneVsId[iEneBin]  -> Fill(twrID, eTwr);
          hTwrEtaVsId[iEneBin]  -> Fill(twrID, hTwr);
          hTwrPhiVsId[iEneBin]  -> Fill(twrID, fTwr);
          hTwrPhiVsEta[iEneBin] -> Fill(hTwr, fTwr);
        }
      }

      // E/p and p/E calculation
      const Bool_t hasOneMatch = (nMatch == 1);
      const Bool_t indexIsGood = ((iMatch >= 0) && (iMatch < nTrks));
      if (hasOneMatch && indexIsGood) {
        // match info
        const UInt_t   nFit  = PrimaryTrackArray_nHitsFit[iMatch];
        const UInt_t   nPos  = PrimaryTrackArray_nHitsPoss[iMatch];
        const Double_t rFit  = (Double_t) nFit / (Double_t) nPos;
        const Double_t dca   = PrimaryTrackArray_dcag[iMatch];
        const Double_t dEdX  = PrimaryTrackArray_dEdx[iMatch];
        const Double_t dEdXn = dEdX * dEdXunits;
        const Double_t pXmat = PrimaryTrackArray_pX[iMatch];
        const Double_t pYmat = PrimaryTrackArray_pY[iMatch];
        const Double_t pZmat = PrimaryTrackArray_pZ[iMatch];
        const Double_t pMat  = TMath::Sqrt((pXmat * pXmat) + (pYmat * pYmat) + (pZmat * pZmat));
        const Double_t ePtwr = eTwr / pMat;
        const Double_t pEtwr = pMat / eTwr;

        // QA cuts
        const Bool_t isInFitCut  = (nFit >= nFitMin);
        const Bool_t isInRatCut  = (rFit >= rFitMin);
        const Bool_t isInDcaCut  = (dca < dcaMax);
        const Bool_t isInEtaCut  = (TMath::Abs(hTrk) < etaMax);
        const Bool_t isInPmatCut = (pMat > pMatMin);
        const Bool_t isInElecCut = ((dEdXn > elecMin) && (dEdXn < elecMax));
        if (isInFitCut && isInRatCut && isInDcaCut && isInEtaCut) {
          for (Long64_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {
            if (eTwr > eneCut[iEneBin]) {
              hTwrPe[iEneBin] -> Fill(pEtwr);
              if (isInPmatCut && isInElecCut)
                hTwrEp[iEneBin][0] -> Fill(ePtwr);
              if (isInPmatCut && !isInElecCut)
                hTwrEp[iEneBin][1] -> Fill(ePtwr);
            }
          }
        }  // end QA cut
      }  // end E/p, p/E calculation

    }  // end tower loop

  }  // end event loop

  cout << "    Event loop finished!" << endl;


  fOutput        -> cd();
  hVtxZ          -> Write();
  hVtxR          -> Write();
  hTrgEta        -> Write();
  hClustEta      -> Write();
  hTrgPhi        -> Write();
  hTrgEt         -> Write();
  hTrgPhiVsEta   -> Write();
  hClustVsTrgEta -> Write();
  for (Long64_t iTrkHist = 0; iTrkHist < nTrkHist; iTrkHist++) {
    hTrkP[iTrkHist]        -> Write();
    hTrkPtH1[iTrkHist]     -> Write();
    hTrkPtH05[iTrkHist]    -> Write();
    hTrkEta[iTrkHist]      -> Write();
    hTrkPhi[iTrkHist]      -> Write();
    hTrkDeDx[iTrkHist]     -> Write();
    hTrkPhiVsEta[iTrkHist] -> Write();
    hTrkDeDxVsP[iTrkHist]  -> Write();
  }
  for (Long64_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {
    hTwrId[iEneBin]       -> Write();
    hTwrEne[iEneBin]      -> Write();
    hTwrEta[iEneBin]      -> Write();
    hTwrPhi[iEneBin]      -> Write();
    hTwrPe[iEneBin]       -> Write();
    hTwrEp[iEneBin][0]    -> Write();
    hTwrEp[iEneBin][1]    -> Write();
    hTwrEneVsId[iEneBin]  -> Write();
    hTwrEtaVsId[iEneBin]  -> Write();
    hTwrPhiVsId[iEneBin]  -> Write();
    hTwrPhiVsEta[iEneBin] -> Write();
  }
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
