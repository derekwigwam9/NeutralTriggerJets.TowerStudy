// 'ListCheckerWhoChecksGoodAndBadValuesSperately.C'
// Derek Anderson
// 10.27.2017
//
// This takes a CALIBRATIONS list from the STAR EMC
// database, and plots the 5 ADC values for towers
// with good and bad calibration statuses seperately.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;


// global constants
static const UInt_t  nLists(2);
static const UInt_t  nLines(4800);
static const UInt_t  nAverage(100);
static const UInt_t  nAvgValues(nLines / nAverage);
static const UInt_t  nColumns(5);
static const TString sListA("lists/pp200r9ofl.calibrationsPruned.d15m12y2008.txt");
static const TString sListB("lists/pp200r12ofl.calibrationsPruned.d20m12y2011.txt");
static const TString sOutput("bemcAdcCheck.m12y2008oflXm12y2011ofl.onlyGoodStatuses.d27m10y2017.root");



void ListCheckerWhoChecksGoodAndBadValuesSeperately() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Checking lists..." << endl;


  // create output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr <<"PANIC: couldn't create output file!" << endl;
    return;
  }
  fOutput -> cd();
  cout << "    Output file created." << endl;

  // create histograms
  TH1D *hGudVal[nColumns][nLists];
  TH1D *hGudAvg[nColumns][nLists];
  TH1D *hBadVal[nColumns][nLists];
  TH1D *hBadAvg[nColumns][nLists];
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    TString sGudValA("hGoodVal");
    TString sGudAvgA("hGoodAvg");
    TString sGudValB("hGoodVal");
    TString sGudAvgB("hGoodAvg");
    TString sBadValA("hBadVal");
    TString sBadAvgA("hBadAvg");
    TString sBadValB("hBadVal");
    TString sBadAvgB("hBadAvg");
    sGudValA += iColumn;
    sGudAvgA += iColumn;
    sGudValB += iColumn;
    sGudAvgB += iColumn;
    sBadValA += iColumn;
    sBadAvgA += iColumn;
    sBadValB += iColumn;
    sBadAvgB += iColumn;
    sGudValA += "A";
    sGudAvgA += "A";
    sGudValB += "B";
    sGudAvgB += "B";
    sBadValA += "A";
    sBadAvgA += "A";
    sBadValB += "B";
    sBadAvgB += "B";
    hGudVal[iColumn][0] = new TH1D(sGudValA.Data(), "", nLines, 0., (Double_t) nLines);
    hGudVal[iColumn][1] = new TH1D(sGudValB.Data(), "", nLines, 0., (Double_t) nLines);
    hGudAvg[iColumn][0] = new TH1D(sGudAvgA.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
    hGudAvg[iColumn][1] = new TH1D(sGudAvgB.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
    hBadVal[iColumn][0]  = new TH1D(sBadValA.Data(), "", nLines, 0., (Double_t) nLines);
    hBadVal[iColumn][1]  = new TH1D(sBadValB.Data(), "", nLines, 0., (Double_t) nLines);
    hBadAvg[iColumn][0]  = new TH1D(sBadAvgA.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
    hBadAvg[iColumn][1]  = new TH1D(sBadAvgB.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
  }
  cout << "    Created histograms:\n"
       << "      " << nColumns << " quantities to compare.\n"
       << "      " << nAvgValues << " averages to calculate."
       << endl;


  // open streams
  ifstream listA(sListA.Data());
  ifstream listB(sListB.Data());
  if (!listA || !listB) {
    cerr << "PANIC: couldn't open one of the streams!" << endl;
    return;
  }
  cout << "    Opened streams, processesing..." << endl;

  // read lists, calculate averages
  UInt_t   iVal(0);
  UInt_t   iAvg(0);
  UInt_t   iLine(0);
  UInt_t   iGudAvg[nLists] = {0, 0};
  UInt_t   iBadAvg[nLists] = {0, 0};
  UInt_t   nGud[nLists]    = {0, 0};
  UInt_t   nBad[nLists]    = {0, 0};
  Double_t towerID[nLists] = {0., 0.};
  Double_t status[nLists]  = {0., 0.};
  Double_t gudVal[nLists][nLines][nColumns];
  Double_t gudAvg[nLists][nAvgValues][nColumns];
  Double_t badVal[nLists][nLines][nColumns];
  Double_t badAvg[nLists][nAvgValues][nColumns];
  while (listA && listB) {

    // read line, check status
    listA >> towerID[0];
    listB >> towerID[1];
    listA >> status[0];
    listB >> status[1];
    const Bool_t statIsGood[nLists] = {(status[0] == 1.), (status[1] == 1.)};
    for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
      gudVal[0][iLine][iColumn] = -666.;
      gudVal[1][iLine][iColumn] = -666.;
      badVal[0][iLine][iColumn] = -666.;
      badVal[1][iLine][iColumn] = -666.;
      // fill list A values
      if (statIsGood[0])  {
        listA >> gudVal[0][iLine][iColumn];
        gudAvg[0][iAvg][iColumn] += gudVal[0][iLine][iColumn];
      } else {
        listA >> badVal[0][iLine][iColumn];
        badAvg[0][iAvg][iColumn] += badVal[0][iLine][iColumn];
      }
      // fill list B values
      if (statIsGood[1])  {
        listB >> gudVal[1][iLine][iColumn];
        gudAvg[1][iAvg][iColumn] += gudVal[1][iLine][iColumn];
      } else {
        listB >> badVal[1][iLine][iColumn];
        badAvg[1][iAvg][iColumn] += badVal[1][iLine][iColumn];
      }
    }  // end column loop
    if (statIsGood[0]) {
      iGudAvg[0]++;
      nGud[0]++;
    } else {
      iBadAvg[0]++;
      nBad[0]++;
    }
    if (statIsGood[1]) {
      iGudAvg[1]++;
      nGud[1]++;
    } else {
      iBadAvg[1]++;
      nBad[1]++;
    }
    iVal++;
    iLine++;

    // check if time to exit loop
    const Bool_t loopIsFinished       = (iLine == nLines);
    const Bool_t averageIsCalculated  = (iVal == nAverage);
    const Bool_t gudAvgIsGood[nLists] = {(iGudAvg[0] > 0), (iGudAvg[1] > 0)};
    const Bool_t badAvgIsGood[nLists] = {(iBadAvg[0] > 0), (iBadAvg[1] > 0)};
    if (loopIsFinished) {
      for (UInt_t iList = 0; iList < nLists; iList++) {
        for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
          if (gudAvgIsGood[iList]) {
            gudAvg[iList][iAvg][iColumn] /= iGudAvg[iList];
          } else {
            gudAvg[iList][iAvg][iColumn] = -666.;
          }
          if (badAvgIsGood[iList]) {
            badAvg[iList][iAvg][iColumn] /= iBadAvg[iList];
          } else {
            badAvg[iList][iAvg][iColumn] = -666.;
          }
        }  // end column loop
      }  // end list loop
      break;
    }
    else if (averageIsCalculated) {
      for (UInt_t iList = 0; iList < nLists; iList++) {
        for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
          if (gudAvgIsGood[iList]) {
            gudAvg[iList][iAvg][iColumn] /= iGudAvg[iList];
          } else {
            gudAvg[iList][iAvg][iColumn] = -666.;
          }
          if (badAvgIsGood[iList]) {
            badAvg[iList][iAvg][iColumn] /= iBadAvg[iList];
          } else {
            badAvg[iList][iAvg][iColumn] = -666.;
          }
        }  // end column loop
      }  // end list loop
      iGudAvg[0] = 0;
      iGudAvg[1] = 0;
      iBadAvg[0] = 0;
      iBadAvg[1] = 0;
      iVal       = 0;
      iAvg++;
    }

  }   // end line loop

  // close streams
  listA.close();
  listB.close();
  cout << "    Done processing.\n"
       << "      List A: " << nGud[0] << " good towers, " << nBad[0] << " bad towers.\n"
       << "      List B: " << nGud[1] << " good towers, " << nBad[1] << " bad towers."
       << endl;


  // fill histograms
  const Double_t fracError = 0.2;
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    for (iLine = 0; iLine < nLines; iLine++) {
      const Double_t gudValA       = gudVal[0][iLine][iColumn];
      const Double_t gudValB       = gudVal[1][iLine][iColumn];
      const Double_t badValA       = badVal[0][iLine][iColumn];
      const Double_t badValB       = badVal[1][iLine][iColumn];
      const Double_t gudValErrA    = gudValA * fracError;
      const Double_t gudValErrB    = gudValB * fracError;
      const Double_t badValErrA    = badValA * fracError;
      const Double_t badValErrB    = badValB * fracError;
      const Bool_t   gudValAisGood = (gudValA != -666.);
      const Bool_t   gudValBisGood = (gudValB != -666.);
      const Bool_t   badValAisGood = (badValA != -666.);
      const Bool_t   badValBisGood = (badValB != -666.);
      const UInt_t   iBinVal       = iLine + 1;
      if (gudValAisGood) {
        hGudVal[iColumn][0] -> SetBinContent(iBinVal, gudValA);
        hGudVal[iColumn][0] -> SetBinError(iBinVal, gudValErrA);
      }
      if (gudValBisGood) {
        hGudVal[iColumn][1] -> SetBinContent(iBinVal, gudValB);
        hGudVal[iColumn][1] -> SetBinError(iBinVal, gudValErrB);
      }
      if (badValAisGood) {
        hBadVal[iColumn][0] -> SetBinContent(iBinVal, badValA);
        hBadVal[iColumn][0] -> SetBinError(iBinVal, badValErrA);
      }
      if (badValBisGood) {
        hBadVal[iColumn][1] -> SetBinContent(iBinVal, badValB);
        hBadVal[iColumn][1] -> SetBinError(iBinVal, badValErrB);
      }
    }  // end line loop
    for (iAvg = 0; iAvg < nAvgValues; iAvg++) {
      const Double_t gudAvgA       = gudAvg[0][iAvg][iColumn];
      const Double_t gudAvgB       = gudAvg[1][iAvg][iColumn];
      const Double_t badAvgA       = badAvg[0][iAvg][iColumn];
      const Double_t badAvgB       = badAvg[1][iAvg][iColumn];
      const Double_t gudAvgErrA    = gudAvgA * fracError;
      const Double_t gudAvgErrB    = gudAvgB * fracError;
      const Double_t badAvgErrA    = badAvgA * fracError;
      const Double_t badAvgErrB    = badAvgB * fracError;
      const Bool_t   gudAvgAisGood = (gudAvgA != -666.);
      const Bool_t   gudAvgBisGood = (gudAvgB != -666.);
      const Bool_t   badAvgAisGood = (badAvgA != -666.);
      const Bool_t   badAvgBisGood = (badAvgB != -666.);
      const UInt_t   iBinVal       = iAvg + 1;
      if (gudAvgAisGood) {
        hGudAvg[iColumn][0] -> SetBinContent(iBinVal, gudAvgA);
        hGudAvg[iColumn][0] -> SetBinError(iBinVal, gudAvgErrA);
      }
      if (gudAvgBisGood) {
        hGudAvg[iColumn][1] -> SetBinContent(iBinVal, gudAvgB);
        hGudAvg[iColumn][1] -> SetBinError(iBinVal, gudAvgErrB);
      }
      if (badAvgAisGood) {
        hBadAvg[iColumn][0] -> SetBinContent(iBinVal, badAvgA);
        hBadAvg[iColumn][0] -> SetBinError(iBinVal, badAvgErrA);
      }
      if (badAvgBisGood) {
        hBadAvg[iColumn][1] -> SetBinContent(iBinVal, badAvgB);
        hBadAvg[iColumn][1] -> SetBinError(iBinVal, badAvgErrB);
      }
    }  // end average loop
  }  // end column loop
  cout << "    Histograms filled." << endl;

  // set styles
  const UInt_t   col[nLists] = {810, 890};
  const UInt_t   mar[nLists] = {7, 4};
  const UInt_t   lin[nLists] = {1, 2};
  const UInt_t   wid[nLists] = {1, 1};
  const UInt_t   fil[nLists] = {3004, 3005};
  const Double_t lab[nLists] = {0.02, 0.02};
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    hGudVal[iColumn][0] -> SetMarkerColor(col[0]);
    hGudVal[iColumn][0] -> SetMarkerStyle(mar[0]);
    hGudVal[iColumn][0] -> SetLineColor(col[0]);
    hGudVal[iColumn][0] -> SetLineStyle(lin[0]);
    hGudVal[iColumn][0] -> SetLineWidth(wid[0]);
    hGudVal[iColumn][0] -> SetFillColor(col[0]);
    hGudVal[iColumn][0] -> SetFillStyle(fil[0]);
    hGudVal[iColumn][0] -> GetXaxis() -> SetLabelSize(lab[0]);
    hGudVal[iColumn][0] -> GetYaxis() -> SetLabelSize(lab[0]);
    hGudVal[iColumn][1] -> SetMarkerColor(col[1]);
    hGudVal[iColumn][1] -> SetMarkerStyle(mar[1]);
    hGudVal[iColumn][1] -> SetLineColor(col[1]);
    hGudVal[iColumn][1] -> SetLineStyle(lin[1]);
    hGudVal[iColumn][1] -> SetLineWidth(wid[1]);
    hGudVal[iColumn][1] -> SetFillColor(col[1]);
    hGudVal[iColumn][1] -> SetFillStyle(fil[1]);
    hGudVal[iColumn][1] -> GetXaxis() -> SetLabelSize(lab[1]);
    hGudVal[iColumn][1] -> GetYaxis() -> SetLabelSize(lab[1]);
    hGudAvg[iColumn][0] -> SetMarkerColor(col[0]);
    hGudAvg[iColumn][0] -> SetMarkerStyle(mar[0]);
    hGudAvg[iColumn][0] -> SetLineColor(col[0]);
    hGudAvg[iColumn][0] -> SetLineStyle(lin[0]);
    hGudAvg[iColumn][0] -> SetLineWidth(wid[0]);
    hGudAvg[iColumn][0] -> SetFillColor(col[0]);
    hGudAvg[iColumn][0] -> SetFillStyle(fil[0]);
    hGudAvg[iColumn][0] -> GetXaxis() -> SetLabelSize(lab[0]);
    hGudAvg[iColumn][0] -> GetYaxis() -> SetLabelSize(lab[0]);
    hGudAvg[iColumn][1] -> SetMarkerColor(col[1]);
    hGudAvg[iColumn][1] -> SetMarkerStyle(mar[1]);
    hGudAvg[iColumn][1] -> SetLineColor(col[1]);
    hGudAvg[iColumn][1] -> SetLineStyle(lin[1]);
    hGudAvg[iColumn][1] -> SetLineWidth(wid[1]);
    hGudAvg[iColumn][1] -> SetFillColor(col[1]);
    hGudAvg[iColumn][1] -> SetFillStyle(fil[1]);
    hGudAvg[iColumn][1] -> GetXaxis() -> SetLabelSize(lab[1]);
    hGudAvg[iColumn][1] -> GetYaxis() -> SetLabelSize(lab[1]);
    hBadVal[iColumn][0] -> SetMarkerColor(col[0]);
    hBadVal[iColumn][0] -> SetMarkerStyle(mar[0]);
    hBadVal[iColumn][0] -> SetLineColor(col[0]);
    hBadVal[iColumn][0] -> SetLineStyle(lin[0]);
    hBadVal[iColumn][0] -> SetLineWidth(wid[0]);
    hBadVal[iColumn][0] -> SetFillColor(col[0]);
    hBadVal[iColumn][0] -> SetFillStyle(fil[0]);
    hBadVal[iColumn][0] -> GetXaxis() -> SetLabelSize(lab[0]);
    hBadVal[iColumn][0] -> GetYaxis() -> SetLabelSize(lab[0]);
    hBadVal[iColumn][1] -> SetMarkerColor(col[1]);
    hBadVal[iColumn][1] -> SetMarkerStyle(mar[1]);
    hBadVal[iColumn][1] -> SetLineColor(col[1]);
    hBadVal[iColumn][1] -> SetLineStyle(lin[1]);
    hBadVal[iColumn][1] -> SetLineWidth(wid[1]);
    hBadVal[iColumn][1] -> SetFillColor(col[1]);
    hBadVal[iColumn][1] -> SetFillStyle(fil[1]);
    hBadVal[iColumn][1] -> GetXaxis() -> SetLabelSize(lab[1]);
    hBadVal[iColumn][1] -> GetYaxis() -> SetLabelSize(lab[1]);
    hBadAvg[iColumn][0] -> SetMarkerColor(col[0]);
    hBadAvg[iColumn][0] -> SetMarkerStyle(mar[0]);
    hBadAvg[iColumn][0] -> SetLineColor(col[0]);
    hBadAvg[iColumn][0] -> SetLineStyle(lin[0]);
    hBadAvg[iColumn][0] -> SetLineWidth(wid[0]);
    hBadAvg[iColumn][0] -> SetFillColor(col[0]);
    hBadAvg[iColumn][0] -> SetFillStyle(fil[0]);
    hBadAvg[iColumn][0] -> GetXaxis() -> SetLabelSize(lab[0]);
    hBadAvg[iColumn][0] -> GetYaxis() -> SetLabelSize(lab[0]);
    hBadAvg[iColumn][1] -> SetMarkerColor(col[1]);
    hBadAvg[iColumn][1] -> SetMarkerStyle(mar[1]);
    hBadAvg[iColumn][1] -> SetLineColor(col[1]);
    hBadAvg[iColumn][1] -> SetLineStyle(lin[1]);
    hBadAvg[iColumn][1] -> SetLineWidth(wid[1]);
    hBadAvg[iColumn][1] -> SetFillColor(col[1]);
    hBadAvg[iColumn][1] -> SetFillStyle(fil[1]);
    hBadAvg[iColumn][1] -> GetXaxis() -> SetLabelSize(lab[1]);
    hBadAvg[iColumn][1] -> GetYaxis() -> SetLabelSize(lab[1]);
  }
  cout << "    Styles set." << endl;

  // create canvases
  const UInt_t width(700);
  const UInt_t height(700);
  const UInt_t grid(0);
  TCanvas *cGudVal[nColumns];
  TCanvas *cGudAvg[nColumns];
  TCanvas *cBadVal[nColumns];
  TCanvas *cBadAvg[nColumns];
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    TString sGudVal("cGoodVal");
    TString sGudAvg("cGoodAvg");
    TString sBadVal("cBadVal");
    TString sBadAvg("cBadAvg");
    sGudVal += iColumn;
    sGudAvg += iColumn;
    sBadVal  += iColumn;
    sBadAvg  += iColumn;
    cGudVal[iColumn] = new TCanvas(sGudVal.Data(), "", width, height);
    cGudVal[iColumn]    -> cd();
    cGudVal[iColumn]    -> SetGrid(grid, grid);
    hGudVal[iColumn][0] -> Draw("E5");
    hGudVal[iColumn][1] -> Draw("E5 same");
    cGudVal[iColumn]    -> Write();
    cGudVal[iColumn]    -> Close();
    cGudAvg[iColumn] = new TCanvas(sGudAvg.Data(), "", width, height);
    cGudAvg[iColumn]    -> cd();
    cGudAvg[iColumn]    -> SetGrid(grid, grid);
    hGudAvg[iColumn][0] -> Draw("E5");
    hGudAvg[iColumn][1] -> Draw("E5 same");
    cGudAvg[iColumn]    -> Write();
    cGudAvg[iColumn]    -> Close();
    cBadVal[iColumn] = new TCanvas(sBadVal.Data(), "", width, height);
    cBadVal[iColumn]    -> cd();
    cBadVal[iColumn]    -> SetGrid(grid, grid);
    hBadVal[iColumn][0] -> Draw("E5");
    hBadVal[iColumn][1] -> Draw("E5 same");
    cBadVal[iColumn]    -> Write();
    cBadVal[iColumn]    -> Close();
    cBadAvg[iColumn] = new TCanvas(sBadAvg.Data(), "", width, height);
    cBadAvg[iColumn]    -> cd();
    cBadAvg[iColumn]    -> SetGrid(grid, grid);
    hBadAvg[iColumn][0] -> Draw("E5");
    hBadAvg[iColumn][1] -> Draw("E5 same");
    cBadAvg[iColumn]    -> Write();
    cBadAvg[iColumn]    -> Close();
  }
  cout << "    Canvases created." << endl;


  // close output file
  fOutput -> cd();
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    hGudVal[iColumn][0] -> Write();
    hGudVal[iColumn][1] -> Write();
    hGudAvg[iColumn][0] -> Write();
    hGudAvg[iColumn][1] -> Write();
    hBadVal[iColumn][0]  -> Write();
    hBadVal[iColumn][1]  -> Write();
    hBadAvg[iColumn][0]  -> Write();
    hBadAvg[iColumn][1]  -> Write();
  }
  fOutput -> Close();
  cout << "  Finished checking!\n" << endl;

}

// End ------------------------------------------------------------------------
