// 'ListChecker.C'
// Derek Anderson
// 10.19.2017
//
// This opens a pair of lists (preferably WSD) and
// checks them.  If a particular line differs between
// the 2 lists, it spits out a notification.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;


// global constants
static const UInt_t  nLines(4800);
static const UInt_t  nAverage(100);
static const UInt_t  nAvgValues(nLines / nAverage);
static const UInt_t  nColumns(7);
static const TString sListA("lists/pp200r9ofl.calibrationsPruned.d15m12y2008.txt");
static const TString sListB("lists/pp200r12ofl.calibrationsPruned.d20m12y2011.txt");
static const TString sOutput("bemcAdcCheck.m12y2008oflXm12y2011ofl.d27m10y2017.root");



UInt_t ListChecker() {

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
  TH1D *hValue[nColumns][2];
  TH1D *hAverage[nColumns][2];
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    TString sValueA("hValues");
    TString sAverageA("hAverages");
    TString sValueB("hValues");
    TString sAverageB("hAverages");
    sValueA   += iColumn;
    sAverageA += iColumn;
    sValueB   += iColumn;
    sAverageB += iColumn;
    sValueA   += "A";
    sAverageA += "A";
    sValueB   += "B";
    sAverageB += "B";
    hValue[iColumn][0]   = new TH1D(sValueA.Data(), "", nLines, 0., (Double_t) nLines);
    hValue[iColumn][1]   = new TH1D(sValueB.Data(), "", nLines, 0., (Double_t) nLines);
    hAverage[iColumn][0] = new TH1D(sAverageA.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
    hAverage[iColumn][1] = new TH1D(sAverageB.Data(), "", nAvgValues, 0., (Double_t) nAvgValues);
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
  UInt_t   nDiff(0);
  Double_t valuesA[nLines][nColumns];
  Double_t valuesB[nLines][nColumns];
  Double_t averageA[nAvgValues][nColumns];
  Double_t averageB[nAvgValues][nColumns];
  while (listA && listB) {

    // read line, check if same
    Int_t  iDiffValue(-1);
    Bool_t valueIsDifferent(false);
    for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
      listA >> valuesA[iLine][iColumn];
      listB >> valuesB[iLine][iColumn];
      valueIsDifferent = (valuesA[iLine][iColumn] != valuesB[iLine][iColumn]);
      if (valueIsDifferent) {
        iDiffValue = iColumn;
      }
      averageA[iAvg][iColumn] += valuesA[iLine][iColumn];
      averageB[iAvg][iColumn] += valuesB[iLine][iColumn];
    }  // end column loop

    // if not same, notify
    const Bool_t lineIsDifferent = (iDiffValue != -1);
    if (lineIsDifferent) {
      cout << "      Line is different!\n"
           << "        (line, column)   = (" << iLine << ", " << iDiffValue << ")\n"
           << "        (valueA, valueB) = (" << valuesA[iLine][iDiffValue] << ", " << valuesB[iLine][iDiffValue] << ")"
           << endl;
       nDiff++;
    }
    iVal++;
    iLine++;

    // check if time to exit loop
    const Bool_t loopIsFinished      = (iLine == nLines);
    const Bool_t averageIsCalculated = (iVal == nAverage);
    if (loopIsFinished) {
      for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
        averageA[iAvg][iColumn] /= iVal;
        averageB[iAvg][iColumn] /= iVal;
      }
      break;
    }
    else if (averageIsCalculated) {
      for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
        averageA[iAvg][iColumn] /= nAverage;
        averageB[iAvg][iColumn] /= nAverage;
      }
      iVal = 0;
      iAvg++;
    }

  }   // end line loop

  // close streams
  listA.close();
  listB.close();
  cout << "    Done processing.\r"
       << "      " << nDiff << " lines differed!"
       << endl;


  // fill histograms
  const Double_t fracError = 0.2;
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    for (iLine = 0; iLine < nLines; iLine++) {
      const Double_t valA    = valuesA[iLine][iColumn];
      const Double_t valB    = valuesB[iLine][iColumn];
      const Double_t valErrA = valA * fracError;
      const Double_t valErrB = valB * fracError;
      const UInt_t   iBinVal = iLine + 1;
      hValue[iColumn][0] -> SetBinContent(iBinVal, valA);
      hValue[iColumn][1] -> SetBinContent(iBinVal, valB);
      hValue[iColumn][0] -> SetBinError(iBinVal, valErrA);
      hValue[iColumn][1] -> SetBinError(iBinVal, valErrB);
    }  // end line loop
    for (iAvg = 0; iAvg < nAvgValues; iAvg++) {
      const Double_t avgA    = averageA[iAvg][iColumn];
      const Double_t avgB    = averageB[iAvg][iColumn];
      const Double_t avgErrA = avgA * fracError;
      const Double_t avgErrB = avgB * fracError;
      const UInt_t   iBinAvg = iAvg + 1;
      hAverage[iColumn][0] -> SetBinContent(iBinAvg, avgA);
      hAverage[iColumn][1] -> SetBinContent(iBinAvg, avgB);
      hAverage[iColumn][0] -> SetBinError(iBinAvg, avgErrA);
      hAverage[iColumn][1] -> SetBinError(iBinAvg, avgErrB);
    }  // end average loop
  }  // end column loop
  cout << "    Histograms filled." << endl;

  // set styles
  const UInt_t   col[2] = {810, 890};
  const UInt_t   mar[2] = {7, 4};
  const UInt_t   lin[2] = {1, 2};
  const UInt_t   wid[2] = {1, 1};
  const UInt_t   fil[2] = {3004, 3005};
  const Double_t lab[2] = {0.02, 0.02};
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    hValue[iColumn][0]   -> SetMarkerColor(col[0]);
    hValue[iColumn][0]   -> SetMarkerStyle(mar[0]);
    hValue[iColumn][0]   -> SetLineColor(col[0]);
    hValue[iColumn][0]   -> SetLineStyle(lin[0]);
    hValue[iColumn][0]   -> SetLineWidth(wid[0]);
    hValue[iColumn][0]   -> SetFillColor(col[0]);
    hValue[iColumn][0]   -> SetFillStyle(fil[0]);
    hValue[iColumn][0]   -> GetXaxis() -> SetLabelSize(lab[0]);
    hValue[iColumn][0]   -> GetYaxis() -> SetLabelSize(lab[0]);
    hValue[iColumn][1]   -> SetMarkerColor(col[1]);
    hValue[iColumn][1]   -> SetMarkerStyle(mar[1]);
    hValue[iColumn][1]   -> SetLineColor(col[1]);
    hValue[iColumn][1]   -> SetLineStyle(lin[1]);
    hValue[iColumn][1]   -> SetLineWidth(wid[1]);
    hValue[iColumn][1]   -> SetFillColor(col[1]);
    hValue[iColumn][1]   -> SetFillStyle(fil[1]);
    hValue[iColumn][1]   -> GetXaxis() -> SetLabelSize(lab[1]);
    hValue[iColumn][1]   -> GetYaxis() -> SetLabelSize(lab[1]);
    hAverage[iColumn][0] -> SetMarkerColor(col[0]);
    hAverage[iColumn][0] -> SetMarkerStyle(mar[0]);
    hAverage[iColumn][0] -> SetLineColor(col[0]);
    hAverage[iColumn][0] -> SetLineStyle(lin[0]);
    hAverage[iColumn][0] -> SetLineWidth(wid[0]);
    hAverage[iColumn][0] -> SetFillColor(col[0]);
    hAverage[iColumn][0] -> SetFillStyle(fil[0]);
    hAverage[iColumn][0] -> GetXaxis() -> SetLabelSize(lab[0]);
    hAverage[iColumn][0] -> GetYaxis() -> SetLabelSize(lab[0]);
    hAverage[iColumn][1] -> SetMarkerColor(col[1]);
    hAverage[iColumn][1] -> SetMarkerStyle(mar[1]);
    hAverage[iColumn][1] -> SetLineColor(col[1]);
    hAverage[iColumn][1] -> SetLineStyle(lin[1]);
    hAverage[iColumn][1] -> SetLineWidth(wid[1]);
    hAverage[iColumn][1] -> SetFillColor(col[1]);
    hAverage[iColumn][1] -> SetFillStyle(fil[1]);
    hAverage[iColumn][1] -> GetXaxis() -> SetLabelSize(lab[1]);
    hAverage[iColumn][1] -> GetYaxis() -> SetLabelSize(lab[1]);
  }
  cout << "    Styles set." << endl;

  // create canvases
  const UInt_t width(700);
  const UInt_t height(700);
  const UInt_t grid(0);
  TCanvas *cValues[nColumns];
  TCanvas *cAverages[nColumns];
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    TString sValue("cValues");
    TString sAverage("cAverage");
    sValue   += iColumn;
    sAverage += iColumn;
    cValues[iColumn] = new TCanvas(sValue.Data(), "", width, height);
    cValues[iColumn]   -> cd();
    cValues[iColumn]   -> SetGrid(grid, grid);
    hValue[iColumn][0] -> Draw("E5");
    hValue[iColumn][1] -> Draw("E5 same");
    cValues[iColumn]   -> Write();
    cValues[iColumn]   -> Close();
    cAverages[iColumn] = new TCanvas(sAverage.Data(), "", width, height);
    cAverages[iColumn]   -> cd();
    cAverages[iColumn]   -> SetGrid(grid, grid);
    hAverage[iColumn][0] -> Draw("E5");
    hAverage[iColumn][1] -> Draw("E5 same");
    cAverages[iColumn]   -> Write();
    cAverages[iColumn]   -> Close();
  }
  cout << "    Canvases created." << endl;


  // close output file
  fOutput -> cd();
  for (UInt_t iColumn = 0; iColumn < nColumns; iColumn++) {
    hValue[iColumn][0]   -> Write();
    hValue[iColumn][1]   -> Write();
    hAverage[iColumn][0] -> Write();
    hAverage[iColumn][1] -> Write();
  }
  fOutput -> Close();
  cout << "  Finished checking!\n" << endl;
  return nDiff;

}

// End ------------------------------------------------------------------------
