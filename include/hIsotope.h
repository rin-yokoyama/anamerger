// // // // // // // // // // // // // // // // // // 
// 
// Definition of an structure which gathers
// all the histograms associated with ONE isotope
// 
// VERSION 2.0. NEW HISTOGRAMS ADDED!!
// 
// // // // // // // // // // // // // // // // // // 


#ifndef __h_isotopeFastID_hpp__
#define __h_isotopeFastID_hpp__

#include "TFile.h"
#include "TCutG.h"
#include "TList.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TKey.h"
#include "TIterator.h"
#include "TMath.h"
#include "data2Tree.hpp"

#include <string>
#include <iostream>
#include <vector>

//#define BETA
//#define DECAYCURVE
//#define GAMMA_ET
//#define GAMMA_GAMMA
#define ISOMER

/// Watch out the time-offsets!

//Beta-neutron time correlation cuts (DT=200us)
const Long_t DTbnlow=-2e4;
const Long_t DTbnhigh=18e4;

//Beta-gamma time correlation cuts (DT=7.5us) 
const Long_t DTgblow=1.25e4;
const Long_t DTgbhigh=2e4;

//Beta-implant time gate for gammas (DT=5s)
const Long_t DTibhigh=5e9;

//Beta-F11 time cut window (DT=100us)
const Long_t DTfblow=-30000;
const Long_t DTfbhigh=70000;

//Neutron-F11 time cut window (DT=250us)
const Long_t DTfnlow=-50e3;
const Long_t DTfnhigh=2e5;

//Neutron energy cut
const Double_t Enlow=175;
const Double_t Enhigh=850;
const Double_t Eplow=1100;

//Beta-implant max. X,Y distance
const Double_t DSib=3.5;
 
//Beta energy cut
const Double_t Eblow=0;
   
class hIsotope
{
  public:
    hIsotope(){}
    // Constructor with ROOT TCutG
    hIsotope( TCutG * aCut, TList* outputList );
    // Constructor with parameters for ellipse cut.
    hIsotope( std::string isoname, Double_t a, Double_t b, Double_t x0, Double_t y0, TList* outputList );
    ~hIsotope();
    void initializeCut(TCutG* aCut);
    void initializeHistos(TList* outputList);
    bool IsInside(Double_t x, Double_t y) const;

    /// double numberOfImplants
    int32_t numberOfImplants;
    /// double numberOfImplants
    int32_t numberOfPID;
    /// number of correlated implants in an event
    int32_t numberOfCorrelatedImp;
    /// myCutG: cut for PID
    TCutG * myCutG;
    Double_t ellipseA;
    Double_t ellipseB;
    Double_t ellipseX0;
    Double_t ellipseY0;

    TObjArray *fHistArray;
    std::string isotopeName;
};


#endif //

