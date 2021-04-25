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
#include <map>

/// Watch out the time-offsets!

//Beta-neutron time correlation cuts (DT=200us)
const Long_t DTbnlow=-20e4;
const Long_t DTbnhigh=20e4;

//Beta-gamma time correlation cuts (DT=7.5us) 
const Long_t DTgblow=1.25e4;
const Long_t DTgbhigh=2e4;

//Beta-implant time gate for gammas (DT=5s)
const Long_t DTibhigh=5e9;

//Beta-F11 time cut window (DT=100us)
const Long_t DTfblow=5000;
const Long_t DTfbhigh=150000;

//Neutron-F11 time cut window (DT=250us)
const Long_t DTfnlow=-50e3;
const Long_t DTfnhigh=2e5;

//Neutron energy cut
const Double_t Enlow=175;
const Double_t Enhigh=1000;
const Double_t Eplow=1100;

//Beta energy cut
const Double_t Eblow=100;
   
/**
 * @brief hIsotope manages all the histograms related to an isotope
 * 
 * This class defines all isotope-by-isotope histograms.
 * It also provides the PID gate of the isotope which can be
 * defined either as a ROOT TCutG object or an ellipse with
 * four parameters, a, b, x0, y0 as in the formula
 * a*(x-x0)^2 + b*(y-y0)^2 < 1
 */
class hIsotope
{
  public:
    /**
     * @brief Default constructor
     * 
     */
    hIsotope(){}
    
    /**
     * @brief Construct a new h Isotope object with a ROOT TCutG
     * 
     * @param aCut a pointer to the TCutG object
     * @param outputList a pointer to the fOutputList of the TSelector
     * @param hist_group_map a map of histogram group names
     */
    hIsotope( TCutG * aCut, TList* outputList, const std::map<std::string, bool> hist_group_map );

    /**
     * @brief Construct a new h Isotope object with an ellipse cut
     * 
     * @param isoname name of the isotope
     * @param a  1/r^2 of the x-radius of an ellipse as in the formula, a*(x-x0)^2 + b*(y-y0)^2 < 1
     * @param b  1/r^2 of the y-radius of an ellipse as in the formula, a*(x-x0)^2 + b*(y-y0)^2 < 1
     * @param x0  x-center of an ellipse as in the formula, a*(x-x0)^2 + b*(y-y0)^2 < 1
     * @param y0  y-center of an ellipse as in the formula, a*(x-x0)^2 + b*(y-y0)^2 < 1
     * @param outputList a pointer to the fOutputList of the TSelector
     * @param hist_group_map a map of histogram group names
     */
    hIsotope( std::string isoname, Double_t a, Double_t b, Double_t x0, Double_t y0, TList* outputList, const std::map<std::string, bool> hist_group_map );

    /**
     * @brief Destroy the h Isotope object
     * 
     */
    ~hIsotope();

    /**
     * @brief Initialization with a TCutG object
     * 
     * @param aCut a pointer to the TCutG object
     */
    void initializeCut(TCutG* aCut);

    /**
     * @brief Initialization of histograms
     * 
     * @param outputList a pointer to the fOutputList
     */
    void initializeHistos(TList* outputList);

    /**
     * @brief checks if the given x, y is inside the cut
     * 
     * @param x AoQ
     * @param y Z
     * @return true if (x,y) is inside the TCutG or the ellipse
     * @return false otherwise
     */
    bool IsInside(Double_t x, Double_t y) const;

    int32_t numberOfImplants;
    int32_t numberOfPID;
    int32_t numberOfCorrelatedImp;
    TCutG * myCutG; /// Cut for PID
    
    Double_t ellipseA;  /// 1/r_x^2 of the PID ellipse
    Double_t ellipseB;  /// 1/r_y^2 of the PID ellipse
    Double_t ellipseX0;  /// x center of the PID ellipse
    Double_t ellipseY0;  /// y center of the PID ellipse

    TObjArray *fHistArray; /// pointer to the histogram array
    std::string isotopeName; /// name of the isotope
    std::map<std::string, bool> hist_group_map_; /// map of the histogram group names
};


#endif // end of __h_isotopeFastID_hpp__
