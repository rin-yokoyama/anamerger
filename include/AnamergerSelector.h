// AnamergerSelector for PROOF created by Rin Yokoyama on 7/21/2017

#ifndef ANAMERGER_SELECTOR_H
#define ANAMERGER_SELECTOR_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TSelector.h"
#include "TProofServ.h"
#include "TMath.h"
#include "TParameter.h"
#include "data2Tree.hpp"
#include "hIsotope.h"

/**
 * @brief TSelector class for the anamerger main analysis
 * 
 */
class AnamergerSelector : public TSelector {
public:

	/**
	 * @brief Construct a new Anamerger Selector object
	 * 
	 */
	AnamergerSelector(TTree* = 0);

	/**
	 * @brief Destroy the Anamerger Selector object
	 * 
	 */
	virtual ~AnamergerSelector();

	/// Overloaded from TSelector
	virtual Int_t   Version() const { return 1; }
	virtual void    Init(TTree* mergedData);
	virtual void    Begin(TTree* mergedData);
	virtual void    SlaveBegin(TTree* mergedData);
	virtual Bool_t  Notify() { return kTRUE; }
	virtual Bool_t  Process(Long64_t entry);
	virtual void    SetOption(const char* option) { fOption = option; }
	virtual void    SetObject(TObject* obj) { fObject = obj; }
	virtual void    SetInputList(TList* input) { fInput = input; }
	virtual TList* GetOutputList() const { return fOutput; }
	virtual void    SlaveTerminate() { tree_reader_.SetTree((TTree*)nullptr); }
	virtual void    Terminate();

	/**
	 * @brief Set the Hist Group object
	 * 
	 * @param group_name 
	 */
	void SetHistGroup(const std::string &group_name) { hist_group_map_[group_name] = true; }

	/**
	 * @brief Set the Reference Cut Name
	 * 
	 * @param name name of the reference cut file
	 */
	void SetReferenceCutName(const std::string &name) { ref_cut_name_ = name; }

	/**
	 * @brief Set the Output File Name
	 * 
	 * @param file_name output file name
	 */
	void SetOutputFileName(const std::string& file_name) {
		output_file_name_ = file_name;
	}

	/**
	 * @brief Set the Correlation Radius
	 * 
	 * @param radius correlation raidus
	 */
	void SetCorrelationRadius(const double &radius) { correlation_radius_ = radius; }

	/**
	 * @brief Set the gamma-gamma tdiff
	 * 
	 * @param tdiff time window of the gamma-gamma correlation
	 */
	void SetGGTdiff(const double &tdiff) { gg_tdiff_ = tdiff; }

protected:
	/**
	 * @brief creates isotope cut
	 * 
	 * @param icutname name of the reference cut file name
	 * @return int 0: success, 1: failed
	 */
	int loadCUTG(std::string icutname);

	std::vector<hIsotope> vectorIsotopes; /// vector of hIsotope objects

	TTreeReader tree_reader_;
	TTreeReaderValue <brData2TTree>    bigrips;
	TTreeReaderValue <impData2TTree>   implant;
	TTreeReaderValue <betaData2TTree>  beta;
	TTreeReaderValue <neuData2TTree>   neutron;
	TTreeReaderValue <gammaData2TTree> gamma;
	TTreeReaderValue <ancData2TTree>   ancillary;
	ULong64_t total_entry_; /// Event counter

	TObjArray* fHistArray = nullptr; /// a pointer to the histogram array
	TFile* fOutputFile = nullptr; /// a pointer to the output file object
	std::string output_file_name_; /// output file name
	std::string ref_cut_name_; /// reference cut file name
	std::map<std::string, bool> hist_group_map_; /// a map of histogram groups
	double correlation_radius_; /// wasabi correlation distance (pixel)
	double gg_tdiff_; /// gamma gamma Tib time window (sec)

	ClassDef(AnamergerSelector, 1)
};

#endif // End of ANAMERGER_SELECTOR_H
