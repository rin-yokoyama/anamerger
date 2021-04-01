// AnamergerSelector for PROOF created by Rin Yokoyama on 7/21/2017

#ifndef ANAMERGER_SELECTOR_H
#define ANAMERGER_SELECTOR_H

#include <iostream>
#include <vector>
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

class AnamergerSelector : public TSelector {
public:

	AnamergerSelector(TTree* = 0);
	virtual ~AnamergerSelector();

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
	void SetReferenceCutName(const std::string &name) { ref_cut_name_ = name; }
	void SetOutputFileName(const std::string& file_name) {
		output_file_name_ = file_name;
	}

protected:
	int loadCUTG(std::string icutname);
	std::vector<hIsotope> vectorIsotopes;

	TTreeReader tree_reader_;
	TTreeReaderValue <brData2TTree>    bigrips;
	TTreeReaderValue <impData2TTree>   implant;
	TTreeReaderValue <betaData2TTree>  beta;
	TTreeReaderValue <neuData2TTree>   neutron;
	TTreeReaderValue <gammaData2TTree> gamma;
	TTreeReaderValue <ancData2TTree>   ancillary;
	ULong64_t total_entry_;

	// array for histograms
	TObjArray* fHistArray = nullptr;
	// output file
	TFile* fOutputFile = nullptr;
	std::string output_file_name_;
	std::string ref_cut_name_;

	ClassDef(AnamergerSelector, 1)
};

#endif
