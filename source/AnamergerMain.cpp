#include <iostream>
#include <fstream>
#include <string>
#include "TChain.h"
#include "TProof.h"
#include "TROOT.h"
#include "YamlParameter.hpp"
#include "YamlReader.hpp"
#include "AnamergerSelector.h"
#include "LibraryConfig.h"

void usage(char const* arg)
{
	std::cout << arg << "-c [yaml_file]" << std::endl;
}

int main(int argc, char** argv)
{
	std::string input_file_name;
	std::string config_file_name;
	std::string output_file_name = "anamerger_output.root";

	/** parsing commandline arguments **/
	if (argc < 3) {
		usage(argv[0]);
		return 1;
	}
	int opt = 0;
	while ((opt = getopt(argc, argv, "i:o:c:")) != -1) {
		switch (opt) {
		case 'i':
			input_file_name = optarg;
			break;
		case 'o':
			output_file_name = optarg;
			break;
		case 'c':
			config_file_name = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
			break;
		}
	}

	gROOT->SetBatch();
	TProof* pr;

	try {
		/** Create YamlParameter instance **/
		YamlParameter::Create(config_file_name);
		YamlReader* yaml_reader_ = new YamlReader("AnamergerMain");
		/** Read run configurations from the yaml file**/
		const std::string tree_name = yaml_reader_->GetString("TreeName");
		const std::string merger_list_name = yaml_reader_->GetString("MergerListName");
		const std::string ref_cut_name = yaml_reader_->GetString("ReferenceCutsName");
		const bool use_proof = yaml_reader_->GetBoolean("UseProof", false, false);
		if (use_proof) {
			const std::string num_workers = yaml_reader_->GetString("NumWorkers");
			const std::string proof_arg = "workers=" + num_workers;
			/** Open a proof lite server **/
			pr = TProof::Open("lite://", proof_arg.c_str());
			/** Absolute paths for the required libraries **/
			std::vector<std::string> libs =
			{
			  getYamlcppLibDir() + "/libyaml-cpp.so",
			  getMergerLibDir() + "/libAnamergerLib.so",
			};
			/** Load libraries to the proof server **/
			for (const auto& lib : libs) {
				pr->Load(lib.c_str());
			}
		}
		/** Add files to the input TChain **/
		TChain* chain = new TChain(tree_name.c_str());
		std::ifstream mergerListFiles(merger_list_name.c_str());
		std::string tempName("");

		while (getline(mergerListFiles, tempName))
		{
			chain->AddFile(tempName.c_str());
			std::cout << tempName.c_str() << std::endl;
		}

		mergerListFiles.close();
		std::cout << "[AnamergerMain]: Number of TTrees added to the chain: " << chain->GetNtrees() << std::endl;

		/** Number of entries to scan and the first entry to start scanning **/
		const unsigned long long n_entries = yaml_reader_->GetULong64("NumEntries", false, chain->GetEntries());
		const unsigned long long first_entry = yaml_reader_->GetULong64("FirstEntry", false, 0);

		if (use_proof) {
			chain->SetProof();
			std::cout << "SetProof to the chain: " << chain->GetName() << std::endl;
			/** Add parameters to the proof server **/
			pr->AddInput(new TNamed("output_file_name", output_file_name.c_str()));
			pr->AddInput(new TNamed("ref_cut_name", ref_cut_name.c_str()));
			chain->Process("AnamergerSelector", "", n_entries, first_entry);
		}
		else {
			std::cout << "Start Processing (Proof OFF)..." << std::endl;
			AnamergerSelector* selector = new AnamergerSelector(chain);
			selector->SetOutputFileName(output_file_name);
			selector->SetReferenceCutName(ref_cut_name);
			chain->Process(selector, "", n_entries, first_entry);
		}

		if (use_proof)
			pr->Close();

		/** destroys YamlParameter instance **/
		YamlParameter::Destroy();
	}
	catch (std::string msg) {
		std::cout << msg << std::endl;
		std::cout << "[AnamergerMain]: exiting from main() due to the error" << std::endl;
		return 1;
	}
	catch (std::bad_alloc) {
		std::cout << "[AnamergerMain]: bad_alloc occured while setting up." << std::endl;
		std::cout << "[AnamergerMain]: exiting from main() due to the error" << std::endl;
		return 1;
	}

	return 0;
}
