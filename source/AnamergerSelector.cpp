#include "AnamergerSelector.h"
#include <fstream>
#include <iostream>
#include <string>
#include <TList.h>

ClassImp(AnamergerSelector);

AnamergerSelector::AnamergerSelector(TTree *mergedData) : tree_reader_(mergedData),
														  bigrips(tree_reader_, "bigrips."),
														  implant(tree_reader_, "implantation."),
														  beta(tree_reader_, "beta."),
														  neutron(tree_reader_, "neutron."),
														  gamma(tree_reader_, "gamma."),
														  ancillary(tree_reader_, "ancillary."),
														  output_file_name_("anamerger_output.root")
{
}

AnamergerSelector::~AnamergerSelector()
{
	if (fHistArray)
	{
		delete fHistArray;
		fHistArray = nullptr;
	}
	if (fOutputFile)
	{
		delete fOutputFile;
		fOutputFile = nullptr;
	}
}

void AnamergerSelector::Begin(TTree *mergedData)
{
	GetOutputList()->Clear();
	if (fInput)
	{
		TNamed *named = (TNamed *)fInput->FindObject("output_file_name");
		if (named)
			output_file_name_ = named->GetTitle();
	}
}

void AnamergerSelector::SlaveBegin(TTree *mergedData)
{
	if (fHistArray)
	{
		delete fHistArray;
		fHistArray = nullptr;
	}
	fHistArray = new TObjArray();

	if (fInput)
	{
		TNamed *named = (TNamed *)fInput->FindObject("ref_cut_name");
		if (named)
			ref_cut_name_ = named->GetTitle();
		auto param = (TParameter<Double_t>*)fInput->FindObject("correlation_radius");
		if (param)
			correlation_radius_ = param->GetVal();
		param = (TParameter<Double_t>*)fInput->FindObject("gg_tdiff");
		if (param)
			gg_tdiff_ = param->GetVal();


		TIter next(fInput);
		while (auto *obj = (TObject *)next())
		{
			const std::string name(obj->GetName());
			if (name == "output_file_name" || name == "reference_cut_name" || name == "gg_tdiff" || name == "correlation_radius")
				continue;
			auto n = (TNamed *)obj;
			const std::string title(n->GetTitle());
			if (title == "True" || title == "true")
				hist_group_map_[n->GetName()] = true;
		}
	}

	fHistArray->Add(new TH1F("countPID", "countPID", 100, 0, 100));
	if (hist_group_map_.count("TEST"))
	{
		fHistArray->Add(new TH1F("betaEx", "betaEx", 1000, 0, 1000));
	}

	//adding histograms to output list
	TIter next(fHistArray);
	while (auto hist = next())
	{
		GetOutputList()->Add(hist);
	}

	loadCUTG(ref_cut_name_);

	if (gProofServ)
	{
		const TString msg = TString::Format("SlaveBegin() of Ord = %s called. %d histograms are initialized.",
											gProofServ->GetOrdinal(), GetOutputList()->GetEntries());
		gProofServ->SendAsynMessage(msg);
	}
	else
	{
		std::cout << "SalveBegin() called. " << GetOutputList()->GetEntries() << " histograms are initialized." << std::endl;
	}

	return;
}

void AnamergerSelector::Init(TTree *mergedData)
{
	tree_reader_.SetTree(mergedData);
	return;
}

Bool_t AnamergerSelector::Process(Long64_t entry)
{
	tree_reader_.SetLocalEntry(entry);

	// test
	if (hist_group_map_.count("TEST"))
	{
		((TH1F *)fHistArray->FindObject("betaEx"))->Fill((*beta).Ex);
	}
	// Implant events
	if (!(*implant).vectorOfPid.empty())
	{
		Int_t i = 0;
		for (const auto &vit : vectorIsotopes)
		{
			i++;
			if (vit.IsInside((*implant).aoq, (*implant).zet))
			{
				break;
			}
		}
		((TH1F *)fHistArray->FindObject("countPID"))->Fill(i);

		if (hist_group_map_.count("ISOMER"))
		{
			for (auto vit : vectorIsotopes)
			{
				if (vit.IsInside((*implant).aoq, (*implant).zet))
				{
					std::vector<GammaData> d4_vec;
					std::vector<GammaData> g7_vec;
					for (auto gm : (*implant).vectorOfGamma)
					{
						const Double_t tdiff = (gm.TIME - (*implant).T) * 1E-9;
						((TH2F *)vit.fHistArray->FindObject(std::string("hETgI" + vit.isotopeName).c_str()))->Fill(gm.EN, tdiff);

						if (tdiff > 1e-6 && tdiff < 5e-6)
						{
							if (gm.INDEX1 == 1)
								d4_vec.push_back(gm);
							if (gm.INDEX1 == 2)
								g7_vec.push_back(gm);
						}
					}

					for (auto d4 : d4_vec)
					{
						for (auto g7 : g7_vec)
						{
							const Double_t tdiff = d4.TIME - g7.TIME;
							if (tdiff > -100. && tdiff < 100.)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEEgI" + vit.isotopeName).c_str()))->Fill(d4.EN, g7.EN);
						}
					}
				}
			}
		}
	}

	if (hist_group_map_.count("BETA"))
	{
		//if ((*beta).z > 3) // WAS3ABi
		//	return kTRUE;

		//for (const auto &vbeta : (*beta).vectorOfBeta)
		//{
		//	const Double_t tdiff = ((double)(*beta).T - (double)vbeta.TIME) * 1.E-9;
		//	if (tdiff > 0 && tdiff < 50.E-6)
		//		return kTRUE;
		//}

		for (const auto &imp : (*beta).vectorOfImp)
		{
			//if (imp.Z > 3) // WAS3ABi
			//	continue;
			//bool yso = false;
			//for (const auto &imp2 : (*beta).vectorOfImp)
			//{
			//	const Double_t tdiff = (imp.TIME - imp2.TIME) * 1.E-9;
			//	if (imp2.Z > 3 && tdiff > -10E-6 && tdiff < 10E-6)
			//		yso = true;
			//}
			//	continue;
			if ((*beta).z != imp.Z)
				continue;
			{
				const Double_t pdist = correlation_radius_ * correlation_radius_; // WAS3ABi
				if (TMath::Power((*beta).x - imp.X, 2) + TMath::Power((*beta).y - imp.Y, 2) > pdist)
					continue;
			}
			Int_t nmult = 0;
			for (const auto &neu : (*beta).vectorOfNeu)
			{
				Double_t tdiff = neu.TIME - (*beta).T;
				if (DTbnlow < tdiff && tdiff < DTbnhigh && neu.EN > Enlow && neu.EN < Enhigh)
					nmult++;
			}
			for (const auto &vit : vectorIsotopes)
			{
				if (vit.IsInside(imp.AOQ, imp.ZET))
				{

					const Double_t tdiff = ((*beta).T - imp.TIME) * 1.E-9;

					if (hist_group_map_.count("DECAYCURVE"))
					{
						((TH1F *)vit.fHistArray->FindObject(std::string("hTib" + vit.isotopeName).c_str()))->Fill(tdiff);
						if (!nmult)
							((TH1F *)vit.fHistArray->FindObject(std::string("hTib0n" + vit.isotopeName).c_str()))->Fill(tdiff);
						if (nmult == 1)
							((TH1F *)vit.fHistArray->FindObject(std::string("hTibn" + vit.isotopeName).c_str()))->Fill(tdiff);
						if (nmult == 2)
							((TH1F *)vit.fHistArray->FindObject(std::string("hTibnn" + vit.isotopeName).c_str()))->Fill(tdiff);
						if (nmult == 3)
							((TH1F *)vit.fHistArray->FindObject(std::string("hTib3n" + vit.isotopeName).c_str()))->Fill(tdiff);
					}

					std::vector<GammaData> d4_vec;
					std::vector<GammaData> g7_vec;
					std::vector<GammaData> d4gg_vec;
					std::vector<GammaData> g7gg_vec;
					for (auto gm : (*beta).vectorOfGamma)
					{
						const Double_t bg_tdiff = ((*beta).T - gm.TIME);

						if (hist_group_map_.count("GAMMA_ET"))
						{
							if (!nmult)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTbg0n" + vit.isotopeName).c_str()))->Fill(gm.EN, bg_tdiff);
							if (nmult == 1)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTbgn" + vit.isotopeName).c_str()))->Fill(gm.EN, bg_tdiff);
							if (nmult == 2)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTbgnn" + vit.isotopeName).c_str()))->Fill(gm.EN, bg_tdiff);
						}

						//if (bg_tdiff < 2000. || bg_tdiff > 5000.)
						//	continue;
						if (gm.INDEX1 == 1)
							d4_vec.push_back(gm);
						if (gm.INDEX1 == 2)
							g7_vec.push_back(gm);

						if (tdiff > 0. && tdiff < gg_tdiff_)
						{
							if (gm.INDEX1 == 1)
								d4gg_vec.push_back(gm);
							if (gm.INDEX1 == 2)
								g7gg_vec.push_back(gm);
						}
					}

					// energy sum of four clover leaves
					Double_t esum_d4 = 0;
					for (const auto &d4 : d4_vec)
					{
						esum_d4 += d4.EN;
					}

					if (hist_group_map_.count("GAMMA_ET"))
					{
						if (esum_d4 > 0)
						{
							if (!nmult)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTib0n" + vit.isotopeName).c_str()))->Fill(esum_d4, tdiff);
							if (nmult == 1)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTibn" + vit.isotopeName).c_str()))->Fill(esum_d4, tdiff);
							if (nmult == 2)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTibnn" + vit.isotopeName).c_str()))->Fill(esum_d4, tdiff);
						}
					}

					Double_t esum_g7 = 0;
					for (const auto &g7 : g7_vec)
					{
						esum_g7 += g7.EN;
					}

					if (hist_group_map_.count("GAMMA_ET"))
					{
						if (esum_g7 > 0)
						{
							if (!nmult)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTib0n" + vit.isotopeName).c_str()))->Fill(esum_g7, tdiff);
							if (nmult == 1)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTibn" + vit.isotopeName).c_str()))->Fill(esum_g7, tdiff);
							if (nmult == 2)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEgTibnn" + vit.isotopeName).c_str()))->Fill(esum_g7, tdiff);
						}
					}

					if (hist_group_map_.count("GAMMA_GAMMA"))
					{
						if (tdiff > 0 && tdiff < 0.5)
						{
							if (!nmult)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEEg0n" + vit.isotopeName).c_str()))->Fill(esum_d4, esum_g7);
							if (nmult == 1)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEEgn" + vit.isotopeName).c_str()))->Fill(esum_d4, esum_g7);
							if (nmult == 2)
								((TH2F *)vit.fHistArray->FindObject(std::string("hEEgnn" + vit.isotopeName).c_str()))->Fill(esum_d4, esum_g7);
						}
					}
				}
			}
		}
	} //end loop through the mergedData TTree

	return kTRUE;
}

void AnamergerSelector::Terminate()
{
	fOutputFile = new TFile(output_file_name_.c_str(), "recreate");
	std::cout << "[AnamergerSelector::Terminate()]: output file: " << output_file_name_ << std::endl;
	// write the histograms
	TIter next(GetOutputList());
	while (TObject *obj = next())
	{
		std::cout << "[AnamergerSelector::Terminate]: writing " << obj->GetName() << " to file." << std::endl;
		obj->Write();
	}

	// write the parameters
	{
		std::string parString("Beta-neutron time correlation cuts. Long_t DTbnlow= ");
		parString += std::to_string(DTbnlow);
		parString += " ns. Long_t DTbnhigh= ";
		parString += std::to_string(DTbnhigh);
		parString += " ns.";
		TNamed parTNamed("par1", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	{
		std::string parString("Beta-gamma time correlation cuts. Long_t DTgblow= ");
		parString += std::to_string(DTgblow);
		parString += " ns. Long_t DTgbhigh= ";
		parString += std::to_string(DTgbhigh);
		parString += " ns.";
		TNamed parTNamed("par2", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	{
		std::string parString("Beta-implant time gate for gammas. Long_t DTibhigh= ");
		parString += std::to_string(DTibhigh);
		parString += " ns.";
		TNamed parTNamed("par3", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	{
		std::string parString("Beta-F11 time cut window. Long_t DTfblow= ");
		parString += std::to_string(DTfblow);
		parString += " ns. Long_t DTfbhigh= ";
		parString += std::to_string(DTfbhigh);
		parString += " ns.";
		TNamed parTNamed("par4", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	{
		std::string parString("Neutron-F11 time cut window. Long_t DTfnhigh= ");
		parString += std::to_string(DTfnhigh);
		parString += " ns.";
		TNamed parTNamed("par5", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	{
		std::string parString("Neutron energy cuts. Enlow= ");
		parString += std::to_string(Enlow);
		parString += ". keV Enhigh= ";
		parString += std::to_string(Enhigh);
		parString += ". keV. Eplow= ";
		parString += std::to_string(Eplow);
		parString += ". keV.";
		TNamed parTNamed("par6", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	//{
	//	std::string parString("Beta-implant max. X,Y distance. Double_t DSib= ");
	//	parString += std::to_string(DSib);
	//	parString += " pixels (square 2*DSib x 2*DSib).";
	//	TNamed parTNamed("par7", parString.c_str());
	//	parTNamed.Write(0, 2, 0);
	//}
	{
		std::string parString("Beta energy cut. Double_t Eblow= ");
		parString += std::to_string(Eblow);
		parString += " keV.";
		TNamed parTNamed("par7", parString.c_str());
		parTNamed.Write(0, 2, 0);
	}
	fOutputFile->Close();

	if (fOutputFile)
	{
		delete fOutputFile;
		fOutputFile = nullptr;
	}

	return;
}

int AnamergerSelector::loadCUTG(std::string icutname)
{
	if (icutname.find(".root") == std::string::npos)
	{
		// If the file name extension is not .root, assume a text file and read
		// paremeters for ellipse cuts.
		std::ifstream fcut(icutname);
		if (!fcut)
		{
			std::cerr << "File " << icutname.c_str() << " not found. Selector empty.\n";
			return 1;
		}
		while (!fcut.eof())
		{
			std::string isoname;
			fcut >> isoname;
			Double_t ellipse_a, ellipse_b, ellipse_x0, ellipse_y0;
			fcut >> ellipse_a;
			fcut >> ellipse_b;
			fcut >> ellipse_x0;
			fcut >> ellipse_y0;
      if (fcut.eof())
        break;
			vectorIsotopes.push_back(hIsotope(isoname, ellipse_a, ellipse_b, ellipse_x0, ellipse_y0, GetOutputList(), hist_group_map_));
      std::cout << "hIsotope " << isoname << " is created." << std::endl;
		}
		fcut.close();
	}
	else
	{
		TFile *const fcut = new TFile(icutname.c_str());

		if (!fcut)
		{
			std::cerr << "File " << icutname.c_str() << " not found. Selector empty.\n";
			return 1;
		}
		TKey *key = 0;
		TIter keyNext(fcut->GetListOfKeys());
		TCutG *temp = nullptr;
		while ((key = (TKey *)keyNext()))
		{
			std::string tempName(key->GetName());

			if (tempName.find("CUTG") != std::string::npos || tempName.find("cutg") != std::string::npos)
			{
				fcut->GetObject(key->GetName(), temp);

				vectorIsotopes.push_back(hIsotope(temp, GetOutputList(), hist_group_map_));
			}
		}
		fcut->Close();
		delete fcut;
	}
	return 0;
}
