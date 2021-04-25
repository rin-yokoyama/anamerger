// // // // // // // // // // // // // // // // // //
//
// Definition of an structure which gathers
// all the histograms associated with ONE isotope
//
// VERSION 2.0. NEW HISTOGRAMS ADDED!!
//
// // // // // // // // // // // // // // // // // //

#include "hIsotope.h"

hIsotope::hIsotope(TCutG *aCut, TList *outputList, const std::map<std::string, bool> hist_group_map) : numberOfImplants(0),
                                                                                                       numberOfPID(0),
                                                                                                       isotopeName(""),
                                                                                                       hist_group_map_(hist_group_map)
{
  myCutG = nullptr;
  isotopeName = "";
  ellipseA = 0;
  ellipseB = 0;
  ellipseX0 = 0;
  ellipseY0 = 0;

  initializeCut(aCut);
  initializeHistos(outputList);
}

hIsotope::hIsotope(std::string isoname, Double_t a, Double_t b, Double_t x0, Double_t y0, TList *outputList, const std::map<std::string, bool> hist_group_map) : numberOfImplants(0),
                                                                                                                                                                 numberOfPID(0),
                                                                                                                                                                 isotopeName(""),
                                                                                                                                                                 hist_group_map_(hist_group_map)
{
  myCutG = nullptr;
  isotopeName = isoname;
  ellipseA = a;
  ellipseB = b;
  ellipseX0 = x0;
  ellipseY0 = y0;

  initializeHistos(outputList);
}

hIsotope::~hIsotope()
{
  //  if(fHistArray)
  //    delete fHistArray;
}

void hIsotope::initializeCut(TCutG *aCut)
{
  std::string tempName4Cut(aCut->GetName());
  tempName4Cut = "i" + tempName4Cut;
  myCutG = (TCutG *)(aCut->Clone(tempName4Cut.c_str()));
  return;
}

void hIsotope::initializeHistos(TList *outputList)
{
  if (myCutG)
  {
    // get the name of the isotope
    // it must be at the end
    std::string tempName4Cut(myCutG->GetName());
    size_t lastindex = tempName4Cut.find_last_of("CUTG");
    if (lastindex == size_t(-1))
      lastindex = tempName4Cut.find_last_of("cutg");
    // the name of TCutG could be "blablaCUTG82ga"
    // and the name of the isotope will be "82ga"
    std::cout << lastindex << " " << (tempName4Cut.substr(3 + 1, tempName4Cut.size())).c_str() << std::endl;
    isotopeName = tempName4Cut.substr(4 + 1, tempName4Cut.size());
  }

  try
  {
    fHistArray = new TObjArray();

    if (hist_group_map_.count("PID"))
    {
      {
        /// PID plot
        std::string histName = "hPID" + isotopeName;
        TH2F *h = new TH2F(histName.c_str(), "PID plot", 200, 2.5, 2.8, 200, 50, 60);
        h->SetXTitle("A/Q");
        h->SetYTitle("Z");
        fHistArray->Add(h);
      }
    }
    /// Histogram definitions
    if (hist_group_map_.count("ISOMER"))
    {
      {
        /// E_gamma vs (T_gamma - T_implant)
        std::string histName = "hETgI" + isotopeName;
        TH2F *h = new TH2F(histName.c_str(), "ET plot for isomers", 4000, 0, 4000, 2000, -100E-6, 100E-6);
        h->SetXTitle("Gamma-ray energy (keV)");
        h->SetYTitle("Tgamma - Timplant (s)");
        fHistArray->Add(h);
      }
      {
        /// E_gamma vs E_gamma
        std::string histName = "hEEgI" + isotopeName;
        TH2F *h = new TH2F(histName.c_str(), "gamma-gamma plot for isomer (1 - 5 us)", 4000, 0, 4000, 4000, 0, 4000);
        h->SetXTitle("D4 energy");
        h->SetYTitle("G7 energy");
        fHistArray->Add(h);
      }
    }

    if (hist_group_map_.count("BETA"))
    {
      if (hist_group_map_.count("DECAYCURVE"))
      {
        {
          /// Total decay curve
          std::string histName = "hTib" + isotopeName;
          fHistArray->Add(new TH1F(histName.c_str(), "decay curve (total)", 2000, -10, 10));
        }
        {
          /// 0n decay curve
          std::string histName = "hTib0n" + isotopeName;
          fHistArray->Add(new TH1F(histName.c_str(), "decay curve (0n)", 2000, -10, 10));
        }
        {
          /// 1n decay curve
          std::string histName = "hTibn" + isotopeName;
          fHistArray->Add(new TH1F(histName.c_str(), "decay curve (1n)", 2000, -10, 10));
        }
        {
          /// 2n decay curve
          std::string histName = "hTibnn" + isotopeName;
          fHistArray->Add(new TH1F(histName.c_str(), "decay curve (2n)", 2000, -10, 10));
        }
        {
          /// 3n decay curve
          std::string histName = "hTib3n" + isotopeName;
          fHistArray->Add(new TH1F(histName.c_str(), "decay curve (3n)", 2000, -10, 10));
        }
      }
      if (hist_group_map_.count("GAMMA_ET"))
      {
        // E_gamma vs (T_beta - T_implant)
        {
          std::string histName = "hEgTib0n" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma vs Tbeta-Timp (0n)", 2000, 0, 2000, 2000, -10, 10));
        }
        {
          std::string histName = "hEgTibn" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma vs Tbeta-Timp (1n)", 2000, 0, 2000, 2000, -10, 10));
        }
        {
          std::string histName = "hEgTibnn" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma vs Tbeta-Timp (2n)", 2000, 0, 2000, 2000, -10, 10));
        }
        // E_gamma_addback vs (T_beta - T_implant)
        {
          std::string histName = "hEgaTib0n" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma addback vs Tbeta-Timp (0n)", 2000, 0, 2000, 2000, -10, 10));
        }
        {
          std::string histName = "hEgaTibn" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma addback vs Tbeta-Timp (1n)", 2000, 0, 2000, 2000, -10, 10));
        }
        {
          std::string histName = "hEgaTibnn" + isotopeName;
          fHistArray->Add(new TH2F(histName.c_str(), "Egamma addback vs Tbeta-Timp (2n)", 2000, 0, 2000, 2000, -10, 10));
        }
      }
      if (hist_group_map_.count("GAMMA_GAMMA"))
      {
        // gamma-gamma
        {
          std::string histName = "hEEg0n" + isotopeName;
          TH2F *h = new TH2F(histName.c_str(), "Gamma-gamma (0n)", 4000, 0, 4000, 4000, 0, 4000);
          h->SetXTitle("D4 energy");
          h->SetYTitle("G7 energy");
          fHistArray->Add(h);
        }
        {
          std::string histName = "hEEgn" + isotopeName;
          TH2F *h = new TH2F(histName.c_str(), "Gamma-gamma (1n)", 4000, 0, 4000, 4000, 0, 4000);
          h->SetXTitle("D4 energy");
          h->SetYTitle("G7 energy");
          fHistArray->Add(h);
        }
        {
          std::string histName = "hEEgnn" + isotopeName;
          TH2F *h = new TH2F(histName.c_str(), "Gamma-gamma (2n)", 4000, 0, 4000, 4000, 0, 4000);
          h->SetXTitle("D4 energy");
          h->SetYTitle("G7 energy");
          fHistArray->Add(h);
        }
      }
    }
  }
  catch (const std::bad_alloc&) {
    std::cout << "bad_alloc at histogram creation in " << isotopeName << std::endl;
    exit(300);
  }
  TIter next(fHistArray);
  while (TObject *obj = next())
  {
    outputList->Add(obj);
  }

  return;
}

bool hIsotope::IsInside(Double_t x, Double_t y) const
{
  if (myCutG)
  {
    return myCutG->IsInside(x, y);
  }
  else
  {
    return (ellipseA * TMath::Power(x - ellipseX0, 2) + ellipseB * TMath::Power(y - ellipseY0, 2) < 1);
  }
}
