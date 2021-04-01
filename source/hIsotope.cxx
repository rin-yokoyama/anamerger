// // // // // // // // // // // // // // // // // // 
// 
// Definition of an structure which gathers
// all the histograms associated with ONE isotope
// 
// VERSION 2.0. NEW HISTOGRAMS ADDED!!
// 
// // // // // // // // // // // // // // // // // // 

#include "hIsotope.h"

hIsotope::hIsotope(TCutG * aCut, TList* outputList):
numberOfImplants(0),
numberOfPID(0),
isotopeName("")
{
  myCutG = nullptr;
  isotopeName = "";
  ellipseA = 0;
  ellipseB = 0;
  ellipseX0 = 0;
  ellipseY0 = 0;

  initializeCut( aCut );
  initializeHistos( outputList );
}

hIsotope::hIsotope( std::string isoname, Double_t a, Double_t b, Double_t x0, Double_t y0, TList* outputList):
numberOfImplants(0),
numberOfPID(0),
isotopeName("")
{
  myCutG = nullptr;
  isotopeName = isoname;
  ellipseA = a;
  ellipseB = b;
  ellipseX0 = x0;
  ellipseY0 = y0;

  initializeHistos( outputList );
}


hIsotope::~hIsotope(){
//  if(fHistArray)
//    delete fHistArray;
}

void hIsotope::initializeCut(TCutG * aCut)
{
  std::string tempName4Cut( aCut->GetName() );   
  tempName4Cut="i"+tempName4Cut;
  myCutG=(TCutG *)(aCut->Clone( tempName4Cut.c_str() ));
  return;
}

void hIsotope::initializeHistos(TList* outputList)
{
    if( myCutG ){
      // get the name of the isotope
      // it must be at the end
      std::string tempName4Cut( myCutG->GetName() );
      size_t lastindex = tempName4Cut.find_last_of("CUTG");
      if( lastindex == size_t(-1) ) lastindex = tempName4Cut.find_last_of("cutg");
      // the name of TCutG could be "blablaCUTG82ga"
      // and the name of the isotope will be "82ga"
      std::cout << lastindex << " " << (tempName4Cut.substr(3+1, tempName4Cut.size())).c_str() << std::endl;
      isotopeName = tempName4Cut.substr(4+1, tempName4Cut.size());
    }

    fHistArray = new TObjArray();
  // initialize histograms
#ifdef ISOMER

    {
      std::string histName = "hETgI"+isotopeName;
      TH2F* h = new TH2F(histName.c_str(),"ET plot for isomers",4000,0,4000,2000,-100E-6,100E-6);
      h->SetXTitle("Gamma-ray energy (keV)");
      h->SetYTitle("Tgamma - Timplant (s)");
      fHistArray->Add(h); //ET histograms
    }
    {
      std::string histName = "hEEgI"+isotopeName;
      TH2F* h = new TH2F(histName.c_str(),"gamma-gamma plot for isomer (1 - 5 us)",4000,0,4000,4000,0,4000);
      h->SetXTitle("D4 energy");
      h->SetYTitle("G7 energy");
      fHistArray->Add(h); //Tib histograms
    }
#endif

#ifdef BETA
#ifdef DECAYCURVE
    // Tbeta - Timp decay curves
    {
      std::string histName = "hTib"+isotopeName;
      fHistArray->Add(new TH1F(histName.c_str(),"",20000,-10,10)); //Tib histograms
    }
    {
      std::string histName = "hTib0n"+isotopeName;
      fHistArray->Add(new TH1F(histName.c_str(),"",20000,-10,10)); //Tib histograms
    }
    {
      std::string histName = "hTibn"+isotopeName;
      fHistArray->Add(new TH1F(histName.c_str(),"",20000,-10,10)); //Tib histograms
    }
    {
      std::string histName = "hTibnn"+isotopeName;
      fHistArray->Add(new TH1F(histName.c_str(),"",20000,-10,10)); //Tib histograms
    }
    {
      std::string histName = "hTib3n"+isotopeName;
      fHistArray->Add(new TH1F(histName.c_str(),"",20000,-10,10)); //Tib histograms
    }
#endif
#ifdef GAMMA_ET
    // Tbeta - Timp vs gamma E
    {
      std::string histName = "hEgTib0n"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10,10)); //ET histograms
    }
    {
      std::string histName = "hEgTibn"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10,10)); //ET histograms
    }
    {
      std::string histName = "hEgTibnn"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10,10)); //ET histograms
    }
    // Tgamma - Tbeta vs gamma E
    {
      std::string histName = "hEgTbg0n"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10000,10000)); //ET histograms
    }
    {
      std::string histName = "hEgTbgn"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10000,10000)); //ET histograms
    }
    {
      std::string histName = "hEgTbgnn"+isotopeName;
      fHistArray->Add(new TH2F(histName.c_str(),"",10000,0,10000,2000,-10000,10000)); //ET histograms
    }

#endif
#ifdef GAMMA_GAMMA
    // gamma-gamma
    {
      std::string histName = "hEEg0n"+isotopeName;
      TH2F* h = new TH2F(histName.c_str(),"",4000,0,8000,4000,0,8000);
      h->SetXTitle("D4 energy");
      h->SetYTitle("G7 energy");
      fHistArray->Add(h); //Tib histograms
    }
    {
      std::string histName = "hEEgn"+isotopeName;
      TH2F* h = new TH2F(histName.c_str(),"",4000,0,8000,4000,0,8000);
      h->SetXTitle("D4 energy");
      h->SetYTitle("G7 energy");
      fHistArray->Add(h); //Tib histograms
    }
    {
      std::string histName = "hEEgnn"+isotopeName;
      TH2F* h = new TH2F(histName.c_str(),"",4000,0,8000,4000,0,8000);
      h->SetXTitle("D4 energy");
      h->SetYTitle("G7 energy");
      fHistArray->Add(h); //Tib histograms
    }
#endif
#endif
    TIter next(fHistArray);
    while( TObject* obj = next() ){ 
      outputList->Add(obj);
    }
   
  return;
}

bool hIsotope::IsInside(Double_t x, Double_t y) const
{
  if(myCutG){
    return myCutG->IsInside(x,y);
  }else{
    return (ellipseA*TMath::Power(x-ellipseX0,2)+ellipseB*TMath::Power(y-ellipseY0,2)<1);
  }
}
