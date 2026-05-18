#include "TGlauberMC.h"

using namespace std;


//---------------------------------------------------------------------------------
// ClassImp(TGlauNucleon)
//---------------------------------------------------------------------------------
void TGlauNucleon::RotateXYZ(Double_t phi, Double_t theta)
{
  TVector3 v(fX,fY,fZ);
  TVector3 vr;
  vr.SetMagThetaPhi(1,theta,phi);
  v.RotateUz(vr);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

//---------------------------------------------------------------------------------
void TGlauNucleon::RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ)
{
  TVector3 v(fX,fY,fZ);
  v.RotateX(psiX);
  v.RotateY(psiY);
  v.RotateZ(psiZ);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn, Double_t omega, Double_t G) 
{ // NN collision profile from https://arxiv.org/abs/1307.0636
  if ((omega<0) || (omega>2))
    return 0;
  Double_t R2 = snn/10./TMath::Pi();
  TF1 *nnprof = new TF1("nnprofgamma","[2]*(1-TMath::Gamma([0],[1]*x^2))",0,5);
  nnprof->SetParameters(1./omega,G/omega/R2,G);
  return nnprof;
}
//---------------------------------------------------------------------------------
TF1 *getNNProfDist(Double_t snn, Double_t omega, Double_t G) 
{ // NN collision profile from https://arxiv.org/abs/1307.0636
  if ((omega<0) || (omega>2))
    return 0;
  Double_t R2 = snn/10./TMath::Pi();
  TF1 *nnprof = new TF1("nnprofgammadist","[2]*x*(1-TMath::Gamma([0],[1]*x^2))",0,5);
  nnprof->SetParameters(1./omega,G/omega/R2,G*TMath::TwoPi());
  return nnprof;
}
//---------------------------------------------------------------------------------
TF1 *getNNHijingDist(Double_t signn, Double_t mu) 
{
  Double_t b02  = 0.5 * signn * 0.1 / TMath::Pi(); //sigS=57mb
  b02 = TMath::Sqrt(b02);
  TF1* eik = new TF1("eikdist", "2*pi*x*(1-exp(-[2]*[0]^2/[1]*([0]*x)^3*TMath::BesselK(3,[0]*x)))", 0, 10);
  eik->SetParameter(0, mu/b02);
  eik->SetParameter(1, 96);
  Double_t k=1;
  while (1) {
    eik->SetParameter(2,k);
    Double_t val=eik->Integral(0,5);
    Double_t diff=val-signn*0.1;
    if (diff<-0.01)
      k*=1.01;
    else if (diff>0.01)
      k*=0.99;
    else break;
  }
  eik->SetLineWidth(4);
  eik->SetLineStyle(9);
  eik->SetLineColor(1);
  return eik;
}

//---------------------------------------------------------------------------------
TF1 *getNNHijing(Double_t signn, Double_t mu) 
{
  Double_t b02  = 0.5 * signn * 0.1 / TMath::Pi(); //sigS=57mb
  b02 = TMath::Sqrt(b02);
  TF1* eik = new TF1("eik", "1-exp(-[2]*[0]^2/[1]*([0]*x)^3*TMath::BesselK(3,[0]*x))", 0, 10);
  eik->SetParameter(0, mu/b02);
  eik->SetParameter(1, 96);
  TF1 *dummy=getNNHijingDist(signn,mu);
  eik->SetParameter(2,dummy->GetParameter(2));
  delete dummy;
  eik->SetLineWidth(4);
  eik->SetLineStyle(9);
  eik->SetLineColor(1);
  return eik;
}

//---------------------------------------------------------------------------------
TF1 *getNNPythiaDist(Double_t signn, Double_t m, Double_t rp) 
{
  TF1* pyt = new TF1("pytdist", "2*pi*x*(1-exp(-[2]*TMath::Exp(-TMath::Power(x/[1],[0]))))",0,10);
  pyt->SetParameters(m,rp,1); //m=1.85 Monash tune
  Double_t k=1;
  while (1) {
    pyt->SetParameter(2,k);
    Double_t val=pyt->Integral(0,5);
    Double_t diff=val-signn*0.1;
    if (diff<-0.01)
      k*=1.01;
    else if (diff>0.01)
      k*=0.99;
    else break;
  }
  pyt->SetLineWidth(4);
  pyt->SetLineStyle(5);
  pyt->SetLineColor(1);
  return pyt;
}

//---------------------------------------------------------------------------------
TF1 *getNNPythia(Double_t signn, Double_t m, Double_t rp) 
{
  TF1* pyt = new TF1("pyt", "1-exp(-[2]*TMath::Exp(-TMath::Power(x/[1],[0])))",0,10);
  pyt->SetParameters(m,rp,1); //m=1.85 Monash tune
  TF1 *dummy=getNNPythiaDist(signn,m);
  pyt->SetParameter(2,dummy->GetParameter(2));
  delete dummy;
  pyt->SetLineWidth(4);
  pyt->SetLineStyle(5);
  pyt->SetLineColor(1);
  return pyt;
}

//---------------------------------------------------------------------------------
TF1 *getNNTrentoDist(Double_t signn, Double_t w) 
{
  const Double_t maxb = 6;
  TF1* t = new TF1("trentodist","2*pi*x*(1-exp(-[0]*exp(-x^2/4/[1]^2)))",0,maxb);
  t->SetParameters(1,w); 
  Double_t k=1;
  while (1) {
    t->SetParameter(0,k);
    Double_t val=t->Integral(0,maxb);
    Double_t diff=val-signn*0.1;
    if (diff<-0.01)
      k*=1.01;
    else if (diff>0.01)
      k*=0.99;
    else break;
  }
  t->SetLineWidth(4);
  t->SetLineStyle(9);
  t->SetLineColor(1);
  return t;
}

//---------------------------------------------------------------------------------
TF1 *getNNTrento(Double_t signn, Double_t w) 
{
  TF1* t = new TF1("trentodist","1-exp(-[0]*exp(-x^2/4/[1]^2))",0,6);
  t->SetParameters(1,w); 
  TF1 *dummy=getNNTrentoDist(signn,w);
  t->SetParameter(0,dummy->GetParameter(0));
  delete dummy;
  t->SetLineWidth(4);
  t->SetLineStyle(9);
  t->SetLineColor(1);
  return t;
}

//---------------------------------------------------------------------------------
TF1 *getSigmaNNvsEnergy()
{
  TF1 *lnNs1 = new TF1("lnNs","[0]+[1]*pow(log(x*x),[2])",5,100000, TF1::EAddToList::kDefault);
  lnNs1->SetFillColor(19);
  lnNs1->SetFillStyle(0);
  lnNs1->SetLineColor(2);
  lnNs1->SetLineWidth(3);
  lnNs1->SetChisquare(17.60074);
  lnNs1->SetNDF(24);
  lnNs1->GetXaxis()->SetLabelFont(42);
  lnNs1->GetXaxis()->SetTitleOffset(1);
  lnNs1->GetXaxis()->SetTitleFont(42);
  lnNs1->GetYaxis()->SetLabelFont(42);
  lnNs1->GetYaxis()->SetTitleFont(42);
  lnNs1->SetParameter(0,28.84374);
  lnNs1->SetParError(0,0.5163632);
  lnNs1->SetParLimits(0,0,0);
  lnNs1->SetParameter(1,0.04584121);
  lnNs1->SetParError(1,0.0167335);
  lnNs1->SetParLimits(1,0,0);
  lnNs1->SetParameter(2,2.374257);
  lnNs1->SetParError(2,0.1231245);
  lnNs1->SetParLimits(2,1,4);
  return lnNs1;
}

//---------------------------------------------------------------------------------
Double_t getSigmaNN(Double_t energy)
{
  if (energy < 10) {
    cerr << "getSigmaNN: energy out of range (10-100000 GeV): " << energy << endl;
    return 0.;
  }
  if (energy < 100) {
    cerr << "getSigmaNN: warning, provided result may not be accuruate at energies lower than 100 GeV: " << energy << endl;
  }
  TF1 *lnNs1 = getSigmaNNvsEnergy();
  Double_t sigma = lnNs1->Eval(energy);
  delete lnNs1;
  return sigma;
}

//---------------------------------------------------------------------------------
TGraph *getSigmaNPvsEnergy_Bystricky()
{
  //Tkin (GeV)  sigma_inel(np) (mb) Error (mb)
  // from Bystricky_JPhysique48  Table VI
  const Double_t x[] = {
    0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550,
    0.575, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10,
    1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.60, 1.70, 1.80, 1.90,
    2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.20,
    3.40, 3.60, 3.80, 4.00, 4.20
  };
  const Double_t y[] = {
    0.003, 0.046, 0.175, 0.419, 0.783, 1.26, 1.84, 2.50, 3.23, 4.00, 4.81,
    5.64, 6.48, 8.15, 9.77, 11.31, 12.76, 14.10, 15.33, 16.46, 17.48, 18.42,
    19.26, 20.03, 20.73, 21.36, 21.93, 22.44, 22.91, 23.34, 23.72, 24.39,
    24.95, 25.42, 25.81, 26.14, 26.43, 26.68, 26.90, 27.10, 27.28, 27.45,
    27.60, 27.75, 27.89, 28.03, 28.30, 28.57, 28.84, 29.13, 29.42, 29.73
  };
  const Double_t yerr[] = {
    0.001, 0.002, 0.005, 0.007, 0.009, 0.010, 0.013, 0.015, 0.018, 0.022, 0.025,
    0.029, 0.033, 0.041, 0.049, 0.058, 0.069, 0.080, 0.092, 0.11, 0.12, 0.13,
    0.15, 0.16, 0.17, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.26, 0.28, 0.30,
    0.31, 0.33, 0.34, 0.36, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.52,
    0.57, 0.63, 0.69, 0.76, 0.83, 0.90
  };
  TGraphErrors *gr = new TGraphErrors(sizeof(x)/sizeof(x[0]),x,y,0,yerr);
  gr->SetName("sigma_NP_Bystricky");
  return gr;
}

//---------------------------------------------------------------------------------
Double_t getSigmaNP_Bystricky(Double_t energy)
{
  if ((energy < 0.3) || (energy > 4.2)) {
    cerr << "getSigmaNP_Bystricky: energy out of range (0.3-4.2 GeV): " << energy << endl;
    return 0.;
  }
  TGraph *gr = getSigmaNPvsEnergy_Bystricky();
  Double_t sigma = gr->Eval(energy);
  delete gr;
  return sigma;
}

//---------------------------------------------------------------------------------
TGraph *getSigmaPPvsEnergy_Bystricky()
{
  //Tkin (GeV)  sigma_inel(pp) (mb) Error (mb)
  //from Bystricky_JPhysique48  Table VI
  const Double_t x[] = {
    0.28, 0.29, 0.30, 0.31, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.92, 0.96, 1.00, 1.04, 1.08, 1.12, 1.16, 1.20, 1.24, 1.28, 1.32, 1.36, 1.40, 1.44, 1.48, 1.52, 1.56, 1.64, 1.68, 1.72, 1.76, 1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0
  };
  const Double_t y[] = {
    0.027, 0.051, 0.082, 0.120, 0.171, 0.321, 0.563, 0.921, 1.41, 2.05, 2.83, 3.74, 4.76, 5.87, 7.04, 9.48, 11.89, 14.14, 16.15, 17.89, 19.34, 20.53, 21.51, 22.33, 23.03, 23.81, 24.39, 24.89, 25.30, 25.66, 25.95, 26.21, 26.42, 26.59, 26.74, 26.86, 26.96, 27.04, 27.11, 27.17, 27.21, 27.28, 27.30, 27.32, 27.33, 27.34, 27.35, 27.36, 27.36, 27.36, 27.36, 27.37, 27.37, 27.38, 27.38, 27.40, 27.41, 27.43, 27.46, 27.49, 27.52, 27.97, 28.48, 28.92, 29.28, 29.55, 29.75, 29.89, 30.04, 30.10, 30.09, 30.05, 30.00, 29.88, 29.80, 29.78, 29.91, 30.13, 30.38, 30.66, 30.93, 31.19, 31.78, 32.23, 32.56, 32.79, 32.92, 32.97, 32.96, 32.89, 32.78, 32.63, 32.44, 32.24, 32.01
  };
  TGraph *gr = new TGraph(sizeof(x)/sizeof(x[0]),x,y);
  gr->SetName("sigma_PP _Bystricky");
  return gr;
}

//---------------------------------------------------------------------------------
Double_t getSigmaPP_Bystricky(Double_t energy)
{
  if ((energy < 0.28) || (energy > 425.)) {
    cerr << "getSigmaPP_Bystricky: energy out of range (0.28-425 GeV): " << energy << endl;
    return 0.;
  }
  TGraph *gr = getSigmaPPvsEnergy_Bystricky();
  Double_t sigma = gr->Eval(energy);
  delete gr;
  return sigma;
}

//---------------------------------------------------------------------------------
TGraph *getSigmaHardvsEnergy()
{
  const Double_t x[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};
  const Double_t y[] = {4.41642385E-03, 0.188464478, 0.649059832, 1.25692534, 1.92962086, 2.63335800, 3.34105325, 4.05063772, 4.73708200, 5.42901564, 11.6562014, 16.9493980, 21.6152058, 25.8751488, 29.7382126, 33.4025536, 36.8270454, 40.0877686, 43.2026596, 69.1051025, 89.9556122, 107.744072, 124.087372, 138.940948, 153.014984, 166.308075, 178.939285, 190.986938, 294.216583, 379.379303, 455.911163, 525.982910, 592.784119, 656.306885, 717.265808, 774.954590, 830.319702};
  TGraph *gr = new TGraph(sizeof(x)/sizeof(x[0]),x,y);
  gr->SetFillColor(19);
  gr->SetFillStyle(0);
  gr->SetLineColor(2);
  gr->SetLineWidth(3);
  gr->SetName("sigma_hard_HIJING");
  return gr;
} 

//---------------------------------------------------------------------------------
Double_t getSigmaHard(Double_t energy)
{
  if (energy < 10)
    return 0;
  if (energy>100000) {
    cerr << "getSigmaHard: warning provided value may not be accurate as energy out of range (10-100000 GeV): " << energy << endl;
  }
  TGraph *gr = getSigmaHardvsEnergy();
  Double_t sigma = gr->Eval(energy);
  delete gr;
  return sigma;
}

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA,
                      const char *sysB,
                      const Double_t signn,
                      const Double_t sigwidth,
                      const Double_t mind,
		                  const Double_t omega,
                      const Double_t noded,
                      const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->SetNodeDistance(noded);
  mcg->SetCalcLength(0);
  mcg->SetCalcArea(0);
  mcg->SetCalcCore(0);
  mcg->SetDetail(99);
  TString om;
  if ((omega>=0) && (omega<=11)) {
    mcg->SetOmega(omega);
    om=Form("-om%.1f",omega);
  }
  TString name;
  if (fname) 
    name = fname; 
  else {
    TString nd;
    if (noded>0) 
      nd=Form("-nd%.1f",noded);
    name = Form("%s%s%s.root",mcg->Str(),om.Data(),nd.Data());
  }
  mcg->Run(n);
  TFile out(name,"recreate",name,9);
  TNtuple  *nt=mcg->GetNtuple();
  nt->Write();
  out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA,           
                        const char *sysB,           
                        const Double_t signn,
                        const Double_t sigwidth,
                        const Double_t mind,
                        const Bool_t verbose,
			                  const Double_t bmin,
			                  const Double_t bmax,
			                  const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  TFile *out=0;
  if (fname) 
    out=new TFile(fname,"recreate",fname,9);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    //get an event with at least one collision
    mcg->Run(1);
    if (ievent%100==0)
      cout << "\r" << 100.*ievent/n << "% done" << flush;

    //access, save and (if wanted) print out nucleons
    TObjArray* nucleons=mcg->GetNucleons();
    if (!nucleons) 
      continue;
    if (out)
      nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

    if (verbose) {
      cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
      cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
      printf("Nucleus\t X\t Y\t Z\tNcoll\n");
      Int_t nNucls=nucleons->GetEntries();
      for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
        TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
        Char_t nucleus='A';
        if (nucl->IsInNucleusB()) 
	        nucleus='B';
        Double_t x=nucl->GetX();
        Double_t y=nucl->GetY();
        Double_t z=nucl->GetZ();
        Int_t ncoll=nucl->GetNColl();
        printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
      }
    }
  }
  cout << endl << "Done!" << endl;
  if (out) {
    TNtuple *nt = mcg->GetNtuple();
    nt->Write();
    if (verbose)
      out->ls();
    delete out;
  }
}

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
		                   const Double_t bmin,
		                   const Double_t bmax,
                       const char *fname)
{
  // Run Glauber and store ntuple with smeared eccentricities in file.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out)
    return;
  TNtuple *nt = new TNtuple("nt","nt",
      "Npart:Ncoll:B:Psi1P:Ecc1P:Psi2P:Ecc2P:Psi3P:Ecc3P:Psi4P:Ecc4P:Psi5P:Ecc5P:Psi1G:Ecc1G:Psi2G:Ecc2G:Psi3G:Ecc3G:Psi4G:Ecc4G:Psi5G:Ecc5G:Sx2P:Sy2P:Sx2G:Sy2G");
  nt->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }
      for (Int_t i = 0; i<BN; ++i) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
        if (!nucleonB->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }

      Double_t MeanX  = 0;
      Double_t MeanY  = 0;
      Double_t MeanX2 = 0;
      Double_t MeanY2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
        MeanX  += xvals[i];
        MeanY  += yvals[i];
        MeanX2 += xvals[i]*xvals[i];
        MeanY2 += yvals[i]*yvals[i];
      }
      MeanX  /= ni;
      MeanY  /= ni;
      MeanX2 /= ni;
      MeanY2 /= ni;
      sx2g        += MeanX2-MeanX*MeanX;
      sy2g        += MeanY2-MeanY*MeanY;

      for (Int_t j = 1; j<9; ++j) {
        for (Int_t i = 0; i<ni; ++i) {
          Double_t x   = xvals[i] - MeanX;
          Double_t y   = yvals[i] - MeanY;
          Double_t r   = TMath::Sqrt(x*x+y*y);
          Double_t phi = TMath::ATan2(y,x);
          Double_t w = j;
          if (j==1)
            w = 3; // use r^3 weighting for Ecc1/Psi1
          cosphi[j] += TMath::Power(r,w)*TMath::Cos(j*phi);
          sinphi[j] += TMath::Power(r,w)*TMath::Sin(j*phi);
          rn[j]     += TMath::Power(r,w);
        }
      }
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    Float_t v[27]; Int_t i=0;
    v[i++] = mcg->GetNpart();
    v[i++] = mcg->GetNcoll();
    v[i++] = mcg->GetB();
    v[i++] = mcg->GetPsi(1); // point-like calculation values
    v[i++] = mcg->GetEcc(1);
    v[i++] = mcg->GetPsi(2);
    v[i++] = mcg->GetEcc(2);
    v[i++] = mcg->GetPsi(3);
    v[i++] = mcg->GetEcc(3);
    v[i++] = mcg->GetPsi(4);
    v[i++] = mcg->GetEcc(4);
    v[i++] = mcg->GetPsi(5);
    v[i++] = mcg->GetEcc(5);
    v[i++] = psi[1];         // Gaussian smeared values
    v[i++] = ecc[1];
    v[i++] = psi[2];
    v[i++] = ecc[2];
    v[i++] = psi[3];
    v[i++] = ecc[3];
    v[i++] = psi[4];
    v[i++] = ecc[4];
    v[i++] = psi[5];
    v[i++] = ecc[5];
    v[i++] = mcg->GetSx2();
    v[i++] = mcg->GetSy2();
    v[i++] = sx2g/NSAMP;
    v[i++] = sy2g/NSAMP;
    nt->Fill(v);
  }

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
void runAndOutputLemonTree(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
		                   const Double_t bmin,
		                   const Double_t bmax,
		                   const Bool_t   ogrid,
                       const char *fname)
{
  // Run Glauber and store Lemon TTree in format needed for IP-Jazma input 

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);
  
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out) return;

  // create new TTree with MC Glauber information for input to IP-Jazma
  const Int_t lemonmaxNucleons = 500;
  Int_t       lemonnpart;
  Int_t       lemonncoll;
  Int_t       lemonnparta;
  Int_t       lemonnpartb;  
  Float_t     lemonb;                        // collision impact parameter
  Float_t     lemoneccgaus[10];
  Float_t     lemoneccpoint[10];  
  Int_t       lemonnproj;
  Int_t       lemonntarg;
  Float_t     lemonxproj[lemonmaxNucleons];  // x,y,z coordinates for all nucleons
  Float_t     lemonyproj[lemonmaxNucleons];  // note these must be in the global coordinate frame
  Float_t     lemonzproj[lemonmaxNucleons];
  Float_t     lemonxtarg[lemonmaxNucleons];  // x,y,z coordinates for all nucleons
  Float_t     lemonytarg[lemonmaxNucleons];  // note these must be in the global coordinate frame
  Float_t     lemonztarg[lemonmaxNucleons];
  
  TTree *lemon = new TTree("lemon","lemon");
  lemon->Branch("npart",&lemonnpart,"npart/I");
  lemon->Branch("nparta",&lemonnparta,"nparta/I");
  lemon->Branch("npartb",&lemonnpartb,"npartb/I");  
  lemon->Branch("ncoll",&lemonncoll,"ncoll/I");
  lemon->Branch("b",&lemonb,"b/F");
  lemon->Branch("eccgaus",lemoneccgaus,"eccgaus[10]/F");
  lemon->Branch("eccpoint",lemoneccpoint,"eccpoint[10]/F");  
  lemon->Branch("nproj",&lemonnproj,"nproj/I");
  lemon->Branch("ntarg",&lemonntarg,"ntarg/I");  
  lemon->Branch("xproj",lemonxproj,"xproj[500]/F");
  lemon->Branch("yproj",lemonyproj,"yproj[500]/F");
  lemon->Branch("zproj",lemonyproj,"zproj[500]/F");  
  lemon->Branch("xtarg",lemonxtarg,"xtarg[500]/F");
  lemon->Branch("ytarg",lemonytarg,"ytarg[500]/F");
  lemon->Branch("ztarg",lemonytarg,"ztarg[500]/F");    
  lemon->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);
  TF2* smearing_function = new TF2("smear_tf2", "TMath::Exp(-(x*x+y*y)/(2.*[0]*[0]))/(2*TMath::Pi()*[0]*[0])", 0, 10*sigs, 0, 10*sigs);
  smearing_function->SetParameter(0,sigs);
  
  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }
      for (Int_t i = 0; i<BN; ++i) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
        if (!nucleonB->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }

      Double_t MeanX  = 0;
      Double_t MeanY  = 0;
      Double_t MeanX2 = 0;
      Double_t MeanY2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
        MeanX  += xvals[i];
        MeanY  += yvals[i];
        MeanX2 += xvals[i]*xvals[i];
        MeanY2 += yvals[i]*yvals[i];
      }
      MeanX  /= ni;
      MeanY  /= ni;
      MeanX2 /= ni;
      MeanY2 /= ni;
      sx2g        += MeanX2-MeanX*MeanX;
      sy2g        += MeanY2-MeanY*MeanY;

      for (Int_t j = 1; j<9; ++j) {
        for (Int_t i = 0; i<ni; ++i) {
          Double_t x   = xvals[i] - MeanX;
          Double_t y   = yvals[i] - MeanY;
          Double_t r   = TMath::Sqrt(x*x+y*y);
          Double_t phi = TMath::ATan2(y,x);
          Double_t w = j;
          if (j==1)
            w = 3; // use r^3 weighting for Ecc1/Psi1
          cosphi[j] += TMath::Power(r,w)*TMath::Cos(j*phi);
          sinphi[j] += TMath::Power(r,w)*TMath::Sin(j*phi);
          rn[j]     += TMath::Power(r,w);
        }
      }
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    // fill lemon TTree variables for this event
    lemonnpart   = mcg->GetNpart();
    lemonnparta  = mcg->GetNpartA();
    lemonnpartb  = mcg->GetNpartB();    
    lemonncoll   = mcg->GetNcoll();
    lemonb       = mcg->GetB(); 

    lemoneccpoint[0] = 0.0;
    lemoneccpoint[1] = mcg->GetEcc(1);
    lemoneccpoint[2] = mcg->GetEcc(2);
    lemoneccpoint[3] = mcg->GetEcc(3);
    lemoneccpoint[4] = mcg->GetEcc(4);
    lemoneccpoint[5] = mcg->GetEcc(5);
    lemoneccpoint[6] = mcg->GetEcc(6);    
    lemoneccpoint[7] = mcg->GetEcc(7);    
    lemoneccpoint[8] = mcg->GetEcc(8);    
    lemoneccpoint[9] = mcg->GetEcc(9);    

    lemoneccgaus[0] = 0.0;
    lemoneccgaus[1] = ecc[1];    
    lemoneccgaus[2] = ecc[2];
    lemoneccgaus[3] = ecc[3];
    lemoneccgaus[4] = ecc[4];
    lemoneccgaus[5] = ecc[5];
    lemoneccgaus[6] = ecc[6];
    lemoneccgaus[7] = ecc[7];
    lemoneccgaus[8] = ecc[8];
    lemoneccgaus[9] = ecc[9];

    for (Int_t i = 0; i<AN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
      lemonxproj[i] = nucleonA->GetX();
      lemonyproj[i] = nucleonA->GetY();
      lemonzproj[i] = nucleonA->GetZ();      
    }
    for (Int_t i = 0; i<BN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
      lemonxtarg[i] = nucleonB->GetX();
      lemonytarg[i] = nucleonB->GetY();
      lemonztarg[i] = nucleonB->GetZ();      
    }
    
    lemon->Fill();
    
    //======================================================================
    // also include option to write out nucleon smeared energy density map
    //======================================================================

    if (ogrid) {

      const Int_t nbins = 1000; 
      const Int_t nbinsx = nbins;
      const Int_t nbinsy = nbins;
      const Double_t max_x = 7.5;
      
      // now create an energy density distribution (a.u.)
      TH2D* inited_hist = new TH2D(Form("inited_event%i",ievent), ";x;y;E [a.u.]", nbinsx, -max_x, max_x, nbinsy, -max_x, max_x);
      for (Int_t ybin=1; ybin<=nbinsy; ybin++) {
        for (Int_t xbin=1; xbin<=nbinsx; xbin++) {
          // get the center of the bin
          const Double_t xval = inited_hist->GetXaxis()->GetBinCenter(xbin);
          const Double_t yval = inited_hist->GetYaxis()->GetBinCenter(ybin);
          long double content = 0.;  // sum the contributions from all wounded nucleons
          
          for (Int_t i = 0; i<AN; ++i) {
            TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
            if (!nucleonA->IsWounded()) continue;   // skip non-wounded nucleons
            content += smearing_function->Eval(nucleonA->GetX() - xval, nucleonA->GetY() - yval);
          }
          for (Int_t i = 0; i<BN; ++i) {
            TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
            if (!nucleonB->IsWounded()) continue;   // skip non-wounded nucleons
            content += smearing_function->Eval(nucleonB->GetX() - xval, nucleonB->GetY() - yval);
          }
          inited_hist->SetBinContent(xbin, ybin, content);          
        }
      }
      inited_hist->Write();
      if (inited_hist) delete inited_hist;
    }
  } // end loop over events

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
void runAndCalcDens(
        const Int_t n,
		    const Double_t alpha,
		    const char *sysA,
		    const char *sysB,
		    const Double_t signn,
		    const Double_t mind,
		    const char *fname)
{
  // Run Glauber and store per event a density profile in x and y, calculated from participant and binary positions
  // with relative weight given by alpha.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);

  TFile *out = 0;
  TCanvas *c1 = 0;
  if (fname) {
    out = TFile::Open(fname,"recreate",fname,9);
    if (!out)
      return;
  } else {
    c1 = new TCanvas;
  }

  const Int_t NSAMP = 100;
  const Double_t wp = (1-alpha)/2;
  const Double_t wb = alpha;
  const Double_t sigs = TMath::Sqrt(signn/20/TMath::Pi()); //from arXiv::0711.3724

  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  TH2F *h2f = new TH2F("h2f",";x (fm);y (fm)",121,-15.5,15.5,121,-15.5,15.5);
  h2f->SetStats(0);
  
  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}
    h2f->Reset();
    h2f->SetName(Form("Event_%d",ievent));
    h2f->SetTitle(Form("Npart=%d, Ncoll=%d",mcg->GetNpart(), mcg->GetNcoll()));

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    for (Int_t i = 0; i<AN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
      if (!nucleonA->IsWounded())
      	continue;
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	      Double_t sr = rad->GetRandom();
	      Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	      h2f->Fill(xA+sr*TMath::Cos(sp),yA+sr*TMath::Sin(sp),wp);
      }
    }

    for (Int_t j = 0; j<BN; ++j) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
      if (!nucleonB->IsWounded())
	      continue;
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	      Double_t sr = rad->GetRandom();
	      Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	      h2f->Fill(xB+sr*TMath::Cos(sp),yB+sr*TMath::Sin(sp),wp);
      }
    }

    if (alpha>0) {
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t xA=nucleonA->GetX();
        Double_t yA=nucleonA->GetY();
        for (Int_t j = 0; j<BN; ++j) {
          TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
          if (!mcg->IsBC(i,j))
            continue;
          Double_t xB=nucleonB->GetX();
          Double_t yB=nucleonB->GetY();
          Double_t dX=(xA+xB)/2;
          Double_t dY=(yA+yB)/2;
          for (Int_t s=0; s<NSAMP; ++s) {
            Double_t sr = rad->GetRandom();
            Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
            h2f->Fill(dX+sr*TMath::Cos(sp),dY+sr*TMath::Sin(sp),wb);
          }
        }
      }
    }
    h2f->Scale(1./h2f->Integral());
    if (out) {
      h2f->Write();
    } else {
      h2f->Draw("colz");
      c1->Update();
      gSystem->Sleep(1000);
    }
  }
  if (out) {
    out->Write();
    out->Close();
    delete out;
  }
}

//---------------------------------------------------------------------------------
// ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------
TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
  TNamed(iname,""),
  fN(iN),fR(iR),fA(ia),fW(iw),fR2(0),fA2(0),fW2(0),fBeta2(0),fBeta4(0),
  fMinDist(0.4),fNodeDist(0.0),fSmearing(0.0),fRecenter(1),fLattice(0),fSmax(99),
  fF(0),fTrials(0),fNonSmeared(0),fFunc1(ifunc),fFunc2(0),fFunc3(0),fNucleons(0),
  fPhiRot(0),fThetaRot(0),fNucArr(0),fNucCounter(-1),fIsUsed(0),fMaxR(15),fWeight(0),
  fInputDens(0)
{
  if (fN==0) {
    cout << "Setting up nucleus " << iname << endl;
    Lookup(iname);
  }
}

TGlauNucleus::~TGlauNucleus()
{
  if (fIsUsed)
    delete fIsUsed;
  if (fNucleons)
    delete fNucleons;
  delete fFunc1;
  delete fFunc2;
  delete fFunc3;
  delete fInputDens;
  FreeNucArr();
}

void TGlauNucleus::AllocateNucleons()
{
  if (fN<=0) 
    return;
  if (fNucleons==0) {
    fNucleons=new TObjArray(fN);
    fNucleons->SetOwner();
  } else {
    if (fNucleons->GetEntries()>=fN)
      return;
    fNucleons->Clear();
    fNucleons->Expand(fN);
  }
  for (Int_t i=0; i<fN; ++i) {
    TGlauNucleon *nucleon=new TGlauNucleon(); 
    nucleon->SetType(0);
    if (i<fZ) 
      nucleon->SetType(1);
    fNucleons->Add(nucleon); 
  } 
}

void TGlauNucleus::RandomizeNucleons()
{
  for (Int_t i=0,iz=0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)fNucleons->At(i);
    Double_t frac=double(fZ-iz)/double(fN-i);
    Double_t rn=gRandom->Uniform(0,1);
    if (rn<frac) {
      nucleon->SetType(1);
      ++iz;
    } else {
      nucleon->SetType(0);
    }
  }
}

void TGlauNucleus::AllocateNucArr()
{
  FreeNucArr();
  if (fNucCounter<=0) return;
  fNucArr = new Float_t**[fNucCounter];
  for (Int_t i = 0; i < fNucCounter; i++) {
    fNucArr[i] = new Float_t*[fN];
    for (Int_t j = 0; j < fN; j++) {
      fNucArr[i][j] = new Float_t[4]; 
    }
  }
}

void TGlauNucleus::FreeNucArr()
{
  if (!fNucArr) return;
  for (Int_t i = 0; i < fNucCounter; i++) {
    for (Int_t j = 0; j < fN; j++) {
      delete[] fNucArr[i][j];
    }
    delete[] fNucArr[i];
  }
  delete[] fNucArr;
  fNucArr = 0;
  fNucCounter = -1;
}

void TGlauNucleus::ReadNucArr(const char *fname=0)
{
  if (fNucCounter>0) {
    cout << "Warning: Should only call this function once (fNucCounter=" << fNucCounter << ")" << endl;
    return;
  }
  TString filename;
  Int_t extraInfo = 0;
  if (fname) {
    filename = fname;
    extraInfo = 1;
  } else {
    TString tmpname = GetName();
    if (tmpname=="He3") {
      filename = "he3_plaintext.dat";
      fNucCounter = 13699;
      extraInfo = 2;
    } else if (tmpname=="H3") {
      filename = "h3_plaintext.dat";
      fNucCounter = 13729;
      extraInfo = 2;
    } else if (tmpname=="He4") {
      filename = "he4_plaintext.dat";
      fNucCounter = 50626;
    } else if (tmpname=="C") {
      filename = "carbon_plaintext.dat";
      fNucCounter = 9999;
      extraInfo = 3;
    } else if (tmpname=="O") {
      filename = "oxygen_plaintext.dat";
      fNucCounter = 10000;
    } else {
      filename = GetName();
    }
  }

  ifstream myfile;
  myfile.open(filename);
  if (!myfile.is_open()) {
    TString datfile = Form("dat/%s",filename.Data());
    myfile.open(datfile.Data());
  }
  if (!myfile) {
    cerr << "ERROR:  no file for nucleon configurations found with name = " << filename << endl;
    gSystem->Exit(123);
  }
  if (fN==0) { // get all information from the file
    string line;
    Int_t nprotons = 0;
    Int_t nnucleons = 0;
    Int_t numEvents = 0;
    char nucleusName[99];
    while (getline(myfile, line)) {
      if (line.find("# Number of events:") != string::npos) {
        sscanf(line.c_str(), "# Number of events: %d", &numEvents);
        break; // stop reading after the number of events
      } else if (line.find("# Number of protons:") != string::npos) {
        sscanf(line.c_str(), "# Number of protons: %d", &nprotons);
      } else if (line.find("# Number of nucleons:") != string::npos) {
        sscanf(line.c_str(), "# Number of nucleons: %d", &nnucleons);
      } else if (line.find("# Nucleus name:") != string::npos) {
        sscanf(line.c_str(), "# Nucleus name: %s", nucleusName);
        SetTitle(GetName());
        SetName(nucleusName);
      }
    }
    fN = nnucleons;
    fZ = nprotons;
    if (fN<=0 || fZ<=0 || fZ>fN) {
      cerr << "ERROR: number of nucleons (fN=" << fN << ") or protons (fZ=" << fZ << ") are not consistent" << endl;
      gSystem->Exit(123);
    }
    fNucCounter = numEvents;
    cout << "Found " << numEvents << " nucleon configurations in the file" << endl;
    cout << "Number of protons: " << nprotons << endl;
    cout << "Number of nucleons: " << nnucleons << endl;
  } else if (fNucCounter<=0) { // this is no longer really needed
    Int_t numEvents = 0;
    string line;
    while (getline(myfile, line)) {
      if (!line.empty()) numEvents++;
    }
    myfile.clear();
    myfile.seekg(0);
    cout <<"Found "<<numEvents<<" nucleon configurations in the file"<<endl;
    fNucCounter = numEvents;
  }
  // read in the data from the file
  AllocateNucArr();
  cout << "Reading file " << filename << " for " << fNucCounter << " nucleon configurations with fN = " << fN << endl;
  Int_t inputcounter = 0;
  while (myfile) {
    if (inputcounter >= fNucCounter) break;
    Float_t foo;
    for (Int_t i = 0; i<fN; ++i) {
      if (extraInfo==3) {
        myfile >> foo >> foo; // dummy info not used
      }
      myfile >> fNucArr[inputcounter][i][0] >> fNucArr[inputcounter][i][1] >> fNucArr[inputcounter][i][2];
      fNucArr[inputcounter][i][3] = -1; // no isospin info
      if (extraInfo==1) {
        myfile >> fNucArr[inputcounter][i][3]; // isospin info
      } else if (extraInfo==2) {
        myfile >> foo >> foo >> foo >> foo; // extra isospin info not used
      }
    }
    ++inputcounter;
  }
  myfile.close();
}

void TGlauNucleus::SelectFromNucArr()
{
  if (fNucleons==0) {
    cerr << "ERROR: nucleons not set, should not happen!" << endl;
    gSystem->Exit(123);
  }
  fTrials = 1;
  Bool_t ok = false;
  while (!ok) {
    Int_t nucidx=gRandom->Uniform(0,fNucCounter-1);
    ok = true;
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      Double_t x = fNucArr[nucidx][i][0];
      Double_t y = fNucArr[nucidx][i][1];
      Double_t z = fNucArr[nucidx][i][2];
      Int_t isospin = fNucArr[nucidx][i][3];
      if (isospin>=0)
        nucleon->SetType(isospin);
      nucleon->SetXYZ(x,y,z);
      nucleon->RotateXYZ(fPhiRot,fThetaRot);
      if (!TestMinDist(i,nucleon->GetX(),nucleon->GetY(),nucleon->GetZ())) {
        ++fTrials;
        ok = false;
        break;
      }
    }
  }
}

Double_t TGlauNucleus::CalcMinDist() const
{
  Double_t minr = 1e10;
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleons->At(i));
    Double_t xA = nucleonA->GetX();
    Double_t yA = nucleonA->GetY();
    Double_t zA = nucleonA->GetZ();
    for (Int_t j = i+1; j<fN; ++j) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleons->At(j));
      Double_t xB = nucleonB->GetX();
      Double_t yB = nucleonB->GetY();
      Double_t zB = nucleonB->GetZ();
      Double_t d = (xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB);
      if (d<minr)
        minr = d;
    }
  }
  return TMath::Sqrt(minr);
}

Double_t TGlauNucleus::CalcRmsRadius() const
{
  Double_t rrms = 0;
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    Double_t x = nucleon->GetX();
    Double_t y = nucleon->GetY();
    Double_t z = nucleon->GetZ();
    rrms += x*x+y*y+z*z;
  }
  rrms /= fN;
  return TMath::Sqrt(rrms);
}

void TGlauNucleus::Draw(Double_t xs, Int_t colp, Int_t cols)
{
  Double_t r = 0.5*TMath::Sqrt(xs/TMath::Pi()/10.);
  TEllipse en;
  en.SetLineStyle(1);
  en.SetLineWidth(1);
  en.SetFillStyle(1001);
  for (Int_t i = 0; i<fNucleons->GetEntries(); ++i) {
    TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
    if (!gn->IsSpectator()) {
      en.SetFillColor(colp);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    } else {
      en.SetFillColor(cols);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    }
  }
}

TF1* TGlauNucleus::GetDens(Bool_t n) const
{
  if (!fFunc1)
    return 0;
  TF1 *f = 0;
  if (fInputDens)
    f = new TF1(Form("%s_dens",fFunc1->GetName()),[&](double*x, double *p){return p[0]*fInputDens->Eval(x[0]);}, 0.0, fMaxR, 1);  
  else
    f = new TF1(Form("%s_dens",fFunc1->GetName()),Form("[0]*(%s)/x/x",fFunc1->GetName()),0,fMaxR);
  f->SetNpx(1000);
  f->SetParameter(0,1);
  for (Int_t i=0;i<fFunc1->GetNpar();++i)
    f->SetParameter(i+1,fFunc1->GetParameter(i));
  Double_t norm = fZ/fFunc1->Integral(0,fMaxR,0.001);
  if (n)
    norm /= 4*TMath::Pi();
  f->SetParameter(0,norm);
  return f;
} 

Double_t TGlauNucleus::GetSqrtMeanR2() const
{
  if (!fFunc1)
    return 0;
  TF1 *f = 0;
  if (fInputDens)
    f = new TF1(Form("%s_dummy",fFunc1->GetName()),[&](double*x, double *){return fInputDens->Eval(x[0])*pow(x[0],4);}, 0.0, fMaxR, 0);
  else
    f = new TF1(Form("%s_dummy",fFunc1->GetName()),Form("x^2*%s",fFunc1->GetName()),0,fMaxR);
  for (Int_t i=0;i<fFunc1->GetNpar();++i)
    f->SetParameter(i,fFunc1->GetParameter(i));
  Double_t meanr2 = f->Integral(0,fMaxR,0.001);
  delete f;
  Double_t norm = fFunc1->Integral(0,fMaxR,0.001);
  Double_t ret = TMath::Sqrt(meanr2/norm);
  if (TMath::IsNaN(ret))
    ret = 0;
  return ret;
}

void TGlauNucleus::Lookup(const char* name)
{
  SetName(name);
  Double_t r0=0, r1=0, r2=0; // if not specified the parameters are from DeVries, ATOMIC DATA AND NUCLEAR DATA TABLES 36,495536 (1987)
  if      (TString(name) == "p")       {fN = 1;   fR = 0.234;      fA = 0;      fW =  0;       fF = 0;  fZ=1;}
  else if (TString(name) == "pi")      {fN = 1;   fR = 0.234;      fA = 0;      fW =  0;       fF = 0;  fZ=1;}
  else if (TString(name) == "pg")      {fN = 1;   fR = 0.514;      fA = 0;      fW =  0;       fF = 9;  fZ=1;} 
  else if (TString(name) == "pdg")     {fN = 1;   fR = 1;          fA = 0;      fW =  0;       fF = 10; fZ=1;} // from arXiv:1101.5953
  else if (TString(name) == "dpf")     {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 1;  fZ=1;} // deuteron 2pf (tuned to Hulthen)
  else if (TString(name) == "dh")      {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 3;  fZ=1;} // deuteron Hulthen free
  else if (TString(name) == "d")       {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 4;  fZ=1;} // deuteron Hulthen constrained
  else if (TString(name) == "He3")     {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=1;} // read configurations from file
  else if (TString(name) == "H3")      {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "He4")     {fN = 4;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "C")       {fN = 12;  fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=6;} // read configurations from file
  else if (TString(name) == "Npar")    {fN = 14;  fR = 2.570;      fA = 0.0572; fW =  -0.0180; fF = 1;  fZ=7;} 
  else if (TString(name) == "O")       {fN = 16;  fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=8;} // read configurations from file
  else if (TString(name) == "Opar")    {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 1;  fZ=8;  fMaxR=7.5;} // WS parameterization
  else if (TString(name) == "Opar2")   {fN = 16;  fR = 1.850;      fA = 0.497;  fW =  0.912;   fF = 1;  fZ=8;  fMaxR=7.5;} // WS parameterization (from Trajectum OXYGEN_V4)
  else if (TString(name) == "Osat")    {fN = 16;  fR = 0;          fA = 0;      fW =  0;       fF = 90; fZ=8;  fMaxR=7.5; fInputDens=new TGraph("dat/16O_NNLO_sat.txt");}
  else if (TString(name) == "Odat")    {fN = 16;  fR = 2.608,      fA = 0.513;  fW = -0.051;   fF = 16; fZ=8;  fMaxR=7.5;} // WS parameterization (from Nuclear Physics A150 (1970) 631--654)
  else if (TString(name) == "Oho")     {fN = 16;  fR = 1.544;      fA = 1.833;  fW =  0;       fF = 15; fZ=8;  fMaxR=7.5;} // Harmonic oscillator parameterization
  else if (TString(name) == "Oho2")    {fN = 16;  fR = 1.506;      fA = 1.819;  fW =  0;       fF = 15; fZ=8;  fMaxR=7.5;} // Harmonic oscillator fit to data from Nuclear Physics A150 (1970) 631--654
  else if (TString(name) == "Ne")      {fN = 20;  fR = 2.805;      fA = 0.571;  fW =  0;       fF = 1;  fZ=10; fMaxR=8.5;}
  else if (TString(name) == "Ne2")     {fN = 20;  fR = 2.740;      fA = 0.572;  fW =  0;       fF = 1;  fZ=10; fMaxR=8.5;}
  else if (TString(name) == "Ne3")     {fN = 20;  fR = 2.791;      fA = 0.698;  fW = -0.168;   fF = 1;  fZ=10; fMaxR=8.5;}
  else if (TString(name) == "NeTr2")   {fN = 20;  fR = 2.8;        fA = 0.57;   fW =  0;       fF = 8;  fZ=10; fBeta2=0.721; fMaxR=10;} // WS parameterization (TR_NEON_V2 from Trajectum)
  else if (TString(name) == "NeTr3")   {fN = 20;  fR = 2.7243;     fA = 0.4982; fW =  0;       fF = 7;  fZ=10; fBeta2=0.4899; fBeta3=0.2160; fBeta4=0.3055; fGamma=0; fMaxR=10;} // WS parameterization (TR_NEON_V3 from Trajectum)
  else if (TString(name) == "Al")      {fN = 27;  fR = 3.34;       fA = 0.580;  fW =  0.0;     fF = 8;  fZ=13; fBeta2=-0.448; fBeta4=0.239; fMaxR=10;}
  else if (TString(name) == "Si")      {fN = 28;  fR = 3.34;       fA = 0.580;  fW = -0.233;   fF = 1;  fZ=14; fMaxR=10;}
  else if (TString(name) == "Si2")     {fN = 28;  fR = 3.34;       fA = 0.580;  fW =  0;       fF = 8;  fZ=14; fBeta2=-0.478; fBeta4=0.250; fMaxR=10;}
  else if (TString(name) == "S")       {fN = 32;  fR = 2.54;       fA = 2.191;  fW =  0.16;    fF = 2;  fZ=16; fMaxR=10;}
  else if (TString(name) == "Ar")      {fN = 40;  fR = 3.53;       fA = 0.542;  fW =  0;       fF = 1;  fZ=18; fMaxR=10;}
  else if (TString(name) == "Ca")      {fN = 40;  fR = 3.766;      fA = 0.586;  fW = -0.161;   fF = 1;  fZ=20; fMaxR=10;}
  else if (TString(name) == "Ni")      {fN = 58;  fR = 4.309;      fA = 0.517;  fW = -0.1308;  fF = 1;  fZ=28; fMaxR=10;}
  else if (TString(name) == "Cu")      {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 1;  fZ=29; fMaxR=10;}
  else if (TString(name) == "Curw ")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 12; fZ=29; r0=1.00898; r1=-0.000790403; r2=-0.000389897; fMaxR=10;} 
  else if (TString(name) == "Cu2")     {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 8;  fZ=29; fBeta2=0.162; fBeta4=-0.006; fMaxR=10;}  
  else if (TString(name) == "Cu2rw")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 14; fZ=29; fBeta2=0.162; fBeta4=-0.006; r0=1.01269; r1=-0.00298083; r2=-9.97222e-05; fMaxR=10;}  
  else if (TString(name) == "CuHN")    {fN = 63;  fR = 4.28;       fA = 0.5;    fW =  0;       fF = 1;  fZ=29; fMaxR=10;} // from arXiv:0904.4080v1
  // Nb93  <r^2>^(1/2)= 4.324    R= 4.9853(1)  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234  from Landolt-Börnstein
  else if (TString(name) == "Nb93LB")  {fN = 93; fR = 4.9853;      fA = 0.5234; fW =  0;       fF = 1;  fZ=41;}
  // Zr96  <r^2>^(1/2)= 4.349    R= 5.0212(2)  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234  from Landolt-Börnstein
  else if (TString(name) == "Zr96LB")  {fN = 96; fR = 5.0212;      fA = 0.5234; fW =  0;       fF = 1; fZ=40;}
  // Ru96  <r^2>^(1/2)= 4.393    R= 5.0845(1)  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234  from Landolt-Börnstein
  else if (TString(name) == "Ru96LB")  {fN = 96; fR = 5.0845;      fA = 0.5234; fW =  0;       fF = 1; fZ=44;}
  // Ag107  <r^2>^(1/2)= 4.544    R= 5.3006(1)  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234  from Landolt-Börnstein
  else if (TString(name) == "Ag107LB") {fN = 107; fR = 5.3006;     fA = 0.5234; fW =  0;       fF = 1;  fZ=47;}
  // Ag109  <r^2>^(1/2)= 4.565    R= 5.3306(1)  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234  from Landolt-Börnstein
  else if (TString(name) == "Ag109LB") {fN = 109; fR = 5.3306;     fA = 0.5234; fW =  0;       fF = 1; fZ=47;} 
   //from GiBUU Lenske et al.
  else if (TString(name) == "Ag107pn") {fN = 107; fR = 5.2731;     fA = 0.4749;  fW =  0;      fF = 11; fZ=47; fR2=5.4262;fA2=0.4776; fW2=0;}
  else if (TString(name) == "Ag109pn") {fN = 109; fR = 5.2943;     fA = 0.4729;  fW =  0;      fF = 11; fZ=47; fR2=5.4762;fA2=0.4788; fW2=0;}
  //from HFB14
  else if (TString(name) == "Ag107pnHFB14") {fN = 107; fR = 5.2875;fA = 0.4788;  fW =  0;      fF = 11; fZ=47; fR2=5.287; fA2=0.5498; fW2=0;}
  else if (TString(name) == "Ag109pnHFB14") {fN = 109; fR = 5.3160;fA = 0.4776;  fW =  0;      fF = 11; fZ=47; fR2=5.3246;fA2=0.5593; fW2=0;}
  // 50Sn  *stable isotopes
  else if (TString(name) == "Sn112")   {fN = 112; fR = 5.3714; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn114")   {fN = 114; fR = 5.3943; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn116")   {fN = 116; fR = 5.4173; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn117")   {fN = 117; fR = 5.1241; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn118")   {fN = 118; fR = 5.4391; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn119")   {fN = 119; fR = 5.4431; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn120")   {fN = 120; fR = 5.4588; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn122")   {fN = 122; fR = 5.4761; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn124")   {fN = 124; fR = 5.4907; fA = 0.5234; fW =  0;           fF = 1; fZ=50;}
  else if (TString(name) == "Sn112pr3"){fN = 112; fR = 4.962; fA = 2.638/(4.*TMath::Log(3.)); fW =  0.285; fF = 1; fZ=50;}
  else if (TString(name) == "Sn114pr3"){fN = 114; fR = 4.971; fA = 2.636/(4.*TMath::Log(3.)); fW =  0.320; fF = 1; fZ=50;}
  else if (TString(name) == "Sn116pr3"){fN = 116; fR = 5.062; fA = 2.625/(4.*TMath::Log(3.)); fW =  0.272; fF = 1; fZ=50;}
  else if (TString(name) == "Sn117pr3"){fN = 117; fR = 5.058; fA = 2.625/(4.*TMath::Log(3.)); fW =  0.295; fF = 1; fZ=50;}
  else if (TString(name) == "Sn118pr3"){fN = 118; fR = 5.072; fA = 2.623/(4.*TMath::Log(3.)); fW =  0.304; fF = 1; fZ=50;}
  else if (TString(name) == "Sn119pr3"){fN = 119; fR = 5.100; fA = 2.618/(4.*TMath::Log(3.)); fW =  0.290; fF = 1; fZ=50;}
  else if (TString(name) == "Sn120pr3"){fN = 120; fR = 5.110; fA = 2.619/(4.*TMath::Log(3.)); fW =  0.292; fF = 1; fZ=50;}
  else if (TString(name) == "Sn122pr3"){fN = 122; fR = 5.088; fA = 2.611/(4.*TMath::Log(3.)); fW =  0.378; fF = 1; fZ=50;}
  else if (TString(name) == "Sn124pr3"){fN = 124; fR = 5.150; fA = 2.615/(4.*TMath::Log(3.)); fW =  0.311; fF = 1; fZ=50;}
  // 50Sn  non-stable isotopes from Angeli2013  https://inspirehep.net/literature/1611365
  else if (TString(name) == "Sn108")   {fN =108; fR = 5.3274;      fA = 0.5234; fW =  0;       fF = 1;  fZ=50;} // Sn108  <r^2>^(1/2)= 4.5605  --> fR = 5.3274
  else if (TString(name) == "Sn132")   {fN =132; fR = 5.5387;      fA = 0.5234; fW =  0;       fF = 1;  fZ=50;} // Sn132  <r^2>^(1/2)= 4.7093  --> fR = 5.5387
  else if (TString(name) == "Xe")      {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 1;  fZ=54;} // adapted from arXiv:1703.04278
  else if (TString(name) == "Xes")     {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 1;  fZ=54;} // scale from Sb (Antimony, A=122, r=5.32) by 1.019 = (129/122)**0.333
  else if (TString(name) == "Xe2")     {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.161; fBeta4=-0.003;} // adapted from arXiv:1703.04278 and Z. Physik (1974) 270: 113
  else if (TString(name) == "Xe2a")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.18; fBeta4=0;} // ALICE parameters (see public note from 2018 at https://cds.cern.ch/collection/ALICE%20Public%20Notes?ln=en)
  else if (TString(name) == "Xerw")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 12; fZ=54; r0=1.00911; r1=-0.000722999; r2=-0.0002663;}
  else if (TString(name) == "Xesrw")   {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 12; fZ=54; r0=1.0096; r1=-0.000874123; r2=-0.000256708;}
  else if (TString(name) == "Xe2arw")  {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 14; fZ=54; fBeta2=0.18; fBeta4=0; r0=1.01246; r1=-0.0024851; r2=-5.72464e-05;} 
  // W184 <r^2>^(1/2)= 5.373    R= 6.3599  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234   from Landolt-Börnstein
  else if (TString(name) == "W184LB")  {fN = 184; fR = 6.3599;     fA = 0.523;  fW =  0;       fF = 1;  fZ=74;}
  else if (TString(name) == "W")       {fN = 186; fR = 6.58;       fA = 0.480;  fW =  0;       fF = 1;  fZ=74;}
  // W186 <r^2>^(1/2)= 5.381    R= 6.3839  fixed t = 2.3 --> a = t/(4 ln(3)) = 0.5234   from Landolt-Börnstein
  else if (TString(name) == "W186LB")  {fN = 186; fR = 6.3839;     fA = 0.523;  fW =  0;       fF = 1;  fZ=74;}
  else if (TString(name) == "Au")      {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 1;  fZ=79;}
  else if (TString(name) == "Aurw")    {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 12; fZ=79; r0=1.00899; r1=-0.000590908; r2=-0.000210598;}
  else if (TString(name) == "Au2")     {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 8;  fZ=79; fBeta2=-0.131; fBeta4=-0.031; }
  else if (TString(name) == "Au2rw")   {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 14; fZ=79; fBeta2=-0.131; fBeta4=-0.031; r0=1.01261; r1=-0.00225517; r2=-3.71513e-05;}
  else if (TString(name) == "AuHN")    {fN = 197; fR = 6.42;       fA = 0.44;   fW =  0;       fF = 1;  fZ=79;} // from arXiv:0904.4080v1
  else if (TString(name) == "Au197LB")     {fN = 197; fR = 6.5541; fA = 0.523;  fW =  0;       fF = 1;  fZ=79;}
  // from muonic and HBF calc from Landolt-Börnstein
  else if (TString(name) == "Au4pn")       {fN = 197; fR = 6.538;  fA = 0.465;  fW =  0;       fF = 11; fZ=79; fR2=6.794; fA2=0.483; fW2=0;}   //from GiBUU Lenske et al.
  else if (TString(name) == "Au197pnHFB14"){fN = 197; fR = 6.5831; fA = 0.4628; fW =  0;       fF = 11; fZ=79; fR2=6.6604; fA2=0.5464; fW2=0;} //HFB14
  else if (TString(name) == "Pb")      {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbrw")    {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 12; fZ=82; r0=1.00863; r1=-0.00044808; r2=-0.000205872;} //only Pb 207 was tested but should be the same for 208
  else if (TString(name) == "Pb*")     {fN = 208; fR = 6.624;      fA = 0.549;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "PbHN")    {fN = 208; fR = 6.65;       fA = 0.460;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbpn")    {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 11; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  else if (TString(name) == "Pbpnrw")  {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 13; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  else if (TString(name) == "U")       {fN = 238; fR = 6.188;      fA = 0.54;   fW =  0;       fF = 5;  fZ=92; fBeta2=1.77;} // Uranium from Heinz & Kuhlman, nucl-th/0411054, fR is defined as 6.8*0.91, fW=6.8*0.26
  else if (TString(name) == "U2")      {fN = 238; fR = 6.67;       fA = 0.44;   fW =  0;       fF = 8;  fZ=92; fBeta2=0.280; fBeta4=0.093;}
  else if (TString(name).BeginsWith("TR_")) {fN = 0; fR = 0; fA = 0; fW = 0; fF = 99; fZ=0;}  // Trajectum models (needs libtrnucgen)
  else if (TString(name).BeginsWith("input")) {fN = 0; fR = 0; fA = 0; fW = 0; fF = 6; fZ=0;} // Use input file with nucleon configurations
  else {
    cout << "Warning: Could not find nucleus " << name << " in lookup table" << endl;
    return;
  }
  switch (fF) {
    case 0: // Proton exp
      fFunc1 = new TF1(name,"x^2*exp(-x/[0])",0,fMaxR);
      fFunc1->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunc1 = new TF1(name,"x^2*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 2: // 3pG
      fFunc1 = new TF1(name,"x^2*(1+[2]*(x/[0])^2)/(1+exp((x^2-[0]^2)/[1]^2))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 3: // Hulthen (see nucl-ex/0603010)
    case 4: // same but constrain the neutron opposite to the proton event-by-event
      fFunc1 = new TF1(name,"x^2*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;
    case 5: // Ellipsoid (Uranium), box method
      break;
    case 6: // read from configuration file
      fFunc1 = 0;
      fMaxR = 0;
      ReadNucArr();
      break;
    case 7: // Deformed nuclei, box method
#ifndef HAVE_MATHMORE
      cerr << "Need libMathMore.so for deformed nuclei" << endl;
      gSystem->Exit(123);
#endif
      fFunc1 = 0; // no func: only need beta parameters and use uniform box distribution
      break;
    case 8: // Deformed nuclei, TF2 method
      fFunc3 = new TF2(name,"x^2*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4);
      break;
    case 9: // Proton gaus
      fFunc1 = new TF1(name,"x^2*exp(-x*x/[0]/[0]/2)",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 10: // Proton dgaus
      fFunc1 = new TF1(name,"x^2*((1-[0])/[1]^3*exp(-(x/[1])^2)+[0]/(0.4*[1])^3*exp(-x^2/(0.4*[1])^2))",0,5);
      fFunc1->SetParameter(0,0.5);
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // 3pF for proton and neutrons
      fFunc1 = new TF1(name,"x^2*(1+[2]*(x/[0])^2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      fFunc2 = new TF1(name,"x^2*(1+[2]*(x/[0])^2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2);
      break;
    case 12: // reweighted
      fFunc1 = new TF1(name,"x^2*(1+[2]*(x/[0])^2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,r0,r1,r2); 
      fRecenter=1;
      fSmax=0.1;
      break;
    case 13: // Pb for proton and neutrons reweighted
      fFunc1 = new TF1(Form("%s_prot",name),"x^2*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,1.00866,-0.000461484,-0.000203571);
      fFunc2 = new TF1(Form("%s_neut",name),"x^2*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2,1.00866,-0.000461484,-0.000203571);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 14: // Deformed nuclei, TF2 method, reweighted
      fFunc3 = new TF2(name,"x^2*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))/([4]+[5]*x+[6]*x^2)",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4,r0,r1,r2);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 15: // harmonic oscillator model 
      fFunc1 = new TF1(name,"x^2*(1+[0]*(x/[1])^2)*exp(-(x/[1])^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;      
    case 16: // from Nuclear Physics A150 (1970) 631--654;
      fFunc1 = new TF1(name,"x^2*(1-[3]*sin([5]*x)*exp(-([4]*x)**2)/([5]*x)+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,0.102,0.35,2.76);
      break;
    case 90: // via TGraph input
      if (!fInputDens) 
        fInputDens = new TGraph(GetTitle());
      assert(fInputDens!=0);
      fFunc1 = new TF1(name,[&](double*x, double *p){return p[0]*fInputDens->Eval(x[0])*x[0]*x[0];}, 0.0, fMaxR, 1);
      fFunc1->SetParameter(0,1);
      break;
    case 99: // Trajectum models
#ifndef USE_TRNUCGEN 
      cerr << "ERROR: Need libtrnucgen.so for trajectum models" << endl;
      gSystem->Exit(123);
#endif
      break;
    default:
      cerr << "Warning: Could not find function type " << fF << endl;
      return;
  }
  if (fFunc1)   
    fFunc1->SetNpx(1000);
  return;
}

void TGlauNucleus::SetA(Double_t ia, Double_t ia2)
{
  fA  = ia;
  fA2 = ia2;
  switch (fF) {
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 7:  // Box method
    case 15: // harmonic oscillator model 
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(1,fA);
      break;
    case 8:
      fFunc3->SetParameter(1,fA);
      break;
    case 11: //p&n
      fFunc1->SetParameter(1,fA);//proton
      fFunc2->SetParameter(1,fA2);//neutron
      break;
    default:
      cerr << "Warning: fA not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetBeta(Double_t b2, Double_t b4) 
{
  fBeta2=b2; 
  fBeta4=b4;      
  if (fFunc3) {
    fFunc3->SetParameter(2,fBeta2);
    fFunc3->SetParameter(3,fBeta4);
  }
}

void TGlauNucleus::SetBeta(Double_t b2, Double_t b3, Double_t b4, Double_t g) 
{
  fBeta2=b2; 
  fBeta3=b3; 
  fBeta4=b4; 
  fGamma=g;
  fF=7; 
}

void TGlauNucleus::SetGamma(Double_t g) 
{
  fGamma=g; 
  fF=7; 
}

void TGlauNucleus::SetR(Double_t ir, Double_t ir2)
{
  fR  = ir;
  fR2 = ir2;
  switch (fF) {
    case 0:  // Proton exp
    case 9:  // Proton gaus
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 15: // harmonic oscillator model 
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(0,fR);
      break;
    case 8:
      fFunc3->SetParameter(0,fR);
      break;
    case 10: // Proton
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // p&n
      fFunc1->SetParameter(0,fR);//proton
      fFunc2->SetParameter(0,fR2);//neutron
      break;
    default:
      cerr << "Warning: fR not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetW(Double_t iw)
{
  fW = iw;
  switch (fF) {
    case 1: // 3pF
    case 2: // 3pG
      fFunc1->SetParameter(2,fW);
      break;
    default:
      cerr << "Warning: fW not needed for function " << fF <<endl;
  }
}

Bool_t TGlauNucleus::TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const
{
  if (fMinDist<=0)
    return kTRUE;
  const Double_t md2 = fMinDist*fMinDist; 
  for (Int_t j = 0; j<n; ++j) {
    TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
    Double_t xo=other->GetX();
    Double_t yo=other->GetY();
    Double_t zo=other->GetZ();
    Double_t dist2 = (x-xo)*(x-xo)+
		                 (y-yo)*(y-yo)+
		                 (z-zo)*(z-zo);
    if (dist2<md2) {
      return kFALSE;
    }
  }
  return kTRUE;
}

TVector3 &TGlauNucleus::ThrowNucleons(Double_t xshift)
{
  AllocateNucleons();
  RandomizeNucleons();

  cmscheck: /* start over here in case shift was too large */

  fTrials = 0;
  fNonSmeared = 0;
  fPhiRot = gRandom->Rndm()*2*TMath::Pi();
  const Double_t cosThetaRot = 2*gRandom->Rndm()-1;
  fThetaRot = TMath::ACos(cosThetaRot);
  fXRot = gRandom->Rndm()*2*TMath::Pi();
  fYRot = gRandom->Rndm()*2*TMath::Pi();
  fZRot = gRandom->Rndm()*2*TMath::Pi();

  const Bool_t hulthen = (fF==3||fF==4);
  
  if (fF==6) { // select configuration read from file
    SelectFromNucArr();
  } else if (fF==99) { // use trnucgen to generate nucleons
#ifdef USE_TRNUCGEN
    static tr::NucleusGenerator ngen;
    ngen.setMinDistance(fMinDist);
    ngen.setRecenter(0); 
    int nuctype = ngen.getNucleusType(GetName());
    if (nuctype == -1) {
      cerr << "ERROR: Failed to generate nucleons for " << GetName() << " with trnucgen; check if it is supported by looking at the trnucgen/nucleusgenerator.h file" << endl;
      gSystem->Exit(123);
    }
    bool res = false;
    int trials = 0;
    while (!res) {
      if (++trials > 100) {
        cerr << "ERROR: Failed to generate nucleons for " << GetName() << " check if the minimal distance (" << fMinDist << " fm) is supported by looking at the trnucgen/nucleusgenerator.cpp file" << endl;
        gSystem->Exit(123);
      }
      res = ngen.generate(GetName());
      fWeight = ngen.getWeight();
      if (TMath::IsNaN(fWeight)) {
        cerr << "Warning: Weight is NaN" << endl;
        res = false;
      }
    }
    fTrials = trials;
    fN = ngen.getNumNucleons();
    fZ = ngen.getNumProtons();
    AllocateNucleons();
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      nucleon->SetXYZ(ngen.getNucleonPositionX(i),
                      ngen.getNucleonPositionY(i),
                      ngen.getNucleonPositionZ(i));
      if (ngen.getNucleonIsProton(i)) {
        nucleon->SetType(1);
      } else {
        nucleon->SetType(0);
      }
      //rotation is done in the nucleus generator
      //nucleon->RotateXYZ(fPhiRot,fThetaRot);
    }
    SetTitle(ngen.getNucleusDescription().c_str());
#else
    cerr << "ERROR: Should not end here because you do not have libtrnucgen.so" << endl;
#endif
  } else if (fN==1) { //special treatment for proton
    Double_t r = fFunc1->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(0));
    nucleon->Reset();
    nucleon->SetXYZ(r * stheta * TMath::Cos(phi),
		                r * stheta * TMath::Sin(phi),
		                r * ctheta);
    fTrials = 1;
  } else if (fN==2 && hulthen) { //special treatment for Hulten
    Double_t r = fFunc1->GetRandom()/2;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);

    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    nucleon1->Reset();
    nucleon1->SetXYZ(r * stheta * TMath::Cos(phi),
                     r * stheta * TMath::Sin(phi),
                     r * ctheta);
    nucleon2->Reset();
    if (fF==4) { // place opposite of 1
      nucleon2->SetXYZ(-nucleon1->GetX(),
		                   -nucleon1->GetY(),
		                   -nucleon1->GetZ());
    } else {
      r = fFunc1->GetRandom()/2;
      phi = gRandom->Rndm() * 2 * TMath::Pi();
      ctheta = 2*gRandom->Rndm() - 1;
      stheta = TMath::Sqrt(1-ctheta*ctheta);
      nucleon2->SetXYZ(r * stheta * TMath::Cos(phi),
		                   r * stheta * TMath::Sin(phi),
		                   r * ctheta);
    }
    fTrials = 1;
  } else { // all other nuclei 
    const Double_t startingEdge  = 20; // throw nucleons within a cube of this size (fm)
    const Double_t startingEdgeX = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeY = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeZ = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Int_t nslots = 2*startingEdge/fNodeDist+1;
    if (fNodeDist>0) {
      if (fMinDist>fNodeDist) {
        cout << "Minimum distance (nucleon hard core diameter) [" 
          << fMinDist << "] cannot be larger than the nodal spacing of the grid [" 
          << fNodeDist << "]." << endl;
        cout << "Quitting...." << endl;
        gSystem->Exit(123);
      }
      if (!fIsUsed)
        fIsUsed = new TBits(nslots*nslots*nslots);
      else
        fIsUsed->ResetAllBits();
    }
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while (1) {
        ++fTrials;
        Bool_t nucleon_inside = 0;
        Double_t x=999, xsmeared=999;
        Double_t y=999, ysmeared=999;
        Double_t z=999, zsmeared=999;
        if (fF==5||fF==7) { // the extended way, throw in a box and test the weight
          while (!nucleon_inside) {
            x = (fMaxR)*(gRandom->Rndm() * 2 - 1);
            y = (fMaxR)*(gRandom->Rndm() * 2 - 1);
            z = (fMaxR)*(gRandom->Rndm() * 2 - 1);
            Double_t r = TMath::Sqrt(x*x+y*y);
            Double_t theta = TMath::ATan2(r,z);
            Double_t R = TMath::Sqrt(x*x+y*y+z*z);
            Double_t Rtheta = fR;
            if (fF==5)
              Rtheta += fBeta2*TMath::Cos(theta)*TMath::Cos(theta);
            else if (fF==7) {
#ifdef HAVE_MATHMORE
              if (fBeta2!=0) {
                if (fGamma!=0) {
                  Double_t phi = TMath::Abs(y)/y*TMath::ACos(x/r);
                  Double_t v2 = TMath::Sqrt(15./64/TMath::Pi())*TMath::Power(TMath::Sin(theta),2)*TMath::Cos(2*phi);
                  Double_t b2 = TMath::Cos(fGamma)*ROOT::Math::sph_legendre(2,0,theta) + TMath::Sin(fGamma)*v2;
                  Rtheta += fR*fBeta2*b2;
                } else {
                  Rtheta += fR*fBeta2*ROOT::Math::sph_legendre(2,0,theta);
                }
              }
              if (fBeta3!=0) {
                Rtheta += fR*fBeta3*ROOT::Math::sph_legendre(3,0,theta);
              }
              if (fBeta4!=0) {
                Rtheta += fR*fBeta4*ROOT::Math::sph_legendre(4,0,theta);
              }
#else
              cerr << "Should not end here because you do not have libMathMore" << endl;
#endif
            }
            Double_t prob = (1+fW*TMath::Power(R/Rtheta,2))/(1+TMath::Exp((R-Rtheta)/fA));
            if (gRandom->Rndm()<prob) 
              nucleon_inside=1;
          }
        } else if ((fF==8) || (fF==14)) { // use TF2
          Double_t r;
          Double_t theta;
          fFunc3->GetRandom2(r,theta);
          Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
          x = r * TMath::Sin(phi) * TMath::Sin(theta);
          y = r * TMath::Cos(phi) * TMath::Sin(theta);
          z = r *                   TMath::Cos(theta);
        } else { // all other types
          TF1 *ff = fFunc1;
          if ((fFunc2) && (nucleon->GetType()==0) && (fF!=90))
            ff = fFunc2;
          if (fNodeDist<=0) { // "continuous" mode
            Double_t r = ff->GetRandom();
            Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
            Double_t ctheta = 2*gRandom->Rndm() - 1 ;
            Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
            x = r * stheta * TMath::Cos(phi);
            y = r * stheta * TMath::Sin(phi);
            z = r * ctheta;
          } else { // "grid/lattice" mode
            Int_t iNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t jNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t kNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t index=iNode*nslots*nslots+jNode*nslots+kNode;
            if (fIsUsed->TestBitNumber(index))
              continue;
            if (fLattice==1) {       // Primitive cubic system (PCS) -> https://en.wikipedia.org/wiki/Cubic_crystal_system
              x = fNodeDist*(iNode) - startingEdgeX;
              y = fNodeDist*(jNode) - startingEdgeY;
              z = fNodeDist*(kNode) - startingEdgeZ;
            } else if (fLattice==2) { //Body centered cubic (BCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(-iNode+jNode+kNode) - 0.5*startingEdgeX;
              y = 0.5*fNodeDist*(+iNode-jNode+kNode) - 0.5*startingEdgeY;
              z = 0.5*fNodeDist*(+iNode+jNode-kNode) - 0.5*startingEdgeZ;
            } else if (fLattice==3) { //Face Centered Cubic (FCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(jNode+kNode) - startingEdgeX;
              y = 0.5*fNodeDist*(iNode+kNode) - startingEdgeY;
              z = 0.5*fNodeDist*(iNode+jNode) - startingEdgeZ;
            } else {                  //Hexagonal close packing (HCP) -> https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
              x = 0.5*fNodeDist*(2*iNode+((jNode+kNode)%2))          - startingEdgeX;
              y = 0.5*fNodeDist*(TMath::Sqrt(3)*(jNode+(kNode%2)/3)) - startingEdgeY;
              z = 0.5*fNodeDist*(kNode*2*TMath::Sqrt(6)/3)           - startingEdgeZ;
            }
            const Double_t r2 = x*x + y*y + z*z;
            const Double_t r  = TMath::Sqrt(r2);
            if ((r>fMaxR)||(r2*gRandom->Rndm()>ff->Eval(r)))
              continue;
            if (fSmearing>0.0) {
              Int_t nAttemptsToSmear = 0;
              while (1) {
                xsmeared = x*gRandom->Gaus(1.0,fSmearing);
                ysmeared = y*gRandom->Gaus(1.0,fSmearing);
                zsmeared = z*gRandom->Gaus(1.0,fSmearing);
                nAttemptsToSmear++;
                if (TestMinDist(i,xsmeared,ysmeared,zsmeared)) {
                  x = xsmeared;
                  y = ysmeared;
                  z = zsmeared;
                  break;
                }
                if (nAttemptsToSmear>=99) {
                  cerr << "Could not place on this node :: [" << x <<","<< y <<","<< z <<"] r = " << TMath::Sqrt(x*x+y*y+z*z) << " fm; "
                    << "Node (" << iNode << "," << jNode << "," << kNode << ") not smeared !!!" << endl;
                  ++fNonSmeared;
                  break;
                }
              }
            }
            fIsUsed->SetBitNumber(index);
          } /* end "grid/lattice mode" */
        }
        nucleon->SetXYZ(x,y,z);
        if (fF==5||fF==7||fF==8||fF==14) 
          nucleon->RotateXYZ(fPhiRot,fThetaRot); // Uranium etc.
        if (fNodeDist>0) {
          nucleon->RotateXYZ_3D(fXRot,fYRot,fZRot);
          break;
        }
        if (TestMinDist(i,x,y,z))
          break;
      }
    }
  }    

  // calculate center of mass
  Double_t sumx=0;       
  Double_t sumy=0;       
  Double_t sumz=0;       
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    sumx += nucleon->GetX();
    sumy += nucleon->GetY();
    sumz += nucleon->GetZ();
  }
  sumx = sumx/fN;
  sumy = sumy/fN;
  sumz = sumz/fN;

  static TVector3 finalShift;
  finalShift.SetXYZ(sumx,sumy,sumz);
  if (finalShift.Mag()>fSmax)
    goto cmscheck;
  
  Double_t fsumx = 0;
  Double_t fsumy = 0;
  Double_t fsumz = 0;
  if (fRecenter==1) {
    fsumx = sumx;
    fsumy = sumy;
    fsumz = sumz;
  } else if (fRecenter==2) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(fN-1));
    Double_t x = nucleon->GetX() - fN*sumx;
    Double_t y = nucleon->GetY() - fN*sumy;
    Double_t z = nucleon->GetZ() - fN*sumz;
    nucleon->SetXYZ(x,y,z);
  } else if ((fRecenter==3)||(fRecenter==4)) {
    if (finalShift.Mag()>1e-3) {
      TVector3 zVec;
      zVec.SetXYZ(0,0,1);
      TVector3 shiftVec;
      shiftVec.SetXYZ(sumx,sumy,sumz);
      TVector3 orthVec;
      orthVec = shiftVec.Cross(zVec);
      TRotation myRot;
      myRot.Rotate(shiftVec.Angle(zVec),orthVec);
      TVector3 myNuc;
      for (Int_t i = 0; i<fN; ++i) {
        TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
        myNuc.SetXYZ(nucleon->GetX(),nucleon->GetY(),nucleon->GetZ());
        myNuc.Transform(myRot);
        nucleon->SetXYZ(myNuc.X(), myNuc.Y(), myNuc.Z());
      }
      if (fRecenter==3)
        fsumz = shiftVec.Mag();
    }
  } else if (fRecenter==5) {
    fsumx = sumx;
    fsumy = sumy;
  }

  // recenter and shift
  sumx=0;       
  sumy=0;       
  sumz=0;       
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    nucleon->SetXYZ(nucleon->GetX()-fsumx + xshift,
		                nucleon->GetY()-fsumy,
		                nucleon->GetZ()-fsumz);
    sumx += nucleon->GetX();
    sumy += nucleon->GetY();
    sumz += nucleon->GetZ();
  }
  sumx = sumx/fN;
  sumy = sumy/fN;
  sumz = sumz/fN;
  finalShift.SetXYZ(sumx,sumy,sumz);
  return finalShift;
}

//---------------------------------------------------------------------------------
// ClassImp(TGlauberMC)
// ClassImp(TGlauberMC::Event)
Double_t TGlauberMC::gDefOmega=-1;
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect, Double_t xsectsigma, Double_t xsectnp) :
  fANucleus(NA),fBNucleus(NB),
  fXSect(xsect),fXSectNP(xsectnp),fXSectOmega(0),fXSectLambda(0),fXSectEvent(0),
  fNucleonsA(0),fNucleonsB(0),fNucleons(0),
  fAN(0),fBN(0),fNt(0),
  fEvents(0),fTotalEvents(0),fBmin(0),fBmax(20),fHardFrac(0.65),
  fDetail(99),fCalcArea(0),fCalcLength(0), fDoCore(0), fDoAAGG(1),
  fSigH(129.4), fShadow(0), fOmega(gDefOmega),
  fMaxNpartFound(0),f2Cx(0),fPTot(0),fNNProf(0),
  fEv()
{
  if (xsect<0) { // negative cross section, means energy is provided
    Double_t energy = -xsect;
    fSigH = getSigmaHard(energy);
    if (energy < 10) { // low energy, use Bystricky parameterization
      fXSect = getSigmaPP_Bystricky(energy);
      if (xsectnp<0)
        fXSectNP = getSigmaNP_Bystricky(energy);
      cout << "Using sigma_NN=" << fXSect << " mb and sigma_NP=" << fXSectNP << " mb for energy=" << energy << " GeV" << endl;
    } else {
      fXSect = getSigmaNN(energy);
      fXSectNP = 0;
      cout << "Using sigma_NN=" << fXSect << " mb and sigmaHard= " << fSigH << " mb for energy=" << energy << " GeV" << endl;
    }
  }
  if (xsectsigma>0) {
    fXSectOmega = xsectsigma;
    fXSectLambda = 1;
    fPTot = new TF1("fPTot","((x/[2])/(x/[2]+[0]))*exp(-(((x/[2])/[0]-1 )**2)/([1]*[1]))/[2]",0,300);
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    fPTot->SetNpx(1000);
    fXSectLambda = fXSect/fPTot->GetHistogram()->GetMean();
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    Double_t mean = fPTot->GetHistogram()->GetMean();
    cout << "Using fluctuating cross section with <sigma>=" << mean << "mb, using fXSectOmega=" << fXSectOmega << " and lambda=" << fXSectLambda << endl;
  }

  TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("Glauber %s+%s Version %s",fANucleus.GetName(),fBNucleus.GetName(),Version()));
  SetName(name);
  SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
  // calc next event
  if (!fNucleonsA) {
    fNucleonsA = fANucleus.GetNucleons();
    fAN = fANucleus.GetN();
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
    }
  }

  if (!fNucleonsB) {
    fNucleonsB = fBNucleus.GetNucleons();
    fBN = fBNucleus.GetN();
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
    }
  }

  Double_t xsecA[999] = {0};
  Double_t xsecB[999] = {0};
  if (fPTot) {
    fXSectEvent = fPTot->GetRandom();
    if (fDoAAGG) {
      for (Int_t i = 0; i<fAN; ++i)
	      xsecA[i] = fPTot->GetRandom();
      for (Int_t i = 0; i<fBN; ++i)
	      xsecB[i] = fPTot->GetRandom();
    }
  } else 
    fXSectEvent = fXSect;

  if (fOmega>0 && !fNNProf) {
    if (fOmega<2) {
      cout << "Setting up NN profile for omega=" << fOmega << endl;
      fNNProf = getNNProf(fXSect, fOmega);
    } else if (fOmega==7) {
      cout << "Setting up NN profile for HIJING (omega==7)" << endl;
      fNNProf = getNNHijing(fXSect);
    } else if (fOmega==8) {
      cout << "Setting up NN profile for PYTHIA (omega==8)" << endl;
      fNNProf = getNNPythia(fXSect);
    } else if ((fOmega>=9)&&(fOmega<11)) {
      Double_t w = fOmega-9;
      cout << "Setting up NN profile for Trento (w=" << w << ")" << endl;
      fNNProf = getNNTrento(fXSect,w);
    } else {
      cerr << "Value of omega is not supported: " << fOmega << endl;
      gSystem->Exit(1);
    }
  } else if (fOmega<0) {
    cout << "Using hard-sphere approximation (default);\n use SetOmega(omega) or SetDefOmega(omega) in rootlogon.C if you want NN overlap" << endl;
    fOmega = 0;
  }

  // "ball" diameter = distance at which two balls interact
  Double_t d2pp = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
  Double_t d2np = (Double_t)fXSectNP/(TMath::Pi()*10);    // in fm^2
  Double_t d2   = d2pp;
  Double_t bh   = TMath::Sqrt(d2*fHardFrac);
  if (fNNProf) {
    Double_t xmin=0,xmax=0;
    fNNProf->GetRange(xmin,xmax);
    d2 = xmax*xmax;
  }

  Double_t sigS = fXSectEvent; // 57 mb for HIJING;
  const Double_t b02  = 0.5 * sigS * 0.1 / TMath::Pi();
  Double_t sigHS = fSigH;
  if (fShadow) {
    const Double_t rrb   = TMath::Min(1., bgen * bgen / 35.2 / 1.44);
    const Double_t aphx  = 0.1 * 4./3. * 4.92 * TMath::Sqrt(1. - rrb);
    sigHS = fSigH - aphx * 103.65;
  }

  fEv.Reset();
  memset(fBC,0,sizeof(Bool_t)*999*999);
  Int_t nc=0,nh=0,njet=0;
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Bool_t tB=nucleonB->GetType();
    for (Int_t j = 0; j<fAN; ++j) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      Bool_t tA=nucleonA->GetType();
      if (fXSectNP>0) {
        if (tA!=tB) d2 = d2np;
        else        d2 = d2pp;
      } else if (fDoAAGG && fPTot){
        d2 = 0.5*(xsecA[j]+xsecB[i])/(TMath::Pi()*10);
      }
      if (dij>d2) 
        continue;
      Double_t bij = TMath::Sqrt(dij);
      if (fNNProf) {
        Double_t val = fNNProf->Eval(bij);
        Double_t ran = gRandom->Uniform();
        if (ran>val)
          continue;
        //from HIJING
        Double_t tt = 2. * val * sigHS/sigS;
        Double_t xr = - TMath::Log(TMath::Exp(-tt) + gRandom->Rndm()*(1.-TMath::Exp(-tt)));
        while (1) {
          ++njet;
          xr -= TMath::Log(gRandom->Rndm());
          if (xr > tt)
            break;
        }
      }
      fEv.Nmpi = njet;
      nucleonB->Collide();
      nucleonA->Collide();
      fBC[i][j] = 1;
      fEv.BNN  += bij;
      ++nc;
      if (bij<bh)
        ++nh;
      if (tA!=tB)
        ++fEv.Ncollpn;
      else if (tA==1)
        ++fEv.Ncollpp;
      else
        ++fEv.Ncollnn;
      if (nc==1) {
        fEv.X0 = (nucleonA->GetX()+nucleonB->GetX())/2;
        fEv.Y0 = (nucleonA->GetY()+nucleonB->GetY())/2;
      }
    }
  }
  fEv.B = bgen;
  Double_t w1 = fANucleus.GetWeight();
  Double_t w2 = fBNucleus.GetWeight();
  if (w1==0) w1 = 1;
  if (w2==0) w2 = 1;
  Double_t w = w1*w2;
  if (TMath::Abs(w)>1e6) {
    cout << "Warning: Weight is too large: " << w << endl;
    cout << "w1 = " << w1 << ", w2 = " << w2 << endl;
    cout << "fANucleus.GetName() = " << fANucleus.GetName() << endl;
    cout << "fBNucleus.GetName() = " << fBNucleus.GetName() << endl;
    gSystem->Exit(1);
  }
  fTotalEvents += w;
  if (nc>0) {
    fEvents += w;
    fEv.Weight    = w;
    fEv.Ncoll     = nc;
    fEv.Nhard     = nh;
    fEv.BNN      /= nc;
    return CalcResults(bgen);
  }
  return kFALSE;
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  Double_t sumW=0;
  Double_t sumWA=0;
  Double_t sumWB=0;

  Double_t sinphi[10] = {0};
  Double_t cosphi[10] = {0};
  Double_t rn[10]     = {0};

  const Int_t kNc = fDoCore; // used later for core/corona

  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fEv.MeanXSystem  += xA;
    fEv.MeanYSystem  += yA;
    fEv.MeanXA  += xA;
    fEv.MeanYA  += yA;
    if (nucleonA->IsWounded()) {
      Double_t w = nucleonA->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonA->IsNeutron())
        ++fEv.NpartAn;
      if (nucleonA->GetNColl()==1) {
	      ++fEv.Npart0;
        if (nucleonA->IsNeutron())
          ++fEv.Npart0n;
      }
      ++fEv.NpartA;
      sumW   += w;
      sumWA  += w;
      fEv.MeanX  += xA * w;
      fEv.MeanY  += yA * w;
      fEv.MeanX2 += xA * xA * w;
      fEv.MeanY2 += yA * yA * w;
      fEv.MeanXY += xA * yA * w;
    } else { // count spectator neutrons
      if (nucleonA->IsNeutron())
        ++fEv.SpecA;
    }
  }

  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fEv.MeanXSystem  += xB;
    fEv.MeanYSystem  += yB;
    fEv.MeanXB  += xB;
    fEv.MeanYB  += yB;
    if (nucleonB->IsWounded()) {
      Double_t w = nucleonB->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonB->IsNeutron())
        ++fEv.NpartBn;
      if (nucleonB->GetNColl()==1) {
	      ++fEv.Npart0;
        if (nucleonB->IsNeutron())
          ++fEv.Npart0n;
      }
      ++fEv.NpartB;
      sumW   += w;
      sumWB  += w;
      fEv.MeanX  += xB * w;
      fEv.MeanY  += yB * w;
      fEv.MeanX2 += xB * xB * w;
      fEv.MeanY2 += yB * yB * w;
      fEv.MeanXY += xB * yB * w;
    } else { // count spectator neutrons
      if (nucleonB->IsNeutron())
        ++fEv.SpecB;
    }
  }
  if (fEv.Npart>0) {
    fEv.MeanX  /= sumW;
    fEv.MeanY  /= sumW;
    fEv.MeanX2 /= sumW;
    fEv.MeanY2 /= sumW;
    fEv.MeanXY /= sumW;
  } else {
    fEv.MeanX = 0;
    fEv.MeanY  = 0;
    fEv.MeanX2 = 0;
    fEv.MeanY2 = 0;
    fEv.MeanXY = 0;
  }

  if (fAN+fBN>0) {
    fEv.MeanXSystem /= (fAN + fBN);
    fEv.MeanYSystem /= (fAN + fBN);
  } else {
    fEv.MeanXSystem = 0;
    fEv.MeanYSystem = 0;
  }
  if (fAN>0) {
    fEv.MeanXA /= fAN;
    fEv.MeanYA /= fAN;
  } else {
    fEv.MeanXA = 0;
    fEv.MeanYA = 0;
  }
  if (fBN>0) {
    fEv.MeanXB /= fBN;
    fEv.MeanYB /= fBN;
  } else {
    fEv.MeanXB = 0;
    fEv.MeanYB = 0;
  }

  fEv.VarX  = fEv.MeanX2-(fEv.MeanX*fEv.MeanX);
  fEv.VarY  = fEv.MeanY2-(fEv.MeanY*fEv.MeanY);
  fEv.VarXY = fEv.MeanXY-fEv.MeanX*fEv.MeanY;
  Double_t tmpa = fEv.VarX*fEv.VarY-fEv.VarXY*fEv.VarXY;
  if (tmpa<0) 
    fEv.AreaW = -1;
  else 
    fEv.AreaW = TMath::Sqrt(tmpa);

  if (fEv.Npart>0) {
    // do full moments relative to meanX and meanY
    for (Int_t n = 1; n<10; ++n) {
      for (Int_t ia = 0; ia<fAN; ++ia) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(ia));
        if (nucleonA->GetNColl()<=kNc) 
          continue;
        Double_t xA=nucleonA->GetX() - fEv.MeanX;
        Double_t yA=nucleonA->GetY() - fEv.MeanY;
        Double_t r = TMath::Sqrt(xA*xA+yA*yA);
        Double_t phi = TMath::ATan2(yA,xA);
        Double_t w = n;
        if (n==1) 
          w = 3; // use r^3 weighting for Ecc1/Psi1
        Double_t rw = TMath::Power(r,w);
        cosphi[n] += rw*TMath::Cos(n*phi);
        sinphi[n] += rw*TMath::Sin(n*phi);
        rn[n] += rw;
      }
      for (Int_t ib = 0; ib<fBN; ++ib) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(ib));
        if (nucleonB->GetNColl()<=kNc) 
          continue;
        Double_t xB=nucleonB->GetX() - fEv.MeanX;
        Double_t yB=nucleonB->GetY() - fEv.MeanY;
        Double_t r = TMath::Sqrt(xB*xB+yB*yB);
        Double_t phi = TMath::ATan2(yB,xB);
        Double_t w = n;
        if (n==1)
          w = 3; // use r^3 weighting for Ecc1/Psi1
        Double_t rw = TMath::Power(r,w);
        cosphi[n] += rw*TMath::Cos(n*phi);
        sinphi[n] += rw*TMath::Sin(n*phi);
        rn[n] += rw;
      }
      cosphi[n] /= fEv.Npart;
      sinphi[n] /= fEv.Npart;
      rn[n] /= fEv.Npart;
      if (rn[n]>0) {
	      fPsiN[n] = (TMath::ATan2(sinphi[n],cosphi[n]) + TMath::Pi())/n;
	      fEccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
      } else {
	      fPsiN[n] = -1;
	      fEccN[n] = -1;
      }
    }
    if (!kNc) { //silly test but useful to catch errors 
      Double_t t=TMath::Sqrt(TMath::Power(fEv.VarY-fEv.VarX,2)+4.*fEv.VarXY*fEv.VarXY)/(fEv.VarY+fEv.VarX)/fEccN[2];
      if (t<0.99||t>1.01)
        cout << "Warning: Expected t=1 but found t=" << t << endl;
    }
  }

  fEv.B      = bgen;
  fEv.PhiA   = fANucleus.GetPhiRot();
  fEv.ThetaA = fANucleus.GetThetaRot();
  fEv.PhiB   = fBNucleus.GetPhiRot();
  fEv.ThetaB = fBNucleus.GetThetaRot();
  fEv.Psi1   = fPsiN[1];
  fEv.Ecc1   = fEccN[1];
  fEv.Psi2   = fPsiN[2];
  fEv.Ecc2   = fEccN[2];
  fEv.Psi3   = fPsiN[3];
  fEv.Ecc3   = fEccN[3];
  fEv.Psi4   = fPsiN[4];
  fEv.Ecc4   = fEccN[4];
  fEv.Psi5   = fPsiN[5];
  fEv.Ecc5   = fEccN[5];

  if (fCalcArea) {
    const Int_t nbins=200;
    const Double_t ell=10;
    const Double_t da=2*ell*2*ell/nbins/nbins;
    const Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
    const Double_t r2 = d2/4.;
    const Double_t mx = fEv.MeanX;
    const Double_t my = fEv.MeanY;
    TH2D areaA("hAreaA",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    TH2D areaB("hAreaB",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      if (!nucleonA->IsWounded())
        continue;
      if (nucleonA->GetNColl()==kNc)
        continue;
      Double_t x = nucleonA->GetX()-mx;
      Double_t y = nucleonA->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaA.GetBin(xi,yi);
          Double_t val=areaA.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaA.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaA.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaA.SetBinContent(bin,1);
        }
      }
    }
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      if (!nucleonB->IsWounded())
        continue;
      if (nucleonB->GetNColl()==kNc)
        continue;
      Double_t x = nucleonB->GetX()-mx;
      Double_t y = nucleonB->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaB.GetBin(xi,yi);
          Double_t val=areaB.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaB.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaB.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaB.SetBinContent(bin,1);
        }
      }
    }
    Double_t overlap1=0;
    Double_t overlap2=0;
    for (Int_t xi=1; xi<=nbins; ++xi) {
      for (Int_t yi=1; yi<=nbins; ++yi) {
        Int_t bin = areaA.GetBin(xi,yi);
        Double_t vA=areaA.GetBinContent(bin);
        Double_t vB=areaB.GetBinContent(bin);
        if (vA>0&&vB>0)
          ++overlap1;
        if (vA>0||vB>0)
          ++overlap2;
      }
    }
    fEv.AreaO = overlap1*da;
    fEv.AreaA = overlap2*da;
  }

  if (fCalcLength) {
    const Double_t krhs = TMath::Sqrt(fXSectEvent/40./TMath::Pi());
    const Double_t ksg  = krhs/TMath::Sqrt(5);
    const Double_t kDL  = 0.1;
    TF1 rad("rad","2*pi/[0]/[0]*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,5*ksg); 
    rad.SetParameter(0,ksg);
    const Double_t minval = rad.Eval(5*ksg);
    fEv.Phi0         = gRandom->Uniform(0,TMath::TwoPi());
    Double_t kcphi0  = TMath::Cos(fEv.Phi0);
    Double_t ksphi0  = TMath::Sin(fEv.Phi0);
    Double_t x       = fEv.X0;
    Double_t y       = fEv.Y0;
    Double_t i0a     = 0;
    Double_t i1a     = 0;
    Double_t l       = 0;
    Double_t val     = CalcDens(rad,x,y);
    while (val>minval) {
      x     += kDL * kcphi0;
      y     += kDL * ksphi0;
      i0a   += val;
      i1a   += l*val;
      l+=kDL;
      val    = CalcDens(rad,x,y);
    }
    fEv.Length = 2*i1a/i0a;
  }

  if (fEv.Npart > fMaxNpartFound) 
    fMaxNpartFound = fEv.Npart;

  return kTRUE;
}

Double_t TGlauberMC::CalcDens(TF1 &prof, Double_t xval, Double_t yval) const
{
  Double_t rmin=0,rmax=0;
  prof.GetRange(rmin,rmax);
  Double_t r2max = rmax*rmax;
  Double_t ret = 0;
  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    if (!nucleonA->IsWounded())
      continue;
    Double_t x = nucleonA->GetX();
    Double_t y = nucleonA->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    if (!nucleonB->IsWounded())
      continue;
    Double_t x = nucleonB->GetX();
    Double_t y = nucleonB->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  return ret;
}

void TGlauberMC::Draw(Option_t* option)
{
  static TH2F *h2f = new TH2F("hGlauberMC",";x (fm);y(fm)",1,-18,18,1,-12,12);

  h2f->Reset();
  h2f->SetStats(0);
  h2f->Draw();

  TEllipse e;
  e.SetFillColor(0);
  e.SetFillStyle(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
  fANucleus.Draw(fXSect, kMagenta, kYellow);
  fBNucleus.Draw(fXSect, kMagenta, kOrange);

  TString opt(option);
  if (opt.IsNull())
    return;

  Double_t sy2 = GetSy2();
  Double_t sx2 = GetSx2();
  Double_t phase = 0;
  if (sy2<sx2) {
    Double_t d = sx2;
    sx2 = sy2;
    sy2 = d;
    phase = TMath::Pi()/2.;
  }
  Double_t x1 = (0.5*(sy2-sx2)+TMath::Sqrt(TMath::Power(sy2-sx2,2.)-4*TMath::Power(GetSxy(),2)));
  Double_t ang = TMath::ATan2(-GetSxy(),x1)+phase;
  TLine l;
  l.SetLineWidth(3);
  l.DrawLine(-10*TMath::Cos(ang),-10*TMath::Sin(ang),10*TMath::Cos(ang),10*TMath::Sin(ang));
}

Double_t TGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBmax*fBmax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
  if (!fNucleonsA || !fNucleonsB) return 0;
  if (fNucleons) return fNucleons;

  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; ++i) {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    allnucleons->Add(fNucleonsB->At(i));
  }
  fNucleons = allnucleons;
  return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  if (bgen<0) 
    bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

  fANucleus.ThrowNucleons(-bgen/2.);
  Double_t wA = fANucleus.GetWeight();
  fBNucleus.ThrowNucleons(bgen/2.);
  Double_t wB = fBNucleus.GetWeight();
  return CalcEvent(bgen);
}

Bool_t TGlauberMC::ReadNextEvent(Bool_t calc, const char *fname)
{
  static TFile *inf = 0;
  static Int_t iev  = 0;
  if (fname) {
    cout << "ReadNextEvent: Setting up file " << fname << endl;
    delete inf;
    inf = TFile::Open(fname);
    if (!inf) 
      return 0;
    if (!fNucleonsA) {
      fANucleus.ThrowNucleons();
      fNucleonsA = fANucleus.GetNucleons();
      fAN = fANucleus.GetN();
      for (Int_t i = 0; i<fAN; ++i) {
	      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
	      nucleonA->SetInNucleusA();
      }
    }
    if (!fNucleonsB) {
      fBNucleus.ThrowNucleons();
      fNucleonsB = fBNucleus.GetNucleons();
      fBN = fBNucleus.GetN();
      for (Int_t i = 0; i<fBN; ++i) {
	      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
	      nucleonB->SetInNucleusB();
      }
    }
    if (calc)
      return 1;
    fNt = dynamic_cast<TNtuple*>(inf->Get(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName())));
    if (!fNt) {
      cerr << "ReadNextEvent: Could not find ntuple!" << endl;
      inf->ls();
      return 0;
    }
    fNt->SetBranchAddress("Npart",&fEv.Npart);
    fNt->SetBranchAddress("Ncoll",&fEv.Ncoll);
    fNt->SetBranchAddress("B",&fEv.B);
    fNt->SetBranchAddress("BNN",&fEv.BNN);
    fNt->SetBranchAddress("VarX",&fEv.VarX);
    fNt->SetBranchAddress("VarY",&fEv.VarY);
    fNt->SetBranchAddress("VarXY",&fEv.VarXY);
    fNt->SetBranchAddress("NpartA",&fEv.NpartA);
    fNt->SetBranchAddress("NpartB",&fEv.NpartB);
    fNt->SetBranchAddress("Npart0",&fEv.Npart0);
    fNt->SetBranchAddress("Psi1",&fEv.Psi1);
    fNt->SetBranchAddress("Ecc1",&fEv.Ecc1);
    fNt->SetBranchAddress("Psi2",&fEv.Psi2);
    fNt->SetBranchAddress("Ecc2",&fEv.Ecc2);
    fNt->SetBranchAddress("Psi3",&fEv.Psi3);
    fNt->SetBranchAddress("Ecc3",&fEv.Ecc3);
    fNt->SetBranchAddress("Psi4",&fEv.Psi4);
    fNt->SetBranchAddress("Ecc4",&fEv.Ecc4);
    fNt->SetBranchAddress("Psi5",&fEv.Psi5);
    fNt->SetBranchAddress("Ecc5",&fEv.Ecc5);
    return 1;
  }
  if ((!inf)||(!fNt&&!calc)) {
    cerr << "ReadNextEvent was not initialized" <<endl;
    return 0;
  }
  TObjArray *arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  if (!arr) {
    if (iev==0) {
      cerr << "ReadNextEvent could not read nucleon array for event " << iev << endl;
      return 0;
    }
    iev = 0;
    cerr << "ReadNextEvent resetting to first event" << endl;
    arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  }

  Double_t bgenA=0, bgenB=0;
  Int_t inA=0, inB=0;
  const Int_t nNucls = arr->GetEntries();
  for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
    TGlauNucleon *nuclinp = static_cast<TGlauNucleon*>(arr->At(iNucl));
    TGlauNucleon *nuclout = 0;
    if (nuclinp->IsInNucleusB()) { 
      nuclout = static_cast<TGlauNucleon*>(fNucleonsB->At(inB));
      bgenB += nuclinp->GetX();
      ++inB;
    } else {
      nuclout = static_cast<TGlauNucleon*>(fNucleonsA->At(inA));
      bgenA += nuclinp->GetX();
      ++inA;
    }
    nuclout->Reset();
    nuclout->SetXYZ(nuclinp->GetX(),nuclinp->GetY(),nuclinp->GetZ());
    nuclout->SetType(nuclinp->GetType());
    nuclout->SetEnergy(nuclinp->GetEnergy());
    if (!calc)
      nuclout->SetNColl(nuclinp->GetNColl());
  }
  delete arr;
  Double_t bgen = bgenB/inB-bgenA/inA;
  if (calc) {
    Bool_t ret = CalcEvent(bgen);
    if (0) 
      cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
    ++iev;
    return ret;
  }
  Int_t ret = fNt->GetEntry(iev);
  if (ret<=0) 
    return 0;
  fEccN[1]=fEv.Ecc1;
  fEccN[2]=fEv.Ecc2;
  fEccN[3]=fEv.Ecc3;
  fEccN[4]=fEv.Ecc4;
  fEccN[5]=fEv.Ecc5;
  if (0) 
    cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
  if (0) { // test ntuple values vs re-calculated values
    Double_t npart = fEv.Npart;
    Double_t ncoll = fEv.Ncoll;
    Double_t ecc2  = fEv.Ecc2;
    CalcEvent(bgen);
    if (npart!=fEv.Npart) 
      cout << iev << " differ in npart " << npart << " " << fEv.Npart << endl;
    if (ncoll!=fEv.Ncoll) 
      cout << iev << " differ in ncoll " << ncoll << " " << fEv.Ncoll << endl;
    if (TMath::Abs(ecc2-fEv.Ecc2)>0.001) 
      cout << iev << " differ in ecc2 " << ecc2 << " " << fEv.Ecc2 << endl;
  }
  ++iev;
  return 1;
}

void TGlauberMC::Run(Int_t nevents, Double_t b)
{
  if (fNt == 0) {
    TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
    TString vars("Npart:Ncoll:Nhard:Nmpi:B:BNN:Ncollpp:Ncollpn:Ncollnn:VarX:VarY:VarXY:NpartA:NpartB:Npart0:NpartAn:NpartBn:Npart0n:AreaW:SpecA:SpecB:Weight");
    if (fDetail>1)
      vars+=":Psi1:Ecc1:Psi2:Ecc2:Psi3:Ecc3:Psi4:Ecc4:Psi5:Ecc5";
    if (fDetail>2)
      vars+=":AreaO:AreaA:X0:Y0:Phi0:Length";
    if (fDetail>3)
      vars+=":MeanX:MeanY:MeanX2:MeanY2:MeanXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB";
    if (fDetail>4)
      vars+=":PhiA:ThetaA:PhiB:ThetaB";
    fNt = new TNtuple(name,name,vars);
    fNt->SetDirectory(0);
    TObjArray *l = fNt->GetListOfBranches();
    for (Int_t i=0; i<l->GetEntries(); ++i) {
      TBranch *br = dynamic_cast<TBranch*>(l->At(i));
      if (br)
        br->SetCompressionLevel(9);
    }
  }

  for (Int_t i = 0; i<nevents; ++i) {
    while (!NextEvent(b)) {}
    fNt->Fill((Float_t*)(&fEv.Npart));
    if ((i>0)&&(i%100)==0) 
      cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;
  }

  TString title(Form("%s+%s, x-sect=%.1f mb, str=%s, totx-sect=%.2f mb",fANucleus.GetName(),fBNucleus.GetName(),fXSect,Str(),GetTotXSect()*1e3));
  fNt->SetTitle(title);
  if (nevents>99)
    cout << endl << "Done!" << endl;
}

void TGlauberMC::SetOmega(Double_t omega)
{
  if (omega<0) {
    cout << "Using hard-sphere approximation (default)" << endl;
    cout << " If you want NN overlap, use SetOmega(o) where " << endl;
    cout << "  - Gamma dist (0<o<2)" << endl;
    cout << "  - o=7: HIJING" << endl;
    cout << "   -o=8: PYTHIA" << endl;
    cout << "   -Trento (9<o<11), translates into d=o-9" << endl;
    fOmega = 0;
    return;
  }
  fOmega = omega;
}
