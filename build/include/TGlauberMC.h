#ifndef T_GLAUBER_MC_H
#define T_GLAUBER_MC_H

/*
 $Id: runglauber_v3.3.C 199 2025-08-30 16:26:30Z loizides $
 -------------------------------------------------------------------------------------
 Latest documentation: https://arxiv.org/abs/2507.05853
 -------------------------------------------------------------------------------------
 To run the code, you need to have the ROOT v6 (http://root.cern.ch/drupal/)
 environment. It is recommended to use the provided rootlogon.C script. 
 Else at the root prompt, enter
 root [0] gSystem->Load("libMathMore")
 root [1] .L runglauber_X.Y.C+
 (where X.Y denotes the version number).
 If you do not have libMathMore comment out "#define HAVE_MATHMORE" below.
 See the documentation for more information.
 -------------------------------------------------------------------------------------
 v3.3.1: Added possibility to set default value for omega in rootlogon.C (via 
  TGlauberMC::SetDefOmega); increased ranges for O and Ne nuclei, added Osat as
  replacement for Opar2
 -------------------------------------------------------------------------------------
 v3.3: Interface to TrNucGen (https://trnucgen.web.cern.ch), and support to read 
  nucleus configurations from text files, more pre-defined nucleus configurations
  (see Lookup function), calculations using TRENTO, HIJING and PYTHIA based NN profiles 
  (default is still HS approximation) and calculation of the number of multiple parton 
  interactions (MPI), possibity to get sigmaNN from parameterized fit (getSigmaNN),
  use of different NN and NP crossections at Hades energies, 
  several improvements in TGlauNucleus.
  see https://arxiv.org/abs/2507.05853
 -------------------------------------------------------------------------------------
 v3.2: Incorporates changes from v2.7, see https://arxiv.org/abs/1710.07098v3
 -------------------------------------------------------------------------------------
 v3.1:
  Fixes related to spherical nuclei, as well as consistent set of reweighted profiles 
  for Cu, Au and Xe, see https://arxiv.org/abs/1710.07098v2
 -------------------------------------------------------------------------------------
 v3.0:
  Major update to include separate profile for protons and neutrons, placement of nucleon 
  dof on lattice, as well as reweighted profiles for recentering, 
  see https://arxiv.org/abs/1710.07098v1
 -------------------------------------------------------------------------------------
 v2.7:
  New macro "runAndOutputLemonTree" for IP-Jazma input (1808.01276), as well as nucleon 
  configurations for He4, C, and O from wavefunction calculations, clarified use of Hulthen
  for deuteron, harmonic oscillator param for O, and new mode to use GlauberGribov also 
  in AA (enable with SetCalcAAGG)
 -------------------------------------------------------------------------------------
 v2.6:
  Includes runAndCalcDens macro, as well as definition for Al, and fixes beta4 for Si2,
  see https://arxiv.org/abs/1408.2549v8
 -------------------------------------------------------------------------------------
 v2.5:
  Include core/corona determination in Npart, and if requested for area from mc and eccentricity,
  as well as various Xe parameterizations including deformation,
  see https://arxiv.org/abs/1408.2549v7
 -------------------------------------------------------------------------------------
 v2.4: 
  Minor update to include Xenon and fix of the TGlauberMC::Draw function, 
  see https://arxiv.org/abs/1408.2549v4
 -------------------------------------------------------------------------------------
 v2.3: 
  Small bugfixes, see https://arxiv.org/abs/1408.2549v3
 -------------------------------------------------------------------------------------
 v2.2:
  Minor update to provide higher harmonic eccentricities up to n=5, and the average
  nucleon--nucleon impact parameter (bNN) in tree output. 
 -------------------------------------------------------------------------------------
 v2.1: 
  Minor update to include more proton pdfs, see https://arxiv.org/abs/1408.2549v2
 -------------------------------------------------------------------------------------
 v2.0: 
  First major update with inclusion of Tritium, Helium-3, and Uranium, as well as the 
  treatment of deformed nuclei and Glauber-Gribov fluctuations of the proton in p+A 
  collisions, see https://arxiv.org/abs/1408.2549v1
 -------------------------------------------------------------------------------------
 v1.1: 
  First public release of the PHOBOS MC Glauber, see https://arxiv.org/abs/0805.4411
 -------------------------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>
*/

//#define HAVE_MATHMORE

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TBits.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#ifdef HAVE_MATHMORE
#include <Math/SpecFuncMathMore.h>
#endif
#ifdef USE_TRNUCGEN
#define floatingpoint double
#include "nucleusgenerator.h"
#else
//#warning "Need libtrnucgen.so for trajectum models"
#endif
using namespace std;
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_ 3
#endif

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t signn=68.0, Double_t omega=0.6, Double_t G=1);
TF1 *getNNProfDist(Double_t signn=68.0, Double_t omega=0.6, Double_t G=1);
TF1 *getNNHijing(Double_t signn=68.0, Double_t mu=3.9);
TF1 *getNNHijingDist(Double_t signn=68.0, Double_t mu=3.9);
TF1 *getNNPythia(Double_t signn=68.0, Double_t m=1.85, Double_t rp=1);
TF1 *getNNPythiaDist(Double_t signn=68.0, Double_t m=1.85, Double_t rp=1);
TF1 *getNNTrento(Double_t signn=68.0, Double_t w=0.5) ;
TF1 *getNNTrentoDist(Double_t signn=68.0, Double_t w=0.5);
TF1 *getSigmaNNvsEnergy();
Double_t getSigmaNN(Double_t energy=5360);
TGraph *getSigmaNPvsEnergy_Bystricky();
Double_t getSigmaNP_Bystricky(Double_t energy);
TGraph *getSigmaPPvsEnergy_Bystricky();
Double_t getSigmaPP_Bystricky(Double_t energy);
TGraph *getSigmaHardvsEnergy();
Double_t getSigmaHard(Double_t energy=5360);

//---------------------------------------------------------------------------------
class TGlauNucleon : public TObject
{
  protected:
    Double32_t fX;            //Position of nucleon
    Double32_t fY;            //Position of nucleon
    Double32_t fZ;            //Position of nucleon
    Int_t      fType;         //0 = neutron, 1 = proton
    Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
    Int_t      fNColl;        //Number of binary collisions
    Double32_t fEn;           //Energy
  public:
    TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0), fEn(0) {}
    virtual   ~TGlauNucleon() {}
    void       Collide()                                  {++fNColl;}
    Double_t   Get2CWeight(Double_t x) const              {return 2.*(0.5*(1-x)+0.5*x*fNColl);}
    Double_t   GetEnergy()             const              {return fEn;}
    Int_t      GetNColl()              const              {return fNColl;}
    Int_t      GetType()               const              {return fType;}
    Double_t   GetX()                  const              {return fX;}
    Double_t   GetY()                  const              {return fY;}
    Double_t   GetZ()                  const              {return fZ;}
    Bool_t     IsNeutron()             const              {return (fType==0);}
    Bool_t     IsInNucleusA()          const              {return fInNucleusA;}
    Bool_t     IsInNucleusB()          const              {return !fInNucleusA;}
    Bool_t     IsProton()              const              {return (fType==1);}
    Bool_t     IsSpectator()           const              {return !fNColl;}
    Bool_t     IsWounded()             const              {return fNColl>0;}
    void       Reset()                                    {fNColl=0;}
    void       RotateXYZ(Double_t phi, Double_t theta);
    void       RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ);
    void       SetEnergy(Double_t en)                     {fEn = en;}
    void       SetInNucleusA()                            {fInNucleusA=1;}
    void       SetInNucleusB()                            {fInNucleusA=0;}
    void       SetNColl(Int_t nc)                         {fNColl = nc;}
    void       SetType(Bool_t b)                          {fType = b;}
    void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
   // ClassDef(TGlauNucleon,4) // TGlauNucleon class
};

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
  private:
    Int_t      fN;                   //Number of nucleons
    Int_t      fZ;                   //Number of protons
    Double_t   fR;                   //Parameters of 3pF function
    Double_t   fA;                   //Parameters of 3pF function 
    Double_t   fW;                   //Parameters of 3pF function
    Double_t   fR2;                  //Parameters of 3pF function (for p and n separately)
    Double_t   fA2;                  //Parameters of 3pF function (for p and n separately)
    Double_t   fW2;                  //Parameters of 3pF function (for p and n separately)
    Double_t   fBeta2;               //Beta2 (deformed nuclei) 
    Double_t   fBeta3;               //Beta3 (deformed nuclei) 
    Double_t   fBeta4;               //Beta4 (deformed nuclei) 
    Double_t   fGamma;               //Gamma (deformed nuclei) 
    Double_t   fMinDist;             //Minimum separation distance
    Double_t   fNodeDist;            //Average node distance (set to <=0 if you do not want the "crystal lattice")
    Double_t   fSmearing;            //Node smearing (relevant if fNodeDist>0)
    Int_t      fRecenter;            //=1 by default (0=no recentering, 1=recenter all, 2=recenter displacing only one nucleon, 3=recenter by rotating around z and shift in z, 4=recenter by rotation only, 5=recenter in transverse plane)
    Int_t      fLattice;             //=0 use HCP by default (1=PCS, 2=BCC, 3=FCC)
    Double_t   fSmax;                //Maximum magnitude of cms shift tolerated (99, ie all by default) 
    Int_t      fF;                   //Type of radial distribution
    Int_t      fTrials;              //Store trials needed to complete nucleus
    Int_t      fNonSmeared;          //Store number of non-smeared-node nucleons
    Double_t   fWeight;              //Weight of nucleus (for event-by-event weighting)
    TF1*       fFunc1;               //!Probability density function rho(r)
    TF1*       fFunc2;               //!Probability density function rho(r) -> if set 1 is for p, 2 is for n
    TF2*       fFunc3;               //!Probability density function rho(r,theta) for deformed nuclei
    TObjArray* fNucleons;            //!Array of nucleons
    Double_t   fPhiRot;              //!Angle phi for nucleus
    Double_t   fThetaRot;            //!Angle theta for nucleus
    Double_t   fXRot;                //!Angle around X axis for nucleus
    Double_t   fYRot;                //!Angle around Y axis for nucleus
    Double_t   fZRot;                //!Angle around Z axis for nucleus
    Float_t*** fNucArr;              //!Array of events, up to ~20 nucleons (only for small nuclei), 3 coordinates
    Int_t      fNucCounter;          //!Event counter
    TBits     *fIsUsed;              //!Bits for lattice use  
    Double_t   fMaxR;                //!Maximum radius (15fm)
    TGraph    *fInputDens;           //!Input density function (if provided via data points)
    void       AllocateNucleons();
    void       RandomizeNucleons();
    void       AllocateNucArr();
    void       ReadNucArr(const char* fname);
    void       FreeNucArr();
    void       SelectFromNucArr();
    void       Lookup(const char* name);
    Bool_t     TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const;

  public:
    TGlauNucleus(const char* iname="Pb", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
    virtual ~TGlauNucleus();
    Double_t   CalcMinDist()      const;
    Double_t   CalcRmsRadius()    const;
    void       Draw(Option_t* option="") { if (fFunc1) fFunc1->Draw(option); else if (fFunc2) fFunc2->Draw(option); else if (fFunc3) fFunc3->Draw(option); else TObject::Draw(option);}
    void       Draw(Double_t xs, Int_t colp, Int_t cols);
    Double_t   GetA()             const {return fA;}
    TF1*       GetFunc1()         const {return GetFuncP();}
    TF1*       GetFunc2()         const {return GetFuncN();}
    TF2*       GetFunc3()         const {return GetFuncDef();}
    TF1*       GetFuncP()         const {return fFunc1;}
    TF1*       GetFuncN()         const {return fFunc2;}
    TF2*       GetFuncDef()       const {return fFunc3;}
    TF1*       GetDens(Bool_t n=0)const;
    Double_t   GetSqrtMeanR2()    const;
    Double_t   GetMinDist()       const {return fMinDist;}
    Int_t      GetN()             const {return fN;}
    Double_t   GetNodeDist()      const {return fNodeDist;}
    TObjArray *GetNucleons()      const {return fNucleons;}
    Int_t      GetRecenter()      const {return fRecenter;}
    Double_t   GetR()             const {return fR;}
    Double_t   GetPhiRot()        const {return fPhiRot;}
    Double_t   GetThetaRot()      const {return fThetaRot;}
    Int_t      GetTrials()        const {return fTrials;}
    Int_t      GetNonSmeared()    const {return fNonSmeared;}
    Double_t   GetShiftMax()      const {return fSmax;}
    Double_t   GetW()             const {return fW;}
    Double_t   GetWeight()        const {return fWeight;}
    Double_t   GetXRot()          const {return fXRot;}
    Double_t   GetYRot()          const {return fYRot;}
    Int_t      GetZ()             const {return fZ;}
    Double_t   GetZRot()          const {return fZRot;}
    void       SetA(Double_t ia, Double_t ia2=-1);
    void       SetBeta(Double_t b2, Double_t b3, Double_t b4, Double_t g); 
    void       SetBeta(Double_t b2, Double_t b4); 
    void       SetGamma(Double_t g); 
    void       SetLattice(Int_t i)               {fLattice=i;}
    void       SetMinDist(Double_t min)          {fMinDist=min;}
    void       SetN(Int_t in)                    {fN=in;}
    void       SetNodeDist(Double_t nd)          {fNodeDist=nd;}
    void       SetR(Double_t ir, Double_t ir2=-1);
    void       SetRecenter(Int_t b)              {fRecenter=b;}
    void       SetShiftMax(Double_t s)           {fSmax=s;}
    void       SetSmearing(Double_t s)           {fSmearing=s;}
    void       SetW(Double_t iw);
    void       SetWeight(Double_t w)             {fWeight=w;}
    TVector3  &ThrowNucleons(Double_t xshift=0.);
   // ClassDef(TGlauNucleus,7) // TGlauNucleus class
};
//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
  public:
    class Event {
      public:
        Float_t Npart;       //Number of wounded (participating) nucleons in current event
        Float_t Ncoll;       //Number of binary collisions in current event
        Float_t Nhard;       //Number of hard collisions in current event (based on fHardFrac)
        Float_t Nmpi;        //Number of MPI 
        Float_t B;           //Impact parameter (b)
        Float_t BNN;         //Average NN impact parameter
        Float_t Ncollpp;     //Ncoll pp
        Float_t Ncollpn;     //Ncoll pn
        Float_t Ncollnn;     //Ncoll nn
        Float_t VarX;        //Variance of x of wounded nucleons
        Float_t VarY;        //Variance of y of wounded nucleons
        Float_t VarXY;       //Covariance of x and y of wounded nucleons
        Float_t NpartA;      //Number of wounded (participating) nucleons in Nucleus A
        Float_t NpartB;      //Number of wounded (participating) nucleons in Nucleus B
        Float_t Npart0;      //Number of singly-wounded (participating) nucleons
        Float_t NpartAn;     //Number of wounded (participating) neutrons in Nucleus A
        Float_t NpartBn;     //Number of wounded (participating) neutrons in Nucleus B
        Float_t Npart0n;     //Number of singly-wounded (participating) neutrons
        Float_t AreaW;       //Area defined by width of participants
        Float_t SpecA;       //Spectator neutrons in nucleus A
        Float_t SpecB;       //Spectator neutrons in nucleus B
        Float_t Weight;      //Weight of event (needed for e-by-e weighting)
        Float_t Psi1;        //Psi1
        Float_t Ecc1;        //Eps1
        Float_t Psi2;        //Psi2
        Float_t Ecc2;        //Eps2
        Float_t Psi3;        //Psi3
        Float_t Ecc3;        //Eps3
        Float_t Psi4;        //Psi4
        Float_t Ecc4;        //Eps4
        Float_t Psi5;        //Psi5
        Float_t Ecc5;        //Eps5
        Float_t AreaA;       //Area defined by "and" of participants
        Float_t AreaO;       //Area defined by "or" of participants
        Float_t X0;          //Production point in x
        Float_t Y0;          //Production point in y
        Float_t Phi0;        //Direction in phi
        Float_t Length;      //Length in phi0
        Float_t MeanX;       //<x> of wounded nucleons
        Float_t MeanY;       //<y> of wounded nucleons
        Float_t MeanX2;      //<x^2> of wounded nucleons
        Float_t MeanY2;      //<y^2> of wounded nucleons
        Float_t MeanXY;      //<xy> of wounded nucleons
        Float_t MeanXSystem; //<x> of all nucleons
        Float_t MeanYSystem; //<y> of all nucleons  
        Float_t MeanXA;      //<x> of nucleons in nucleus A
        Float_t MeanYA;      //<y> of nucleons in nucleus A
        Float_t MeanXB;      //<x> of nucleons in nucleus B
        Float_t MeanYB;      //<y> of nucleons in nucleus B
        Float_t PhiA;        //Phi angle nucleus A
        Float_t ThetaA;      //Theta angle nucleus B
        Float_t PhiB;        //Phi angle nucleus B
        Float_t ThetaB;      //Theta angle nucleus B
        void    Reset()      {Npart=0;Ncoll=0;Nhard=0;Nmpi=0;B=0;BNN=0;Ncollpp=0;Ncollpn=0;Ncollnn=0;VarX=0;VarY=0;VarXY=0;NpartA=0;NpartB=0;Npart0=0;NpartAn=0;NpartBn=0;Npart0n=0;AreaW=0;SpecA=0;SpecB=0;Weight=0;
                              Psi1=0;Ecc1=0;Psi2=0;Ecc2=0;Psi3=0;Ecc3=0;Psi4=0;Ecc4=0;Psi5=0;Ecc5=0;
                              AreaA=0;AreaO=0;X0=0;Y0=0;Phi0=0;Length=0;
                              MeanX=0;MeanY=0;MeanX2=0;MeanY2=0;MeanXY=0;MeanXSystem=0;MeanYSystem=0;MeanXA=0;MeanYA=0;MeanXB=0;MeanYB=0;
                              PhiA=0;ThetaA=0;PhiB=0;ThetaB=0;} // order must match that given in vars below
       // ClassDef(TGlauberMC::Event, 2)
    };

  protected:
    TGlauNucleus  fANucleus;       //Nucleus A
    TGlauNucleus  fBNucleus;       //Nucleus B
    Double_t      fXSect;          //Nucleon-nucleon cross section
    Double_t      fXSectNP;        //Proton-Neutron cross section (needed for Hades energies)
    Double_t      fXSectOmega;     //StdDev of Nucleon-nucleon cross section
    Double_t      fXSectLambda;    //Jacobian from tot to inelastic (Strikman)
    Double_t      fXSectEvent;     //Event value of Nucleon-nucleon cross section
    TObjArray*    fNucleonsA;      //Array of nucleons in nucleus A
    TObjArray*    fNucleonsB;      //Array of nucleons in nucleus B
    TObjArray*    fNucleons;       //Array which joins Nucleus A & B
    Int_t         fAN;             //Number of nucleons in nucleus A
    Int_t         fBN;             //Number of nucleons in nucleus B
    TNtuple*      fNt;             //Ntuple for results (created, but not deleted)
    Double_t      fEvents;         //Number of events with at least one collision
    Double_t      fTotalEvents;    //All events within selected impact parameter range
    Double_t      fBmin;           //Minimum impact parameter to be generated
    Double_t      fBmax;           //Maximum impact parameter to be generated
    Double_t      fHardFrac;       //Fraction of cross section used for Nhard (def=0.65)
    Int_t         fDetail;         //Detail to store (99=all by default)
    Bool_t        fCalcArea;       //If true calculate overlap area via grid (slow, off by default)
    Bool_t        fCalcLength;     //If true calculate path length (slow, off by default)
    Bool_t        fDoCore;         //If true calculate area and eccentricy only for core participants (off by default)
    Bool_t        fDoAAGG;         //If true do Glauber Gribov also for AA
    Double_t      fSigH;           //Sigma hard process
    Bool_t        fShadow;         //If true use shadowed cross section
    Double_t      fOmega;          //Omega parameter for NN profile (default=0)
    Int_t         fMaxNpartFound;  //Largest value of Npart obtained
    Double_t      fPsiN[10];       //Psi N
    Double_t      fEccN[10];       //Ecc N
    Double_t      f2Cx;            //Two-component x
    TF1          *fPTot;           //Cross section distribution
    TF1          *fNNProf;         //NN profile (hard-sphere == 0 by default)
    Event         fEv;             //Glauber event (results of calculation stored in tree)
    Bool_t        fBC[999][999];   //Array to record binary collision
    Bool_t        CalcResults(Double_t bgen);
    Bool_t        CalcEvent(Double_t bgen);
    static Double_t gDefOmega;     //default omega (-1)

  public:
    TGlauberMC(const char* NA = "Pb", const char* NB = "Pb", Double_t xsect = 42, Double_t xsectsigma=0, Double_t xsectnp=0);
    virtual ~TGlauberMC() {delete fNt; fNt=0; delete fNucleons; fNucleons=0; delete fPTot; fPTot=0; if (fOmega!=99) delete fNNProf; fNNProf=0;}

    Double_t            CalcDens(TF1 &prof, Double_t xval, Double_t yval) const;
    void                Draw(Option_t* option="");
    Double_t            GetB()                 const {return fEv.B;}
    Double_t            GetBNN()               const {return fEv.BNN;}
    Double_t            GetBmax()              const {return fBmax;}
    Double_t            GetBmin()              const {return fBmin;}
    Double_t            GetEcc(Int_t i=2)      const {return fEccN[i];}
    Double_t            GetHardFrac()          const {return fHardFrac;}
    Double_t            GetMeanX()             const {return fEv.MeanX;}
    Double_t            GetMeanXParts()        const {return fEv.MeanX;}
    Double_t            GetMeanXSystem()       const {return fEv.MeanXSystem;}
    Double_t            GetMeanY()             const {return fEv.MeanY;}
    Double_t            GetMeanYParts()        const {return fEv.MeanY;}
    Double_t            GetMeanYSystem()       const {return fEv.MeanYSystem;}
    Double_t            GetPsi(Int_t i=2)      const {return fPsiN[i];}
    Double_t            GetSx2()               const {return fEv.VarX;}    
    Double_t            GetSxy()               const {return fEv.VarXY;}    
    Double_t            GetSy2()               const {return fEv.VarY;}    
    Double_t            GetTotXSect()          const;
    Double_t            GetTotXSectErr()       const;
    Double_t            GetXSectEvent()        const {return fXSectEvent;}
    Int_t               GetNcoll()             const {return fEv.Ncoll;}
    Int_t               GetNcollnn()           const {return fEv.Ncollnn;}
    Int_t               GetNcollpn()           const {return fEv.Ncollpn;}
    Int_t               GetNcollpp()           const {return fEv.Ncollpp;}
    Int_t               GetNpart()             const {return fEv.Npart;}
    Int_t               GetNmpi()              const {return fEv.Nmpi;}
    Int_t               GetNpart0()            const {return fEv.Npart0;}
    Int_t               GetNpartA()            const {return fEv.NpartA;}
    Int_t               GetNpartB()            const {return fEv.NpartB;}
    Int_t               GetNpartFound()        const {return fMaxNpartFound;}
    Double_t            GetOmega()             const {return fOmega;}
    Double_t            GetSpecA()             const {return fEv.SpecA;}
    Double_t            GetSpecB()             const {return fEv.SpecB;}
    Double_t            GetWeight()            const {return fEv.Weight;}
    TF1*                GetXSectDist()         const {return fPTot;}
    TGlauNucleus*       GetNucleusA()                {return &fANucleus;}
    TGlauNucleus*       GetNucleusB()                {return &fBNucleus;}
    TNtuple*            GetNtuple()            const {return fNt;}
    TObjArray          *GetNucleons();
    const Event        &GetEvent()             const {return fEv;}
    const Event        *GetEvent()                   {return &fEv;}
    const TGlauNucleus* GetNucleusA()          const {return &fANucleus;}
    const TGlauNucleus* GetNucleusB()          const {return &fBNucleus;}
    Bool_t              IsBC(Int_t i, Int_t j) const {return fBC[i][j];}
    Bool_t              NextEvent(Double_t bgen=-1);
    void                Reset()                      {delete fNt; fNt=0; }
    Bool_t              ReadNextEvent(Bool_t calc=1, const char *fname=0);       
    void                Run(Int_t nevents,Double_t b=-1);
    void                Set2Cx(Double_t x)           {f2Cx = x;}
    void                SetBmax(Double_t bmax)       {fBmax = bmax;}
    void                SetBmin(Double_t bmin)       {fBmin = bmin;}
    void                SetCalcAAGG(Bool_t b)        {fDoAAGG = b;}
    void                SetCalcArea(Bool_t b)        {fCalcArea = b;}
    void                SetCalcCore(Bool_t b)        {fDoCore = b;}
    void                SetCalcLength(Bool_t b)      {fCalcLength = b;}
    void                SetDetail(Int_t d)           {fDetail = d;}
    void                SetHardFrac(Double_t f)      {fHardFrac=f;}
    void                SetLattice(Int_t i)          {fANucleus.SetLattice(i); fBNucleus.SetLattice(i);}
    void                SetMinDistance(Double_t d)   {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
    void                SetNNProf(TF1 *f1)           {fNNProf = f1; fOmega=99;}
    void                SetOmega(Double_t omega=-1);
    void                SetNodeDistance(Double_t d)  {fANucleus.SetNodeDist(d); fBNucleus.SetNodeDist(d);}
    void                SetRecenter(Int_t b)         {fANucleus.SetRecenter(b); fBNucleus.SetRecenter(b);}
    void                SetShiftMax(Double_t s)      {fANucleus.SetShiftMax(s); fBNucleus.SetShiftMax(s);}
    void                SetShadowing(Bool_t b)       {fShadow=b;}
    void                SetSigmaHard(Double_t s)     {fSigH=s;}
    void                SetSmearing(Double_t s)      {fANucleus.SetSmearing(s); fBNucleus.SetSmearing(s);}
    void                SetXSect(Double_t sig)       {fXSect=sig;}
    void                SetXSectNP(Double_t sig)     {fXSectNP=sig;}
    void                SetXSectDist(TF1 *f1)        {fPTot = f1;}
    const char         *Str()                  const {return Form("TGlauberMC%s_%s%s-snn%.1f-md%.1f-om%.1f-nd%.1f-rc%d-smax%.1f",Version(),fANucleus.GetName(),fBNucleus.GetName(),fXSect,fBNucleus.GetMinDist(),fOmega,fBNucleus.GetNodeDist(),fBNucleus.GetRecenter(),fBNucleus.GetShiftMax());}
    static void         PrintVersion()               {cout << "TGlauberMC " << Version() << endl;}
    static Double_t     GetDefOmega()                {return gDefOmega;}
    static void         SetDefOmega(Double_t om)     {gDefOmega=om;}
    static const char  *Version()                    {return "v3.3.1";}
   // ClassDef(TGlauberMC,7) // TGlauberMC class
};
#endif

#endif // T_GLAUBER_MC_H
