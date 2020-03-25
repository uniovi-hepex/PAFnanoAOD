#ifndef MT2_BISECT_H
#define MT2_BISECT_H 1
#include "TLorentzVector.h"
#include "Lepton.h"
#include "Jet.h"


/*The user can change the desired precision below, the larger one of the following two definitions is used. Relative precision less than 0.00001 is not guaranteed to be achievable--use with caution*/ 

#define RELATIVE_PRECISION 0.00001 //defined as precision = RELATIVE_PRECISION * scale, where scale = max{Ea, Eb}
#define ABSOLUTE_PRECISION 0.0     //absolute precision for mt2, unused by default


//Reserved for expert
#define MIN_MASS  0.1   //if ma<MINMASS and mb<MINMASS, use massless code
#define ZERO_MASS 0.000 //give massless particles a small mass
#define SCANSTEP 0.1
class mt2
{  
  public:

    mt2();
    void   mt2_bisect();
    void   mt2_massless();
    void   set_momenta(double *pa0, double *pb0, double* pmiss0);
    void   set_mn(double mn);
    double get_mt2();
    void   print();
    int    nevt;
  private:  

    bool   solved;
    bool   momenta_set;
    double mt2_b;

    int    nsols(double Dsq);
    int    nsols_massless(double Dsq);
    inline int    signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5);
    inline int    signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5);
    int scan_high(double &Deltasq_high);
    int find_high(double &Deltasq_high);
    //data members
    double pax, pay, ma, Ea;
    double pmissx, pmissy;
    double pbx, pby, mb, Eb;
    double mn, mn_unscale;

    //auxiliary definitions
    double masq, Easq;
    double mbsq, Ebsq;
    double pmissxsq, pmissysq;
    double mnsq;

    //auxiliary coefficients
    double a1, b1, c1, a2, b2, c2, d1, e1, f1, d2, e2, f2;
    double d11, e11, f12, f10, d21, d20, e21, e20, f22, f21, f20;

    double scale;
    double precision;

    //ClassDef(mt2,0);
};

int mt2::signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5)
{
  int nsc;
  nsc=0;
  if(t1*t2>0) nsc++;
  if(t2*t3>0) nsc++;
  if(t3*t4>0) nsc++;
  if(t4*t5>0) nsc++;
  return nsc;
}

int mt2::signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5)
{
  int nsc;
  nsc=0;
  if(t1*t2<0) nsc++;
  if(t2*t3<0) nsc++;
  if(t3*t4<0) nsc++;
  if(t4*t5<0) nsc++;
  return nsc;
}

Float_t getMT2(TLorentzVector plep1, TLorentzVector plep2, TLorentzVector pmet, Float_t mass);
Float_t getMT2ll(Lepton l1, Lepton l2, Float_t met, Float_t met_phi);
Float_t getMT2llLepScale(Lepton l1, Lepton l2, Float_t met, Float_t met_phi, int id, int dir);
Float_t getMT2lb(TLorentzVector plep0, TLorentzVector plep1, TLorentzVector pjet0, TLorentzVector pjet1, Float_t met, Float_t met_phi);
Float_t getMT2lb(TLorentzVector plep0, TLorentzVector plep1, TLorentzVector pmet, std::vector<Jet> jets);
Float_t getMT2lb(TLorentzVector plep0, TLorentzVector plep1, TLorentzVector pjet0, TLorentzVector pjet1, TLorentzVector pmet);

#endif
