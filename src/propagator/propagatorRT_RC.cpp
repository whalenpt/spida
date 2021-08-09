
/*
#include "shapecreate.h"
#include "propagators.h"
#include "fileAux.hpp"
#include "report.h"
#include "constants.h"
#include "scales.h"
#include <algorithm>

PropagatorsRT_RC::PropagatorsRT_RC(spida::SpidaRT* ida,pw::ParamBin& bin) 
  : Propagators(ida,bin), spi(ida)
{
  U = new double[spi->getNT()*spi->getNr()];
  V = new double[2*spi->getNST()*spi->getNsr()];
  I = new double[spi->getNT()*spi->getNr()];
  EU = new double[2*spi->getNT()*spi->getNr()];

  double omegRef = bin.getDbl("RefFreq");
  double k0 = bin.getDbl("n0")*omegRef/uppe::LIGHT_SPEED;
  vars.add("k0",k0);
  vars.add("n0",bin.getDbl("n0"));
  vars.add("omeg0",omegRef);
  vars.add("vg",bin.getDbl("vg"));
  vars.add("oldz",0.0);
  vars.addBool("REPORT/FluenceR",bin.getBoolF("REPORT/FluenceR"));
  vars.addBool("REPORT/BeamWidthF",bin.getBoolF("REPORT/BeamWidthF"));
  arrayR_R = new double[spi->getNr()];
  using uppe::LIGHT_SPEED;
  using uppe::EPS0;
  bin.add("INPUT/E0",sqrt(2.0*bin.getDbl("INPUT/Intensity")/(LIGHT_SPEED*EPS0*vars.getDbl("n0"))));

  InitReport(bin);
  InitPulse(bin);
}

PropagatorsRT_RC::~PropagatorsRT_RC()
{
  delete [] U;
  delete [] V;
  delete [] EU;
}

void PropagatorsRT_RC::InitPulse(pw::ParamBin& bin)
{
  const double* r = spi->getR();
  const double* t = spi->getT();
  ShapeInit2D shapeBin(bin);
  shapeBin.getRealShape(r,t,U,spi->getNr(),spi->getNT());
  spi->RT_To_SRST(U,V);
  updateFields(0.0);
}

void PropagatorsRT_RC::updateFields(double z)
{
  using uppe::dcmplx;
  using uppe::ii;
  using uppe::LIGHT_SPEED;
  using uppe::EPS0;
  int nt = spi->getNT();
  int nr = spi->getNr();
  spi->SRST_To_RT(V,U);
  spi->SRST_To_RTc(V,EU);
  double k0 = vars.getDbl("k0");
  double omg = vars.getDbl("omeg0");
  const double* t = spi->getT();
  for(int i = 0; i < nr; i++) { 
    for(int j = 0; j < nt; j++) { 
      dcmplx temp(exp(ii*k0*z - ii*omg*t[j])*dcmplx(EU[2*i*nt + 2*j],EU[2*i*nt + 2*j+1]));
      EU[2*i*nt + 2*j] = temp.real();
      EU[2*i*nt + 2*j+1] = temp.imag();
      I[nt*i+j] = abs(pow(temp,2));
    }
  }

  if(vars.getBool("REPORT/FluenceR") || vars.getBool("REPORT/BeamWidthF"))
    spi->IT(I,arrayR_R);

  // Change units to V/m
  for(int i = 0; i < nr*nt; i++)
    U[i] = sqrt(2.0/(vars.getDbl("n0")*LIGHT_SPEED*EPS0))*U[i];
}

void PropagatorsRT_RC::InitReport(pw::ParamBin& bin)
{
  int nt = spi->getNT();
  int nr = spi->getNr();
  int nkr = spi->getNsr();
  int nst = spi->getNST();
  pw::ParamBin scales(pwScales::ScalesRT(bin));

  //fields::setSQ_SRST(V,rp,spi,bin,scales);
  if(bin.getBool("Report/SQ_SRST")){
    if(!scales.inBin("krmax"))
      scales.addDbl("krmax",spi->getSR()[nkr-1]);

    ReportDef2D* def = new ReportDef2D("SQ_SRST",spi->getSR(),spi->getST(),V,nkr,nst);
    fields::setSQ_SRSTp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setSQ_RT(U,rp,spi,bin,scales);
  if(bin.getBool("Report/SQ_RT")){
    ReportDef2D_R* def = new ReportDef2D_R("SQ_RT",spi->getR(),spi->getT(),I,nr,nt);
    fields::setSQ_RTp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setRT(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/RT")){
    ReportDef2D_R* def = new ReportDef2D_R("RT",spi->getR(),spi->getT(),U,nr,nt);
    fields::setRTp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setSRST(V,rp,spi,bin,scales);
  if(bin.getBool("Report/SRST")){
    if(!scales.inBin("krmax"))
      scales.addDbl("krmax",spi->getSR()[nkr-1]);
    ReportDef2D* def = new ReportDef2D("SRST",spi->getSR(),spi->getST(),V,nkr,nst);
    fields::setSRSTp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setSQ_T(U,rp,spi,bin,scales);
  if(bin.getBool("Report/SQ_T")){
    ReportDef1D_R* def = new ReportDef1D_R("SQ_T",spi->getT(),I,spi->getNr(),spi->getNT());
    def->setOnAxisD2(0);
    fields::setSQ_Tp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setSQ_ST(V,rp,spi,bin,scales);
  if(bin.getBool("Report/SQ_ST")){
    ReportDef1D* def = new ReportDef1D("SQ_ST",spi->getST(),V,spi->getNsr(),spi->getNST());
    def->setOnAxisD2(0);
    fields::setSQ_STp(def,bin,scales);
    rp->addReport(def);
  }

  if(bin.getBool("REPORT/T")){
    ReportDef1D_R* def = new ReportDef1D_R("T",spi->getT(),U,spi->getNr(),spi->getNT());
    def->setOnAxisD2(0);
    fields::setTp(def,bin,scales);
    rp->addReport(def);
  }

  if(bin.inBin("REPORT/T_SnapShot")){
    const double* r = spi->getR();
    int pics = bin.getInt("REPORT/T_SnapShot");
    double w0 = bin.getDbl("INPUT/R/BeamDiameter");
    int maxIndx;
    if(3*w0 < spi->getMaxR())
      maxIndx = spi->RadiToIndx(3*w0);
    else
      maxIndx = nr;
    if(pics > maxIndx)
      maxIndx = pics;
    int count = 0;
    for(int i =0;i<maxIndx;i+=maxIndx/pics){
      std::ostringstream stm;
      stm << count;
      std::string nm = "T_" + stm.str();
      ReportDef1D_R* def = new ReportDef1D_R(nm,spi->getT(),U,nr,nt);
      def->setOnAxisD2(i);
      fields::setTp(def,bin,scales);
      def->addFigsetting("rval",r[i]);
      rp->addReport(def);
      count++;
    }
  }


  //fields::setSQ_R(U,rp,spi,bin,scales);
  if(bin.getBool("Report/SQ_R")){
    ReportDef1D_R* def = new ReportDef1D_R("SQ_R",spi->getR(),I,spi->getNr(),spi->getNT());
    def->setOnAxisD1(spi->TimeToIndx(0.0));
    fields::setSQ_Rp(def,bin,scales);
    rp->addReport(def);
  }

  if(bin.inBin("REPORT/SQ_T_SnapShot")){
    const double* r = spi->getR();
    int pics = bin.getInt("REPORT/SQ_T_SnapShot");
    double w0 = bin.getDbl("INPUT/R/BeamDiameter");
    int maxIndx;
    if(3*w0 < spi->getMaxR())
      maxIndx = spi->RadiToIndx(3*w0);
    else
      maxIndx = nr;
    if(pics > maxIndx)
      maxIndx = pics;
    int count = 0;
    for(int i =0;i<maxIndx;i+=maxIndx/pics){
      std::ostringstream stm;
      stm << i;
      std::string nm = "SQ_T_" + stm.str();
      ReportDef1D_R* def = new ReportDef1D_R(nm,spi->getT(),I,nr,nt);
      def->setOnAxisD2(i);
      fields::setSQ_Tp(def,bin,scales);
      def->addFigsetting("rval",r[i]);
      rp->addReport(def);
      count++;
    }
  }

  //fields::setR(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/R")){
    ReportDef1D* def = new ReportDef1D("R",spi->getR(),U,nr,nt);
    def->setOnAxisD1(spi->TimeToIndx(0.0));
    rp->addReport(def);
  }

  //fields::setOnAxisInte(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/OnAxisInte")){
    ReportMax_R* def = new ReportMax_R("OnAxisIntensity",I,spi->getNr(),spi->getNT());
    def->setOnAxisD2(0);
    fields::setOnAxisIp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setPeakInte(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/PeakInte")){
    ReportMax_R* def = new ReportMax_R("PeakIntensity",I,spi->getNr(),spi->getNT());
    fields::setOnAxisIp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setMinBeamWidth(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/BeamWidth")){
    ReportWidth* def = new ReportWidth("BeamWidth",spi->getR(),I,spi->getNr(),spi->getNT());
    def->setOverAxisD1();
    fields::setMinBWp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setMinBeamWidth(U,rp,spi,bin,scales);
  if(bin.getBoolF("REPORT/BeamWidthF")){
    ReportWidth* def = new ReportWidth("BeamWidthF",spi->getR(),arrayR_R,spi->getNr());
    def->setOverAxisD1();
    fields::setMinBWp(def,bin,scales);
    rp->addReport(def);
  }

  //fields::setP_T(U,rp,spi,bin,scales);
  if(bin.getBool("REPORT/PhaseT")){
    ReportDef1D* def = new ReportDef1D("PhaseT",spi->getT(),EU,spi->getNT());
    def->setOnAxisD1(0);
    fields::setP_Tp(def,bin,scales);
    rp->addReport(def);
  }

  if(bin.getBoolF("REPORT/FluenceR")){
    std::string nm = "FluenceR";
    ReportDef1D_R* def = new ReportDef1D_R(nm,spi->getR(),arrayR_R,nr);
    fields::setFluence_Rp(def,bin,scales);
    rp->addReport(def);
  }

  fields::setDT_R(U,rp,spi,bin,scales);
}
*/


