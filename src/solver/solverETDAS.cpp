
#include "spida/solver/solverETDAS.h"
#include <cmath>
#include <thread>

namespace spida{

SolverAS_ETD::SolverAS_ETD(ModelCV* model,double sf,double qv)
 :  SolverCV_AS(model,sf,qv)
{
    m_mode_cutoff = 0.01;
    m_contour_radi = 1.0;
    m_contourM = 32;
}

ETD34::ETD34(ModelCV* model)
  : SolverAS_ETD(model,0.84,4),m_sz(SolverCV::size()),
    L(model->linOp()),EL(m_sz), EL2(m_sz), 
    N1(m_sz), N2(m_sz), N3(m_sz), N4(m_sz), N5(m_sz), 
    tempK(m_sz), a21(m_sz), a31(m_sz), a32(m_sz),
    a41(m_sz), a43(m_sz), a51(m_sz), a52(m_sz), a54(m_sz), r1(m_sz),
    thmgt(model->threadManager())
{
    c1 = 0.0; c2 = 1.0/2; c3 = 1.0/2; c4 = 1.0; c5 = 1.0;
    N1_init = false;
    statCenter().setHeader("ETD34 STATS");
    statCenter().addCounter("Nonlinear Function Evaluations",4);
}

void ETD34::worker_coeff(double dstep,int tid){

    double dt = dstep;
    std::vector<dcmplx> r(contourPoints(),0.0);
    for(int j = 0; j < contourPoints(); j++){
        double expv = 2.0*PI*(j+0.5)/contourPoints();
        r[j] = contourRadius()*exp(ii*expv);
    }
    double mult = dt/contourPoints();
    for(int i = tid; i < m_sz; i+=thmgt.getNumThreads()){
        dcmplx Lval(dt*L[i]);
        EL[i] = exp(Lval);
        EL2[i] = exp(c2*Lval);
        if(abs(Lval) < modeCutoff()){
            dcmplx a21val(0.0,0.0);
            dcmplx a31val(0.0,0.0); dcmplx a32val(0.0,0.0);
            dcmplx a41val(0.0,0.0); dcmplx a43val(0.0,0.0);
            dcmplx a51val(0.0,0.0); dcmplx a52val(0.0,0.0); dcmplx a54val(0.0,0.0);
            for(int j = 0; j < contourPoints(); j++){
                dcmplx z(Lval + r[j]);
                a21val  = a21val + (-1.0 + exp(z/2.0))/z;
                a31val = a31val + (4.0 + z + exp(z/2.0)*(-4.0 + z))/(pow(z,2));
                a32val = a32val - 2.0*(2.0 + z - 2.0*exp(z/2.0))/(pow(z,2));
                a41val = a41val + (2.0 + z + exp(z)*(-2.0 + z))/pow(z,2);
                a43val = a43val + 2.0*(-1.0 + exp(z) - z)/pow(z,2);
                a51val = a51val + (-4.0 - z + exp(z)*(4.0 - 3.0*z+pow(z,2)))/pow(z,3);
                a52val = a52val + 2.0*(2.0 + z + exp(z)*(-2.0+z))/pow(z,3);
                a54val = a54val + (-4.0-3.0*z-pow(z,2) + exp(z)*(4.0-z))/pow(z,3);
            }
            a21[i] = mult*a21val; a31[i] = mult*a31val; a32[i] = mult*a32val;
            a41[i] = mult*a41val; a43[i] = mult*a43val; a51[i] = mult*a51val;
            a52[i] = mult*a52val; a54[i] = mult*a54val;
        }
        else{
            a21[i]  = dt*(-1.0 + exp(Lval/2.0))/Lval;
            a31[i] = dt*(4.0 + Lval + exp(Lval/2.0)*(-4.0 + Lval))/(pow(Lval,2));
            a32[i] =  -2.0*dt*(2.0 + Lval - 2.0*exp(Lval/2.0))/(pow(Lval,2));
            a41[i] = dt*(2.0 + Lval + exp(Lval)*(-2.0 + Lval))/pow(Lval,2);
            a43[i] = 2.0*dt*(-1.0 + exp(Lval) - Lval)/pow(Lval,2);
            a51[i] = dt*(-4.0 - Lval + exp(Lval)*(4.0 - 3.0*Lval+pow(Lval,2)))/pow(Lval,3);
            a52[i] = 2.0*dt*(2.0 + Lval + exp(Lval)*(-2.0+Lval))/pow(Lval,3);
            a54[i] = dt*(-4.0-3.0*Lval-pow(Lval,2) + exp(Lval)*(4.0-Lval))/pow(Lval,3);
        }
    }
}

void ETD34::updateCoefficients(double dt)
{
    std::vector<std::thread> threads;
    for(unsigned int i = 1; i < thmgt.getNumThreads(); i++)
        threads.push_back(std::thread(&ETD34::worker_coeff,this,dt,i));
    worker_coeff(dt,0);
    for(auto& thread : threads)
        thread.join();
}


void ETD34::updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& ynew,\
        std::vector<dcmplx>& errVec)
{
    if(!N1_init){
        SolverCV::model().nonLinResponse(in,N1);
        N1_init = true;
    }
    else if(SolverCV_AS::accept()){
        for(int i = 0; i < m_sz; i++) 
            N1[i] = N5[i]; 
    }

    for(int i = 0; i < m_sz; i++)
        tempK[i] =  EL2[i]*in[i] + a21[i]*N1[i];

    SolverCV::model().nonLinResponse(tempK,N2);
  
    for(int i = 0; i < m_sz; i++) 
        tempK[i] = EL2[i]*in[i] + a31[i]*N1[i]+a32[i]*N2[i];

    SolverCV::model().nonLinResponse(tempK,N3);
  
    for(int i = 0; i < m_sz; i++) 
        tempK[i] = EL[i]*in[i] + a41[i]*N1[i] + a43[i]*N3[i];

    SolverCV::model().nonLinResponse(tempK,N4);
  
    for(int i = 0; i < m_sz; i++) 
        ynew[i] = EL[i]*in[i] + a51[i]*N1[i] + a52[i]*(N2[i]+N3[i]) + a54[i]*N4[i];
  
    SolverCV::model().nonLinResponse(ynew,N5);

    for (int i = 0; i < m_sz; i++)
        errVec[i] = a54[i]*(N4[i] - N5[i]);
}

ETD35::ETD35(ModelCV* model)
  : SolverAS_ETD(model,0.84,4),m_sz(SolverCV::size()),
    L(model->linOp()),EL(m_sz),EL2(m_sz),EL4(m_sz),EL5(m_sz),
    N1(m_sz),N2(m_sz),N3(m_sz),N4(m_sz),N5(m_sz),N6(m_sz),N7(m_sz),
    tempK(m_sz),a21(m_sz),a31(m_sz),a32(m_sz),a41(m_sz),a43(m_sz),
    a51(m_sz),a52(m_sz),a54(m_sz),a61(m_sz),a62(m_sz),a63(m_sz),a65(m_sz),
    a71(m_sz),a73(m_sz),a74(m_sz),a75(m_sz),a76(m_sz),
    thmgt(model->threadManager())
{
    c1 = 0.0; c2 = 1.0/4.0; c3 = 1.0/4.0; c4 = 1.0/2.0;  
    c5 = 3.0/4.0; c6 = 1.0; c7 = 1.0; 
    N1_init = false;
    statCenter().setHeader("ETD35 STATS");
    statCenter().addCounter("Nonlinear Function Evaluations",6);
}

void ETD35::worker_coeff(double dstep,int tid){

  double ds = dstep;
  std::vector<dcmplx> r(contourPoints(),0.0);
  for(int j = 0; j < contourPoints(); j++){
    double expv = 2.0*PI*(j+0.5)/contourPoints();
    r[j] = contourRadius()*exp(ii*expv);
  }

  for(int i = tid; i < m_sz; i+=thmgt.getNumThreads())
  {
    dcmplx Lval(ds*L[i]);
    EL[i] = exp(Lval);
    EL2[i] = exp(c2*Lval);
    EL4[i] = exp(c4*Lval);
    EL5[i] = exp(c5*Lval);
    if(abs(Lval) < modeCutoff()){
      dcmplx a21val(0.0,0.0);
      dcmplx a31val(0.0,0.0); dcmplx a32val(0.0,0.0);
      dcmplx a41val(0.0,0.0); dcmplx a43val(0.0,0.0);
      dcmplx a51val(0.0,0.0); dcmplx a52val(0.0,0.0); dcmplx a54val(0.0,0.0);
      dcmplx a61val(0.0,0.0); dcmplx a62val(0.0,0.0); dcmplx a63val(0.0,0.0); dcmplx a65val(0.0,0.0);
      dcmplx a71val(0.0,0.0); dcmplx a73val(0.0,0.0); dcmplx a74val(0.0,0.0);
      dcmplx a75val(0.0,0.0); dcmplx a76val(0.0,0.0);
      for(int j = 0; j < contourPoints(); j++){
        dcmplx z(Lval + r[j]);
        a21val+= (-1.0 + exp((z)/4.0))/(z);
        a31val+= (4.0 + exp((z)/4.0)*(-4.0 + z))/(pow(z,2));
        a32val+= -((4.0 - 4.0*exp((z)/4.0) + z)/(pow(z,2)));
        a41val+= (4.0 + z + exp((z)/2.0)*(-4.0 + z))/(pow(z,2));
        a43val+= (-2.0*(2.0 - 2.0*exp((z)/2.0) + z))/(pow(z,2));
        a51val+= (4.0 + z + 2.0*exp((3.0*z)/4.0)*(-2.0 + z))/(2.0*pow(z,2));
        a52val+= -(-1.0 + exp((3.0*z)/4.0))/(2.0*z);
        a54val+= -(4.0 - 4.0*exp((3.0*z)/4.0) + 3.0*z)/(2.0*pow(z,2));
        a61val+= (-118.0 - 41.0*z + exp(z)*(118.0 - 77.0*z))/(42.0*pow(z,2));
        a62val+= (8.0*(-1.0 + exp(z)))/(7.0*z);
        a63val+= (3.0*(58.0 + 21.0*z + exp(z)*(-58.0 + 37.0*z)))/(28.0*pow(z,2));
        a65val+= (-286.0 - 239.0*z + exp(z)*(286.0 - 47.0*z))/(84.0*pow(z,2));
        a71val+= (7.0*(-1620.0 - 626.0*z - 73.0*pow(z,2) + exp(z)*(1620.0 - 994.0*z + 257.0*pow(z,2))))/
          (2700.0*pow(z,3));
        a73val+= (900.0 + 1834.0*z + 287.0*pow(z,2) + exp(z)*(-900.0 - 934.0*z + 1097.0*pow(z,2)))/
          (1350.0*pow(z,3));
        a74val+= (2.0*(810.0 + 412.0*z + 56.0*pow(z,2) + exp(z)*(-810.0 + 398.0*z - 49.0*pow(z,2))))/
          (225.0*pow(z,3));
        a75val+= (540.0 - 1226.0*z - 1183.0*pow(z,2) + exp(z)*(-540.0 + 1766.0*z - 313.0*pow(z,2)))/
          (1350.0*pow(z,3));
        a76val+= (-10980.0 - 6722.0*z - 1741.0*pow(z,2) + exp(z)*(10980.0 - 4258.0*z + 509.0*pow(z,2)))/
          (2700.0*pow(z,3));
      }
      double mult = ds/contourPoints();
      a21[i] = mult*a21val; a31[i] = mult*a31val; a32[i] = mult*a32val; 
      a41[i] = mult*a41val; a43[i] = mult*a43val; 
      a51[i] = mult*a51val; a52[i] = mult*a52val; a54[i] = mult*a54val;
      a61[i] = mult*a61val; a62[i] = mult*a62val;
      a63[i] = mult*a63val; a65[i] = mult*a65val;
      a71[i] = mult*a71val; a73[i] = mult*a73val; a74[i] = mult*a74val; 
      a75[i] = mult*a75val; a76[i] = mult*a76val;
    }
    else{
      a21[i] = ds*(-1.0 + exp((Lval)/4.0))/(Lval);
      a31[i] = ds*(4.0 + exp((Lval)/4.0)*(-4.0 + Lval))/(pow(Lval,2));
      a32[i] = ds*-((4.0 - 4.0*exp((Lval)/4.0) + Lval)/(pow(Lval,2)));
      a41[i] = ds*(4.0 + Lval + exp((Lval)/2.0)*(-4.0 + Lval))/(pow(Lval,2));
      a43[i] = ds*(-2.0*(2.0 - 2.0*exp((Lval)/2.0) + Lval))/(pow(Lval,2));
      a51[i] = ds*(4.0 + Lval + 2.0*exp((3.0*Lval)/4.0)*(-2.0 + Lval))/(2.0*pow(Lval,2));
      a52[i] = ds*-(-1.0 + exp((3.0*Lval)/4.0))/(2.0*Lval);
      a54[i] = ds*-(4.0 - 4.0*exp((3.0*Lval)/4.0) + 3.0*Lval)/(2.0*pow(Lval,2));
      a61[i] = ds*(-118.0 - 41.0*Lval + exp(Lval)*(118.0 - 77.0*Lval))/(42.0*pow(Lval,2));
      a62[i] = ds*(8.0*(-1.0 + exp(Lval)))/(7.0*Lval);
      a63[i] = ds*(3.0*(58.0 + 21.0*Lval + exp(Lval)*(-58.0 + 37.0*Lval)))/(28.0*pow(Lval,2));
      a65[i] = ds*(-286.0 - 239.0*Lval + exp(Lval)*(286.0 - 47.0*Lval))/(84.0*pow(Lval,2));
      a71[i] = ds*(7.0*(-1620.0 - 626.0*Lval - 73.0*pow(Lval,2) + exp(Lval)*(1620.0 - 994.0*Lval + 257.0*pow(Lval,2))))/
        (2700.0*pow(Lval,3));
      a73[i] = ds*(900.0 + 1834.0*Lval + 287.0*pow(Lval,2) + exp(Lval)*(-900.0 - 934.0*Lval + 1097.0*pow(Lval,2)))/
        (1350.0*pow(Lval,3));
      a74[i] = ds*(2.0*(810.0 + 412.0*Lval + 56.0*pow(Lval,2) + exp(Lval)*(-810.0 + 398.0*Lval - 49.0*pow(Lval,2))))/
        (225.0*pow(Lval,3));
      a75[i] = ds*(540.0 - 1226.0*Lval - 1183.0*pow(Lval,2) + exp(Lval)*(-540.0 + 1766.0*Lval - 313.0*pow(Lval,2)))/
        (1350.0*pow(Lval,3));
      a76[i] = ds*(-10980.0 - 6722.0*Lval - 1741.0*pow(Lval,2) + exp(Lval)*(10980.0 - 4258.0*Lval + 509.0*pow(Lval,2)))/
        (2700.0*pow(Lval,3));
    }
  }
}

void ETD35::updateCoefficients(double dt)
{
  std::vector<std::thread> threads;
  for(unsigned int i = 1; i < thmgt.getNumThreads(); i++)
      threads.push_back(std::thread(&ETD35::worker_coeff,this,dt,i));
  worker_coeff(dt,0);
  for(auto& thread : threads)
      thread.join();
}

void ETD35::worker_stage2(const std::vector<dcmplx>& in,int sti,int endi)
{
    for(int i = sti; i < endi; i++)
        tempK[i] = EL2[i]*in[i] + a21[i]*N1[i];
}

void ETD35::worker_stage3(const std::vector<dcmplx>& in,int sti,int endi)
{
    for(int i = sti; i < endi; i++)
        tempK[i] = EL2[i]*in[i] + a31[i]*N1[i] + a32[i]*N2[i];
}

void ETD35::worker_stage4(const std::vector<dcmplx>& in,int sti,int endi)
{
    for(int i = sti; i < endi; i++)
        tempK[i] = EL4[i]*in[i] + a41[i]*N1[i] + a43[i]*N3[i];
}

void ETD35::worker_stage5(const std::vector<dcmplx>& in,int sti,int endi)
{
    for(int i = sti; i < endi; i++)
        tempK[i] = EL5[i]*in[i] + a51[i]*N1[i] + a52[i]*(N2[i]-N3[i]) + a54[i]*N4[i];
}
void ETD35::worker_stage6(const std::vector<dcmplx>& in,int sti,int endi)
{
    for(int i = sti; i < endi; i++){
        tempK[i] = EL[i]*in[i] + a61[i]*N1[i] + a62[i]*(N2[i]-(3.0/2)*N4[i]) \
            + a63[i]*N3[i] + a65[i]*N5[i];
    }
}
void ETD35::worker_stage7(const std::vector<dcmplx>& in,std::vector<dcmplx>& ynew,\
        std::vector<dcmplx>& errVec,int sti,int endi)
{
    for(int i = sti; i < endi; i++){
        ynew[i] = EL[i]*in[i] + a71[i]*N1[i] + a73[i]*N3[i] + a74[i]*N4[i] \
            + a75[i]*N5[i] + a76[i]*N6[i];
        errVec[i] = a75[i]*(-N1[i] + 4.0*N3[i] - 6.0*N4[i] + 4.0*N5[i] - N6[i])/4.0;
    }
}


void ETD35::updateStages(const std::vector<dcmplx>& in,\
        std::vector<dcmplx>& ynew,std::vector<dcmplx>& errVec)
{
    if(!N1_init){
        SolverCV::model().nonLinResponse(in,N1);
        N1_init = true;
    }
    if(SolverCV_AS::accept())
        SolverCV::model().nonLinResponse(ynew,N1);

    unsigned int nthreads = thmgt.getNumThreads();
    std::vector<unsigned int> bounds = thmgt.getBounds(m_sz);
    std::vector<std::thread*> threads;

    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage2,this,std::ref(in),bounds[i],bounds[i+1]));
    worker_stage2(in,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();

    SolverCV::model().nonLinResponse(tempK,N2);

    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage3,this,std::ref(in),bounds[i],bounds[i+1]));
    worker_stage3(in,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();

    SolverCV::model().nonLinResponse(tempK,N3);

    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage4,this,std::ref(in),bounds[i],bounds[i+1]));
    worker_stage4(in,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();

    SolverCV::model().nonLinResponse(tempK,N4);
    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage5,this,std::ref(in),bounds[i],bounds[i+1]));
    worker_stage5(in,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();

    SolverCV::model().nonLinResponse(tempK,N5);

    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage6,this,std::ref(in),bounds[i],bounds[i+1]));
    worker_stage6(in,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();

    SolverCV::model().nonLinResponse(tempK,N6);

    for(unsigned int i = 0; i < nthreads-1; i++)
        threads.push_back(new std::thread(&ETD35::worker_stage7,this,std::ref(in),std::ref(ynew),std::ref(errVec),bounds[i],bounds[i+1]));
    worker_stage7(in,ynew,errVec,bounds[nthreads-1],bounds[nthreads]);
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();
}



}


