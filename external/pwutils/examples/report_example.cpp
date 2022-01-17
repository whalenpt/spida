
#include <vector>
#include <fstream>
#include <iostream>
#include <pwutils/pwdefs.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/json.hpp>
#include <pwutils/report/dataio.hpp>

int main()
{
    int N = 20;

    // Use pw::DataIO for simple vector data file reporting
    pw::DataIO dataio("outfolder");
    std::vector<int> int_vec(N,1);

    dataio.writeFile("int_vec.dat",int_vec);
    std::vector<double> double_vec(N,1.);
    dataio.writeFile("double_vec.dat",double_vec);
    std::vector<pw::dcmplx> cplx_vec(N,0.);
    std::vector<std::complex<float>> cplx_fvec(N,0.);
    dataio.writeFile("complex_vec.dat",cplx_vec);
    std::vector<double> double_vec2(N,5.0);
    dataio.writeFile("two_double_vecs.dat",double_vec,double_vec2);

    // Use ReportData1D for data reporting with x, and y vector variables
    dat::ReportData1D<int,double> data1D("double_data",int_vec,double_vec);
    data1D.setReportMetadata(false);
    std::ofstream os;
    os << std::scientific << std::setprecision(3);
    os << data1D;

    // Use ReportComplexData1D for data reporting with x, and y when y is complex
    dat::ReportComplexData1D<int,double> report_def("complex_data",int_vec,cplx_vec);
    report_def.setReportMetadata(false);
    os << std::fixed << std::setprecision(2);
    os << report_def;

    // Use ReportComplexData1D for data reporting with x, and y when y is a float complex
    dat::ReportComplexData1D<int,float> report_fdef("complex_float_data",int_vec,cplx_fvec);
    report_fdef.setReportMetadata(false);
    os << std::fixed << std::setprecision(2);
    os << report_fdef;

    // Report doubles
    dat::ReportData1D<double,double> report_real("real_data",double_vec,double_vec);
    report_real.setReportMetadata(false);
    os << std::scientific << std::setprecision(3);
    report_real.report(os);

    // Json
    json::ReportData1D<int,double> jdata1D("double_data",int_vec,double_vec);
    jdata1D.setReportMetadata(false);
    os << std::scientific << std::setprecision(3);
    os << jdata1D;

    // Use ReportComplexData1D for data reporting with x, and y when y is complex
    json::ReportComplexData1D<int,double> jreport_def("complex_data",int_vec,cplx_vec);
    jreport_def.setReportMetadata(false);
    os << std::scientific << std::setprecision(3);
    os << jreport_def;

    // Use ReportComplexData1D for data reporting with x, and y when y is complex
    json::ReportComplexData1D<int,float> jreport_fdef("complex_float_data",int_vec,cplx_fvec);
    jreport_fdef.setReportMetadata(false);
    os << std::scientific << std::setprecision(3);
    os << jreport_fdef;

    // Report doubles
    json::ReportData1D<double,double> jreport_real("real_data",double_vec,double_vec);
    jreport_real.setReportMetadata(false);
    os << std::scientific << std::setprecision(3);
    jreport_real.report(os);

    // Report 2D 
    int nx = 10;
    int ny = 20;
    std::vector<double> x(nx,1.3121);
    std::vector<double> y(ny,4.7833203);
    std::vector<double> z(nx*ny,0.0);

    dat::ReportData2D<double,double,double> rp_2D("XY2D",x,y,z);
    rp_2D.setReportMetadata(false);
    // set file output directory
    rp_2D.setDirPath("outfolder2D");
    os << std::scientific << std::setprecision(8);
    rp_2D.report(os);

    // Report 2D complex
    std::vector<pw::dcmplx> z2(nx*ny,pw::dcmplx(1.0,2.0));
    dat::ReportComplexData2D<double,double,double> rp_2DC("XY2DC",x,y,z2);
    // set file output directory
    rp_2DC.setDirPath("outfolder2D");
    rp_2DC.report(os);
    os.close();

    json::ReportData2D<double,double,double> jrp_2D("XY2D",x,y,z);
    jrp_2D.setReportMetadata(false);
    // set file output directory
    jrp_2D.setDirPath("outfolder2D");
    os << std::scientific << std::setprecision(8);
    jrp_2D.report(os);

    // Report 2D complex
    json::ReportComplexData2D<double,double,double> jrp_2DC("XY2DC",x,y,z2);
    // set file output directory
    jrp_2DC.setDirPath("outfolder2D");
    jrp_2DC.report(os);
    os.close();


    return 0;

}







