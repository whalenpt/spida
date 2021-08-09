
#include <vector>
#include <fstream>
#include <iostream>
#include <pwutils/report/dat.hpp>
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
    dataio.writeFile("complex_vec.dat",cplx_vec);
    std::vector<double> double_vec2(N,5.0);
    dataio.writeFile("two_double_vecs.dat",double_vec,double_vec2);

    // Use ReportData1D for data reporting with x, and y vector variables
    dat::ReportData1D<int,double> data1D("double_data",int_vec,double_vec);
    data1D.setReportMetadata(false);
    data1D.setPrecision(3);
    std::ofstream os(data1D.filePath("outfolder").string());
    os << data1D;
    os.close();

    // Use ReportComplexData1D for data reporting with x, and y when y is complex
    dat::ReportComplexData1D<int> report_def("complex_data",int_vec,cplx_vec);
    report_def.setReportMetadata(false);
    report_def.setPrecision(3);
    os.open(report_def.filePath("outfolder").string());
    os << report_def;
    os.close();

    dat::ReportData1D<double,double> report_real("real_data",double_vec,double_vec);
    report_real.setReportMetadata(false);
    report_real.setPrecision(3);
    os.open(report_real.filePath("outfolder").string());
    report_real.report(os);
    os.close();


    return 0;

}







