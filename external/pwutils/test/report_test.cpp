
#include <gtest/gtest.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/json.hpp>

TEST(REPORT_TEST,REPORT_DAT){
    int N = 20;
    std::vector<int> int_vec(N,0);
    std::vector<double> double_vec(N,1.);
    std::vector<double> second_dvec(N,2.);
    std::ofstream os;

    dat::ReportData1D<int,double> data1("test1.json",int_vec,double_vec);
    data1.setReportMetadata(false);
    data1.setDirPath(std::filesystem::temp_directory_path());

    os << data1;
    EXPECT_TRUE(std::filesystem::exists(data1.path()));

    dat::ReportData1D<double,double> data2("double_vector",double_vec,second_dvec);
    data2.setReportMetadata(false);
    data2.setDirPath(std::filesystem::temp_directory_path());

    os << data2;
    EXPECT_TRUE(std::filesystem::exists(data2.path()));

    dat::TrackData<double> track_max("double_vector",pw::TrackType::Max,double_vec);
    track_max.setDirPath(std::filesystem::temp_directory_path());

    os << track_max;
    EXPECT_TRUE(std::filesystem::exists(track_max.path()));
}

TEST(REPORT_TEST,REPORT_JSON){
    int N = 20;
    std::vector<int> int_vec(N,0);
    std::vector<double> double_vec(N,1.);
    std::vector<double> second_dvec(N,2.);
    std::ofstream os;

    json::ReportData1D<int,double> data1("test1",int_vec,double_vec);
    data1.setReportMetadata(false);
    data1.setDirPath(std::filesystem::temp_directory_path());
    os << std::fixed << std::setprecision(4);
    os << data1;
    EXPECT_TRUE(std::filesystem::exists(data1.path()));

    json::ReportData1D<double,double> data2("double_vector",double_vec,second_dvec);
    data2.setReportMetadata(false);
    data2.setDirPath(std::filesystem::temp_directory_path());
    os << std::scientific << std::setprecision(4);
    os << data2;
    EXPECT_TRUE(std::filesystem::exists(data2.path()));

//    json::TrackData<double> track_max("double_vector",pw::TrackType::Max,double_vec);
//    os.open(track_max.filePath(std::filesystem::temp_directory_path()));
//    os << track_max;
//    os.close();

}


TEST(REPORT_TEST,REPORT_DAT_2D){
    int N1 = 10;
    int N2 = 20;
    std::vector<int> x(N1,0);
    std::vector<double> y(N2,5.0);
    std::vector<double> z(N1*N2,1.0);
    std::ofstream os;

    dat::ReportData2D<int,double,double> data1("dat2D.dat",x,y,z);
    data1.setReportMetadata(false);
    data1.setDirPath(std::filesystem::temp_directory_path());
    os << data1;
    EXPECT_TRUE(std::filesystem::exists(data1.path()));

    std::ifstream infile{data1.path()};
    int nx,ny;
    infile >> nx >> ny;
    EXPECT_EQ(nx,N1);
    EXPECT_EQ(ny,N2);
    std::vector<int> xin(nx);
    std::vector<double> yin(ny);
    std::vector<double> zin(nx*ny);
    for(auto i = 0; i < nx; i++)
        infile >> xin[i];
    for(auto j = 0; j < ny; j++)
        infile >> yin[j];
    std::string line;
    int i = 0;
    while(std::getline(infile,line)){
        std::stringstream ss(line);
        double val;
        while(ss >> val){ 
    	    zin[i] = val;
            i++;
        }
    }
    EXPECT_EQ(z.size(),zin.size());
    EXPECT_DOUBLE_EQ(z[10],1.0);
}








