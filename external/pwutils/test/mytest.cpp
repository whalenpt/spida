
#include <gtest/gtest.h>
#include <pwutils/pwstrings.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/json.hpp>
#include <string>
#include <filesystem>

TEST(PWUTILS_TESTS,STRINGS) {
    std::string str("UPPER and lower CaSes");
	EXPECT_EQ(pw::stringLowerCase(str),"upper and lower cases"); 
	EXPECT_EQ(pw::stringUpperCase(str),"UPPER AND LOWER CASES"); 
}

TEST(PWUTILS_TESTS,REPORT_DAT){
    int N = 20;
    std::vector<int> int_vec(N,0);
    std::vector<double> double_vec(N,1.);

    dat::ReportData1D<int,double> data1D("double_vector",int_vec,double_vec);
    data1D.setReportMetadata(false);
    data1D.setPrecision(3);
    std::ofstream os(data1D.filePath(std::filesystem::temp_directory_path()).string());
    os << data1D;
    os.close();

    json::ReportData1D<int,double> jdata1D("double_vector",int_vec,double_vec);
    jdata1D.setReportMetadata(false);
    jdata1D.setPrecision(3);
    os.open(jdata1D.filePath(std::filesystem::temp_directory_path()).string());
    os << jdata1D;
    os.close();

    dat::TrackData<double> track_max("double_vector",pw::TrackType::Max,double_vec);
    os.open(track_max.filePath(std::filesystem::temp_directory_path()));
    os << data1D;
    os.close();

    EXPECT_TRUE(true);
}







