
#include <vector>
#include <fstream>
#include <iostream>
#include <pwutils/read/dat.hpp>
#include <pwutils/read/json.hpp>
#include <pwutils/read/readfile.h>
#include <pwutils/pwdefs.h>
#include <filesystem>

void print_file_sig(std::filesystem::path& path)
{
    pw::FileSignature file_sig = pw::fileSignature(path);
    if(file_sig == pw::FileSignature::DAT)
        std::cout << path.filename().string() << " : " << "is a DAT format file" << std::endl;
    else if(file_sig == pw::FileSignature::JSON)
        std::cout << path.filename().string() << " : " << "is a JSON format file" << std::endl;
    else if(file_sig == pw::FileSignature::UNKNOWN)
        std::cout << path.filename().string() << " : " << "format is unknown" << std::endl;
    else
        std::cout << path.filename().string() << " : " << "has no format specified" << std::endl;
}

void print_data_sig(std::filesystem::path& path,pw::FileSignature file_sig)
{
    pw::DataSignature data_sig = pw::dataSignature(path,file_sig);
    if(data_sig == pw::DataSignature::XY)
        std::cout << path.filename().string() << " : " << "contains data of type XY" << std::endl;
    else if(data_sig == pw::DataSignature::XCVY)
        std::cout << path.filename().string() << " : " << "contains data of type XY_C" << std::endl;
    else if(data_sig == pw::DataSignature::XYZ)
        std::cout << path.filename().string() << " : " << "contains data of type XYZ" << std::endl;
    else if(data_sig == pw::DataSignature::XYCVZ)
        std::cout << path.filename().string() << " : " << "contains data of type XYZ_C" << std::endl;
    else if(data_sig == pw::DataSignature::UNKNOWN)
        std::cout << path.filename().string() << " : " << "contains data of unknown type" << std::endl;
    else
        std::cout << path.filename().string() << " : " << "has no data type specified" << std::endl;
}

void print_op_sig(std::filesystem::path& path,pw::FileSignature file_sig)
{
    pw::OperatorSignature op_sig = pw::operatorSignature(path,file_sig);
    if(op_sig == pw::OperatorSignature::NONE)
        std::cout << path.filename().string() << " : " << "has data that does not need to be processed with an operation" << std::endl;
    else if(op_sig == pw::OperatorSignature::LOGX)
        std::cout << path.filename().string() << " : " << "contains data that needs to be logged in the x-variable" << std::endl;
    else if(op_sig == pw::OperatorSignature::LOGY)
        std::cout << path.filename().string() << " : " << "contains data that needs to be logged in the y-variable" << std::endl;
    else if(op_sig == pw::OperatorSignature::LOGXLOGY)
        std::cout << path.filename().string() << " : " << "contains data that needs to be logged in both the x-variable and y-variable" << std::endl;
    else
        std::cout << path.filename().string() << " : " << "no operation signature specified" << std::endl;
}



int main()
{
    std::filesystem::path input1(std::filesystem::current_path()/std::filesystem::path("data/T_0.json"));
    std::filesystem::path input2(std::filesystem::current_path()/std::filesystem::path("data/T_0.txt"));
    std::filesystem::path input3(std::filesystem::current_path()/std::filesystem::path("data/T_0_sig.txt"));
    std::filesystem::path input4(std::filesystem::current_path()/std::filesystem::path("data/SQ_T_0.dat"));
    std::filesystem::path input5(std::filesystem::current_path()/std::filesystem::path("data/SQ_T_0.txt"));
    std::filesystem::path input6(std::filesystem::current_path()/std::filesystem::path("data/SQ_T_0_sig.txt"));
    print_file_sig(input1);
    print_file_sig(input2);
    print_file_sig(input3);
    print_file_sig(input4);
    print_file_sig(input5);
    print_file_sig(input6);

    std::cout << std::endl << std::endl;
    print_data_sig(input1,pw::FileSignature::JSON);
    print_data_sig(input2,pw::FileSignature::JSON);
    print_data_sig(input3,pw::FileSignature::JSON);
    print_data_sig(input4,pw::FileSignature::DAT);
    print_data_sig(input5,pw::FileSignature::DAT);
    print_data_sig(input6,pw::FileSignature::DAT);
    return 0;
}







