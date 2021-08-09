
#include <gtest/gtest.h>
#include <ParamBin/parambin.hpp>
#include <string>

TEST(PARAMBIN_TESTS,POSITIVE_DOUBLES) {

	ParamBin bin;
	double dbl_val = 1.3;
	bin.set("var1",dbl_val); 
	EXPECT_EQ(bin.getDbl("var1"),dbl_val);

	dbl_val = 0.05;
	bin.set("leading zeros double",dbl_val); 
	EXPECT_EQ(bin.getDbl("leading zeros double"),dbl_val);

    //using template set and get functions
	bin.set<double>("double 100.067",100.067); 
	EXPECT_EQ(bin.get<double>("double 100.067"),100.067); 

	dbl_val = 1e-13;
	bin.set("scientific double",dbl_val); 
	EXPECT_EQ(bin.getDbl("scientific double"),dbl_val);

	dbl_val = 1.0e-13;
	bin.set("scientific double 2",dbl_val);
	EXPECT_EQ(bin.getDbl("scientific double 2"),dbl_val);

	dbl_val = 1.00000000001e-13;
	bin.set("scientific double 3",dbl_val); 
	EXPECT_EQ(bin.getDbl("scientific double 3"),dbl_val);

	dbl_val = 1000001.00000000001e-13;
	bin.set("scientific double 4",dbl_val); 
}

TEST(PARAMBIN_TESTS,NEGATIVE_DOUBLES) {

	ParamBin bin;
	double dbl_val = -1.3;
	bin.set("negative double",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative double"),dbl_val);

	dbl_val = -0.005;
	bin.set("negative leading zeros double",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative leading zeros double"),dbl_val);

	dbl_val = -0.00000000005;
	bin.set("negative many leading zeros double",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative many leading zeros double"),dbl_val);

	dbl_val = -100.067;
	bin.set("double -100.067",dbl_val); 
	EXPECT_EQ(bin.getDbl("double -100.067"),dbl_val);

	dbl_val = -1e-13;
	bin.set("negative scientific double",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative scientific double"),dbl_val);

	dbl_val = -1.0e-13;
	bin.set("negative scientific double 2",dbl_val);
	EXPECT_EQ(bin.getDbl("negative scientific double 2"),dbl_val);

	dbl_val = -1.00000000001e-13;
	bin.set("negative scientific double 3",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative scientific double 3"),dbl_val);

	dbl_val = -1000001.00000000001e-13;
	bin.set("negative scientific double 4",dbl_val); 
	EXPECT_EQ(bin.getDbl("negative scientific double 4"),dbl_val);
}

TEST(PARAMBIN_TESTS,INTS) {
	ParamBin bin;

	int int_val = 0;
	bin.set("zero",int_val); 
	EXPECT_EQ(bin.getInt("zero"),int_val);

	int_val = 1;
	bin.set("one",int_val); 
	EXPECT_EQ(bin.getInt("one"),int_val);
	
	int_val = 100;
	bin.set("ending zeros",int_val); 
	EXPECT_EQ(bin.getInt("ending zeros"),int_val);

	int_val = 1000000000;
	bin.set("many ending zeros",int_val); 
	EXPECT_EQ(bin.getInt("many ending zeros"),int_val);

	int_val = 1924359842;
	bin.set("large int",int_val); 
	EXPECT_EQ(bin.getInt("large int"),int_val);

	int_val = 47;
	bin.set("fourty-seven",int_val); 
	EXPECT_EQ(bin.getInt("fourty-seven"),int_val);

	int_val = -1;
	bin.set("negative one",int_val); 
	EXPECT_EQ(bin.getInt("negative one"),int_val);

	int_val = -149;
	bin.set("negative int",int_val); 
	EXPECT_EQ(bin.getInt("negative int"),int_val);

	int_val = -1924359842;
	bin.set("large negative int",int_val); 
	EXPECT_EQ(bin.getInt("large negative int"),int_val);
}

TEST(PARAMBIN_TESTS,LONGINTS) {
	ParamBin bin;

	long int int_val = 1000000000000005;
	bin.set("long int",int_val); 
	EXPECT_EQ(bin.get<long int>("long int"),int_val);

	int_val = -149000000000000;
	bin.set("long negative int",int_val); 
	EXPECT_EQ(bin.get<long int>("long negative int"),int_val);
}

TEST(PARAMBIN_TESTS,STRINGS) {
	ParamBin bin;
    std::string str_val("string");
	bin.set("var",str_val); // set a string
	EXPECT_EQ(bin.getStr("var"),str_val);
	EXPECT_EQ(bin.getStrU("var"),"STRING");
	EXPECT_EQ(bin.getStrL("var"),"string");

	str_val = "this is a long string";
	bin.set("long string",str_val); 
	EXPECT_EQ(bin.getStr("long string"),str_val);

	str_val = "1253,1444,43434";
	bin.set("string with nums",str_val); 
	EXPECT_EQ(bin.getStr("string with nums"),str_val);

	str_val = "#?&+--@,<>/";
	bin.set("non-alphanumeric string",str_val); 
	EXPECT_EQ(bin.getStr("non-alphanumeric string"),str_val);

	str_val = "string with parenthesis (short_name)";
	bin.set("parenthesis string",str_val); 
	EXPECT_EQ(bin.getStr("parenthesis string"),str_val);
}

TEST(PARAMBIN_TESTS,BOOLS) {
	ParamBin bin;
    bool bool_val = true;
	bin.setBool("var",bool_val); // set a bool
	EXPECT_TRUE(bin.getBool("var"));
	EXPECT_TRUE(bin.getBoolT("non-registered var"));
	EXPECT_FALSE(bin.getBoolF("non-registered var"));
}

TEST(PARAMBIN_TESTS,POSITIVE_FLOATS) {

	ParamBin bin;
	float float_val = 1.3;
	bin.set("var1",float_val); 
	EXPECT_EQ(bin.getFloat("var1"),float_val);

	float_val = 0.05;
	bin.set("leading zeros double",float_val); 
	EXPECT_EQ(bin.getFloat("leading zeros double"),float_val);

	float_val = 100.067;
	bin.set("double 100.067",float_val); 
	EXPECT_EQ(bin.getFloat("double 100.067"),float_val);

	float_val = 100.067;
	bin.set("double 100.067",float_val); 
	EXPECT_EQ(bin.get<float>("double 100.067"),float_val);

	float_val = 1e-13;
	bin.set("scientific double",float_val); 
	EXPECT_EQ(bin.getFloat("scientific double"),float_val);

	float_val = 1.0e-13;
	bin.set("scientific double 2",float_val);
	EXPECT_EQ(bin.getFloat("scientific double 2"),float_val);

	float_val = 1.00000000001e-13;
	bin.set("scientific double 3",float_val); 
	EXPECT_EQ(bin.getFloat("scientific double 3"),float_val);

	float_val = 1000001.00000000001e-13;
	bin.set("scientific double 4",float_val); 
}

TEST(PARAMBIN_TESTS,NEGATIVE_FLOATS) {

	ParamBin bin;
	float float_val = -1.3;
	bin.set("negative double",float_val); 
	EXPECT_EQ(bin.getFloat("negative double"),float_val);

	float_val = -0.005;
	bin.set("negative leading zeros double",float_val); 
	EXPECT_EQ(bin.getFloat("negative leading zeros double"),float_val);

	float_val = -0.00000000005;
	bin.set("negative many leading zeros double",float_val); 
	EXPECT_EQ(bin.getFloat("negative many leading zeros double"),float_val);

	float_val = -100.067;
	bin.set("double -100.067",float_val); 
	EXPECT_EQ(bin.getFloat("double -100.067"),float_val);

	float_val = -1e-13;
	bin.set("negative scientific double",float_val); 
	EXPECT_EQ(bin.getFloat("negative scientific double"),float_val);

	float_val = -1.0e-13;
	bin.set("negative scientific double 2",float_val);
	EXPECT_EQ(bin.getFloat("negative scientific double 2"),float_val);

	float_val = -1.00000000001e-13;
	bin.set("negative scientific double 3",float_val); 
	EXPECT_EQ(bin.getFloat("negative scientific double 3"),float_val);

	float_val = -1000001.00000000001e-13;
	bin.set("negative scientific double 4",float_val); 
	EXPECT_EQ(bin.getFloat("negative scientific double 4"),float_val);
}

TEST(PARAMBIN_TESTS,CHARS) {

	ParamBin bin;
	char val = 'a';
	bin.set("char val",val);
	EXPECT_EQ(bin.getChar("char val"),'a');

	bin.set<char>("char #",'#');
	EXPECT_EQ(bin.getChar("char #"),'#');
	EXPECT_EQ(bin.get<char>("char #"),'#');
}

TEST(PARAMBIN_TESTS,NUMERIC_VECTORS) {
    ParamBin bin;
    std::vector<double> dbl_vals {1.0e-4,1e4,3.14};
    bin.set("double vector",dbl_vals);
    std::vector<double> get_dbl_vec(bin.getDblVec("double vector"));
    ASSERT_EQ(dbl_vals.size(),get_dbl_vec.size());
    for(unsigned int i = 0; i < dbl_vals.size(); i++)
        EXPECT_EQ(dbl_vals[i],get_dbl_vec[i]);

    // Commas seperate data in getVec functions
    bin.set("int vector","1,2,3");
    std::vector<int> get_int_vec(bin.getVec<int>("int vector"));
    ASSERT_EQ(get_int_vec.size(),3);
    EXPECT_EQ(get_int_vec[0],1);
    EXPECT_EQ(get_int_vec[1],2);
    EXPECT_EQ(get_int_vec[2],3);
}










