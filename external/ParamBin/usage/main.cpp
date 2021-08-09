
#include <ParamBin/parambin.hpp>
#include <string>

int main()
{

	ParamBin bin;
	bin.set("var1",1.3); // set a float
	bin.set("var2",10); // set an int
	bin.set("var3","string"); // set a string
	bin.setBool("var4",true); // set a bool

	ParamBin subbin;
	subbin.set("A",5.3e8);
	std::vector<std::string> str_vals;
	str_vals.push_back("first val");
	str_vals.push_back("second val");
	subbin.set("strvec",str_vals); // set a vector (of strings)
	bin.set("SUBBIN",subbin); // set a subbin of bin

	std::cout << bin << std::endl; // print the bin

	double var1 = bin.getDbl("var1"); // access a double
	int var2 = bin.getInt("var2"); // access an int
	std::string var3 = bin.getStr("var3"); // access a string
	ParamBin subbin_access = bin.getBin("SUBBIN"); // access a subbin
	std::vector<std::string> strvec_access = bin.getBin("SUBBIN").getStrVec("strvec"); // access a vector

}
