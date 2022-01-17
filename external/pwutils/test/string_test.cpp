
#include <gtest/gtest.h>
#include <pwutils/pwstrings.h>
#include <string>
#include <filesystem>

TEST(PWUTILS_TESTS,STRINGS) {
    std::string str("UPPER and lower CaSes");
	EXPECT_EQ(pw::stringLowerCase(str),"upper and lower cases"); 
	EXPECT_EQ(pw::stringUpperCase(str),"UPPER AND LOWER CASES"); 
}






