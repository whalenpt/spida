/// @file utils_tests.cpp
/// @brief GoogleTest suite for the utils library: math, strings, report,
///        stats, constants/defs.

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <utils/constants.h>
#include <utils/defs.h>
#include <utils/except.h>
#include <utils/indexing.hpp>
#include <utils/math.hpp>
#include <utils/report.hpp>
#include <utils/stats.h>
#include <utils/strings.h>
#include <utils/threads.h>

namespace fs = std::filesystem;
using dcmplx = detail::dcmplx;

// ============================================================
//  CONSTANTS / DEFS
// ============================================================

TEST(CONSTANTS_TEST, PI_NEAR_STANDARD_VALUE)
{
    EXPECT_NEAR(detail::PI, std::numbers::pi, 1e-15);
}

TEST(CONSTANTS_TEST, IMAGINARY_UNIT_SQUARED_IS_MINUS_ONE)
{
    dcmplx result = detail::ii * detail::ii;
    EXPECT_NEAR(result.real(), -1.0, 1e-15);
    EXPECT_NEAR(result.imag(), 0.0, 1e-15);
}

TEST(CONSTANTS_TEST, AXIS_LABEL_STRINGS)
{
    EXPECT_STREQ(detail::XLABEL, "x");
    EXPECT_STREQ(detail::YLABEL, "y");
    EXPECT_STREQ(detail::ZLABEL, "z");
}

// ============================================================
//  MATH — integer helpers
// ============================================================

TEST(MATH_INTCEIL_TEST, EXACT_DIVISION_INT)
{
    EXPECT_EQ(detail::intceil(6, 3), 2);
    EXPECT_EQ(detail::intceil(9, 3), 3);
    EXPECT_EQ(detail::intceil(1, 1), 1);
}

TEST(MATH_INTCEIL_TEST, INEXACT_DIVISION_INT)
{
    EXPECT_EQ(detail::intceil(7, 3), 3);
    EXPECT_EQ(detail::intceil(1, 2), 1);
    EXPECT_EQ(detail::intceil(8, 3), 3);
}

TEST(MATH_INTCEIL_TEST, EXACT_DIVISION_UNSIGNED)
{
    EXPECT_EQ(detail::intceil(6u, 3u), 2u);
    EXPECT_EQ(detail::intceil(4u, 2u), 2u);
}

TEST(MATH_INTCEIL_TEST, INEXACT_DIVISION_UNSIGNED)
{
    EXPECT_EQ(detail::intceil(7u, 3u), 3u);
    EXPECT_EQ(detail::intceil(1u, 4u), 1u);
}

TEST(MATH_INTCEIL_TEST, EXACT_DIVISION_SIZE_T)
{
    EXPECT_EQ(detail::intceil(std::size_t{10}, std::size_t{5}), std::size_t{2});
}

TEST(MATH_INTCEIL_TEST, INEXACT_DIVISION_SIZE_T)
{
    EXPECT_EQ(detail::intceil(std::size_t{11}, std::size_t{5}), std::size_t{3});
    EXPECT_EQ(detail::intceil(std::size_t{1}, std::size_t{5}), std::size_t{1});
}

TEST(MATH_FACTORIAL_TEST, ZERO_AND_ONE)
{
    EXPECT_EQ(detail::factorial(0), 1);
    EXPECT_EQ(detail::factorial(1), 1);
}

TEST(MATH_FACTORIAL_TEST, SMALL_VALUES)
{
    EXPECT_EQ(detail::factorial(2), 2);
    EXPECT_EQ(detail::factorial(3), 6);
    EXPECT_EQ(detail::factorial(4), 24);
    EXPECT_EQ(detail::factorial(5), 120);
}

// ============================================================
//  MATH — type-check string helpers
// ============================================================

TEST(MATH_ISINTEGER_TEST, VALID_INTEGERS)
{
    EXPECT_TRUE(detail::isInteger("0"));
    EXPECT_TRUE(detail::isInteger("42"));
    EXPECT_TRUE(detail::isInteger("-7"));
    EXPECT_TRUE(detail::isInteger("100"));
}

TEST(MATH_ISINTEGER_TEST, INVALID_INPUTS)
{
    EXPECT_FALSE(detail::isInteger("3.14"));
    EXPECT_FALSE(detail::isInteger("abc"));
    EXPECT_FALSE(detail::isInteger(""));
    EXPECT_FALSE(detail::isInteger("1e3"));
    EXPECT_FALSE(detail::isInteger("12abc"));
}

TEST(MATH_ISDOUBLE_TEST, VALID_DOUBLES)
{
    EXPECT_TRUE(detail::isDouble("3.14"));
    EXPECT_TRUE(detail::isDouble("1e-3"));
    EXPECT_TRUE(detail::isDouble("0.0"));
    EXPECT_TRUE(detail::isDouble("-2.5"));
    EXPECT_TRUE(detail::isDouble("42")); // integers are valid doubles
}

TEST(MATH_ISDOUBLE_TEST, INVALID_INPUTS)
{
    EXPECT_FALSE(detail::isDouble("abc"));
    EXPECT_FALSE(detail::isDouble(""));
    EXPECT_FALSE(detail::isDouble("1.2abc"));
}

TEST(MATH_ISINTEGERS_TEST, SPACE_SEPARATED_ALL_INTEGERS)
{
    EXPECT_TRUE(detail::isIntegers("1 2 3"));
    EXPECT_TRUE(detail::isIntegers("0 -1 42"));
}

TEST(MATH_ISINTEGERS_TEST, SPACE_SEPARATED_MIXED_OR_INVALID)
{
    EXPECT_FALSE(detail::isIntegers("1 2.0 3"));
    EXPECT_FALSE(detail::isIntegers(""));
}

TEST(MATH_ISDOUBLES_TEST, SPACE_SEPARATED_ALL_DOUBLES)
{
    EXPECT_TRUE(detail::isDoubles("1.0 2.5 3.14"));
    EXPECT_TRUE(detail::isDoubles("1e3 -0.5 0"));
}

TEST(MATH_ISDOUBLES_TEST, SPACE_SEPARATED_INVALID)
{
    EXPECT_FALSE(detail::isDoubles("1.0 abc"));
    EXPECT_FALSE(detail::isDoubles(""));
}

// ============================================================
//  MATH — vector math (real)
// ============================================================

TEST(MATH_NORM_REAL_TEST, KNOWN_L2_NORM)
{
    std::vector<double> v{3.0, 4.0};
    EXPECT_NEAR(detail::norm(v), 5.0, 1e-12);
}

TEST(MATH_NORM_REAL_TEST, SINGLE_ELEMENT)
{
    std::vector<double> v{7.0};
    EXPECT_NEAR(detail::norm(v), 7.0, 1e-12);
}

TEST(MATH_MAX_REAL_TEST, FINDS_MAXIMUM)
{
    std::vector<double> v{1.0, 5.0, 3.0, 2.0};
    EXPECT_DOUBLE_EQ(detail::max(v), 5.0);
}

TEST(MATH_MIN_REAL_TEST, FINDS_MINIMUM)
{
    std::vector<double> v{4.0, 1.0, 3.0, 2.0};
    EXPECT_DOUBLE_EQ(detail::min(v), 1.0);
}

TEST(MATH_ARGMAX_REAL_TEST, RETURNS_INDEX_OF_MAXIMUM)
{
    std::vector<double> v{1.0, 5.0, 3.0, 2.0};
    EXPECT_EQ(detail::argmax(v), 1);
}

TEST(MATH_ARGMIN_REAL_TEST, RETURNS_INDEX_OF_MINIMUM)
{
    std::vector<double> v{4.0, 1.0, 3.0, 2.0};
    EXPECT_EQ(detail::argmin(v), 1);
}

TEST(MATH_RELATIVE_ERROR_REAL_TEST, IDENTICAL_VECTORS_GIVE_ZERO)
{
    std::vector<double> v{1.0, 2.0, 3.0};
    EXPECT_NEAR(detail::relative_error(v, v), 0.0, 1e-15);
}

TEST(MATH_RELATIVE_ERROR_REAL_TEST, DIFFERENT_VECTORS_GIVE_POSITIVE)
{
    std::vector<double> v1{1.0, 0.0, 0.0};
    std::vector<double> v2{0.0, 1.0, 0.0};
    EXPECT_GT(detail::relative_error(v1, v2), 0.0);
}

// ============================================================
//  MATH — vector math (complex)
// ============================================================

TEST(MATH_NORM_COMPLEX_TEST, KNOWN_NORM)
{
    // norm = sqrt(|z1|^2 + |z2|^2); |3+4i|=5, |0|=0 → sqrt(25)=5
    std::vector<dcmplx> v{{3.0, 4.0}, {0.0, 0.0}};
    EXPECT_NEAR(detail::norm(v), 5.0, 1e-12);
}

TEST(MATH_MAX_COMPLEX_TEST, FINDS_ELEMENT_WITH_LARGEST_ABS)
{
    std::vector<dcmplx> v{{1.0, 0.0}, {3.0, 4.0}, {2.0, 0.0}};
    dcmplx result = detail::max(v);
    EXPECT_NEAR(std::abs(result), 5.0, 1e-12);
    EXPECT_NEAR(result.real(), 3.0, 1e-12);
    EXPECT_NEAR(result.imag(), 4.0, 1e-12);
}

TEST(MATH_MIN_COMPLEX_TEST, FINDS_ELEMENT_WITH_SMALLEST_ABS)
{
    std::vector<dcmplx> v{{3.0, 4.0}, {1.0, 0.0}, {2.0, 0.0}};
    dcmplx result = detail::min(v);
    EXPECT_NEAR(std::abs(result), 1.0, 1e-12);
}

TEST(MATH_ARGMAX_COMPLEX_TEST, RETURNS_INDEX_OF_LARGEST_ABS)
{
    std::vector<dcmplx> v{{1.0, 0.0}, {3.0, 4.0}, {2.0, 0.0}};
    EXPECT_EQ(detail::argmax(v), 1);
}

TEST(MATH_ARGMIN_COMPLEX_TEST, RETURNS_INDEX_OF_SMALLEST_ABS)
{
    std::vector<dcmplx> v{{3.0, 4.0}, {0.5, 0.0}, {2.0, 0.0}};
    EXPECT_EQ(detail::argmin(v), 1);
}

TEST(MATH_RELATIVE_ERROR_COMPLEX_TEST, IDENTICAL_VECTORS_GIVE_ZERO)
{
    std::vector<dcmplx> v{{1.0, 1.0}, {2.0, -1.0}};
    EXPECT_NEAR(detail::relative_error(v, v), 0.0, 1e-15);
}

TEST(MATH_RELATIVE_ERROR_COMPLEX_TEST, DIFFERENT_VECTORS_GIVE_POSITIVE)
{
    std::vector<dcmplx> v1{{1.0, 0.0}};
    std::vector<dcmplx> v2{{0.0, 1.0}};
    EXPECT_GT(detail::relative_error(v1, v2), 0.0);
}

// ============================================================
//  STRINGS
// ============================================================

TEST(STRINGS_TRIM_TEST, STRIPS_LEADING_AND_TRAILING_SPACES)
{
    EXPECT_EQ(detail::trimString("  hello  "), "hello");
    EXPECT_EQ(detail::trimString("  "), "");
    EXPECT_EQ(detail::trimString("no spaces"), "no spaces");
}

TEST(STRINGS_TRIM_TEST, STRIPS_TABS)
{
    EXPECT_EQ(detail::trimString("\thello\t"), "hello");
}

TEST(STRINGS_EATWHITESPACE_TEST, COLLAPSES_INTERNAL_WHITESPACE)
{
    EXPECT_EQ(detail::eatWhiteSpace("a   b   c"), "a b c");
}

TEST(STRINGS_EATWHITESPACE_TEST, TRIMS_LEADING_TRAILING_FIRST)
{
    EXPECT_EQ(detail::eatWhiteSpace("  a  b  "), "a b");
}

TEST(STRINGS_PARSESTRING_TEST, SPLITS_ON_DELIMITER)
{
    auto parts = detail::parseString("a,b,c", ',');
    ASSERT_EQ(parts.size(), 3u);
    EXPECT_EQ(parts[0], "a");
    EXPECT_EQ(parts[1], "b");
    EXPECT_EQ(parts[2], "c");
}

TEST(STRINGS_PARSESTRING_TEST, NO_DELIMITER_RETURNS_SINGLE_ELEMENT)
{
    auto parts = detail::parseString("hello", ',');
    ASSERT_EQ(parts.size(), 1u);
    EXPECT_EQ(parts[0], "hello");
}

TEST(STRINGS_PARSESTRING_TEST, EMPTY_STRING_RETURNS_EMPTY_VECTOR)
{
    auto parts = detail::parseString("", ',');
    EXPECT_TRUE(parts.empty());
}

TEST(STRINGS_JOINVECTOR_TEST, JOINS_WITH_DELIMITER)
{
    std::vector<std::string> v{"a", "b", "c"};
    EXPECT_EQ(detail::joinVector(v, ','), "a,b,c");
}

TEST(STRINGS_JOINVECTOR_TEST, SINGLE_ELEMENT_NO_DELIMITER)
{
    std::vector<std::string> v{"hello"};
    EXPECT_EQ(detail::joinVector(v, ','), "hello");
}

TEST(STRINGS_JOINVECTOR_TEST, EMPTY_VECTOR_RETURNS_EMPTY)
{
    std::vector<std::string> v;
    EXPECT_EQ(detail::joinVector(v, ','), "");
}

TEST(STRINGS_ISWHITESPACE_TEST, ALL_SPACES_IS_WHITESPACE)
{
    EXPECT_TRUE(detail::isWhitespace("   "));
    EXPECT_TRUE(detail::isWhitespace(""));
}

TEST(STRINGS_ISWHITESPACE_TEST, NON_SPACE_IS_NOT_WHITESPACE)
{
    EXPECT_FALSE(detail::isWhitespace("a"));
    EXPECT_FALSE(detail::isWhitespace("  a  "));
}

TEST(STRINGS_DECOMMENT_TEST, STRIPS_FROM_HASH_TO_END)
{
    EXPECT_EQ(detail::decommentString("hello # comment"), "hello ");
    EXPECT_EQ(detail::decommentString("no comment"), "no comment");
    EXPECT_EQ(detail::decommentString("# full comment"), "");
}

TEST(STRINGS_COUNTWORDS_TEST, COUNTS_SPACE_SEPARATED_WORDS)
{
    EXPECT_EQ(detail::countWords("one two three"), 3);
    EXPECT_EQ(detail::countWords("word"), 1);
    EXPECT_EQ(detail::countWords(""), 0);
    EXPECT_EQ(detail::countWords("  spaced  "), 1);
}

TEST(STRINGS_COUNTCHARACTERS_TEST, COUNTS_OCCURRENCES)
{
    EXPECT_EQ(detail::countCharacters("aababc", 'a'), 3);
    EXPECT_EQ(detail::countCharacters("hello", 'z'), 0);
    EXPECT_EQ(detail::countCharacters("", 'a'), 0);
}

TEST(STRINGS_LOWERCASE_TEST, CONVERTS_TO_LOWER)
{
    EXPECT_EQ(detail::stringLowerCase("Hello World"), "hello world");
    EXPECT_EQ(detail::stringLowerCase("ABC"), "abc");
    EXPECT_EQ(detail::stringLowerCase(""), "");
}

TEST(STRINGS_UPPERCASE_TEST, CONVERTS_TO_UPPER)
{
    EXPECT_EQ(detail::stringUpperCase("hello world"), "HELLO WORLD");
    EXPECT_EQ(detail::stringUpperCase("abc"), "ABC");
    EXPECT_EQ(detail::stringUpperCase(""), "");
}

TEST(STRINGS_FINDSTRING_TEST, FINDS_PRESENT_SUBSTRING)
{
    // findString uses first/last char of str2 as delimiters for a prefix-match
    // approach; test with a single-char pattern contained in str
    EXPECT_TRUE(detail::findString("hello world", "hello"));
}

TEST(STRINGS_FINDSTRING_TEST, RETURNS_FALSE_FOR_ABSENT_PATTERN)
{
    EXPECT_FALSE(detail::findString("hello", "xyz"));
}

TEST(STRINGS_FINDSTRING_TEST, EMPTY_HAYSTACK_RETURNS_FALSE)
{
    EXPECT_FALSE(detail::findString("", "abc"));
}

TEST(STRINGS_SPLITSTRING_TEST, SPLITS_ON_FOUND_PATTERN)
{
    std::string first, last;
    bool found = detail::splitString("abc(xyz)def", "(xyz)", first, last);
    EXPECT_TRUE(found);
    // first is everything before the pattern's first char
    EXPECT_EQ(first, "abc");
    // last is everything after the pattern's last char
    EXPECT_EQ(last, "def");
}

TEST(STRINGS_SPLITSTRING_TEST, RETURNS_FALSE_WHEN_PATTERN_ABSENT)
{
    std::string first, last;
    bool found = detail::splitString("hello", "(xyz)", first, last);
    EXPECT_FALSE(found);
    EXPECT_EQ(first, "hello");
    EXPECT_EQ(last, "");
}

// ============================================================
//  REPORT — filesystem helpers
// ============================================================

TEST(REPORT_FS_TEST, FILEPATH_WITHOUT_REPNUM)
{
    fs::path dir = fs::temp_directory_path();
    fs::path p = detail::filePath(dir, "data", "json");
    EXPECT_EQ(p, dir / "data.json");
}

TEST(REPORT_FS_TEST, FILEPATH_WITH_REPNUM)
{
    fs::path dir = fs::temp_directory_path();
    fs::path p = detail::filePath(dir, "data", 3, "json");
    EXPECT_EQ(p, dir / "data_3.json");
}

TEST(REPORT_FS_TEST, CREATE_THEN_CLEAR_DIRECTORY)
{
    fs::path dir = fs::temp_directory_path() / "spida_create_clear_test";
    detail::createDirectory(dir);
    EXPECT_TRUE(fs::exists(dir));
    detail::clearDirectory(dir);
    EXPECT_FALSE(fs::exists(dir));
}

TEST(REPORT_FS_TEST, CREATE_DIRECTORY_NO_OVERWRITE_PRESERVES_CONTENTS)
{
    fs::path dir = fs::temp_directory_path() / "spida_nooverwrite_test";
    fs::create_directories(dir);
    // Touch a sentinel file
    std::ofstream sentinel{dir / "sentinel.txt"};
    sentinel.close();

    detail::createDirectory(dir, /*overwrite=*/false);
    EXPECT_TRUE(fs::exists(dir / "sentinel.txt"));

    fs::remove_all(dir);
}

// ============================================================
//  REPORT — ReportBase / Report1D
// ============================================================

/// @brief Fixture with a per-test temp directory that is cleaned up automatically.
class ReportTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        m_dir = fs::temp_directory_path() / "spida_utils_test";
        fs::create_directories(m_dir);
    }

    void TearDown() override
    {
        fs::remove_all(m_dir);
    }

    /// @brief Parse the JSON at path, asserting the file exists first.
    nlohmann::json readJson(const fs::path& p)
    {
        EXPECT_TRUE(fs::exists(p)) << "Expected file: " << p;
        std::ifstream ifs{p};
        return nlohmann::json::parse(ifs);
    }

    fs::path m_dir;
};

TEST_F(ReportTest, REPORT_BASE_GET_NAME)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    spida::Report1D<double, double> r{"myname", x, y};
    EXPECT_EQ(r.getName(), "myname");
}

TEST_F(ReportTest, REPORT_BASE_SET_NAME)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    spida::Report1D<double, double> r{"old", x, y};
    r.setName("new");
    EXPECT_EQ(r.getName(), "new");
}

TEST_F(ReportTest, REPORT_BASE_PATH_WITHOUT_REPNUM)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"field", x, y};
    r.setDirPath(m_dir);
    EXPECT_EQ(r.path(), m_dir / "field.json");
}

TEST_F(ReportTest, REPORT_BASE_PATH_WITH_REPNUM)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"field", x, y};
    r.setDirPath(m_dir);
    EXPECT_EQ(r.path(3), m_dir / "field_3.json");
}

TEST_F(ReportTest, REPORT1D_WRITES_XY_JSON)
{
    std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> y{1.0, 4.0, 9.0};
    spida::Report1D<double, double> r{"curve", x, y};
    r.setDirPath(m_dir);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "curve.json");
    EXPECT_EQ(j["type"].get<std::string>(), "xy");
    ASSERT_EQ(j["x"].size(), 3u);
    EXPECT_DOUBLE_EQ(j["x"][2].get<double>(), 2.0);
    ASSERT_EQ(j["y"].size(), 3u);
    EXPECT_DOUBLE_EQ(j["y"][1].get<double>(), 4.0);
}

TEST_F(ReportTest, REPORT1D_SETITEM_STRING_APPEARS_IN_META)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"meta_test", x, y};
    r.setDirPath(m_dir);
    r.setItem("label", std::string{"volts"});

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "meta_test.json");
    EXPECT_EQ(j["meta"]["label"].get<std::string>(), "volts");
}

TEST_F(ReportTest, REPORT1D_SETITEM_DOUBLE_APPEARS_IN_META)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"meta_dbl", x, y};
    r.setDirPath(m_dir);
    r.setItem("t", 1.5);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "meta_dbl.json");
    EXPECT_NEAR(j["meta"]["t"].get<double>(), 1.5, 1e-9);
}

TEST_F(ReportTest, REPORT1D_SETITEM_INT_APPEARS_IN_META)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"meta_int", x, y};
    r.setDirPath(m_dir);
    r.setItem("step", 42);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "meta_int.json");
    EXPECT_EQ(j["meta"]["step"].get<int>(), 42);
}

TEST_F(ReportTest, REPORT1D_REMOVE_ITEM_NOT_PRESENT_IN_META)
{
    std::vector<double> x{0.0};
    std::vector<double> y{0.0};
    spida::Report1D<double, double> r{"meta_rm", x, y};
    r.setDirPath(m_dir);
    r.setItem("key", std::string{"value"});
    r.removeItem("key");

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "meta_rm.json");
    EXPECT_FALSE(j["meta"].contains("key"));
}

TEST_F(ReportTest, REPORT1D_REFLECTS_MUTATED_VECTOR)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 0.0};
    spida::Report1D<double, double> r{"mutate", x, y};
    r.setDirPath(m_dir);

    // Mutate y after constructing the report (it stores a reference)
    y[0] = 99.0;
    y[1] = 100.0;

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "mutate.json");
    EXPECT_DOUBLE_EQ(j["y"][0].get<double>(), 99.0);
    EXPECT_DOUBLE_EQ(j["y"][1].get<double>(), 100.0);
}

TEST_F(ReportTest, REPORT1D_REPNUM_CREATES_NUMBERED_FILE)
{
    std::vector<double> x{0.0};
    std::vector<double> y{1.0};
    spida::Report1D<double, double> r{"numbered", x, y};
    r.setDirPath(m_dir);

    std::ofstream ofs;
    r.report(ofs, 7u);

    EXPECT_TRUE(fs::exists(m_dir / "numbered_7.json"));
}

TEST_F(ReportTest, REPORT1D_SET_DIR_PATH_CHANGES_OUTPUT_LOCATION)
{
    fs::path sub = m_dir / "subdir";
    fs::create_directories(sub);

    std::vector<double> x{0.0};
    std::vector<double> y{1.0};
    spida::Report1D<double, double> r{"moved", x, y};
    r.setDirPath(sub);

    std::ofstream ofs;
    r.report(ofs);

    EXPECT_TRUE(fs::exists(sub / "moved.json"));
}

// ============================================================
//  REPORT — ReportComplex1D
// ============================================================

TEST_F(ReportTest, REPORTCOMPLEX1D_DEFAULT_WRITES_XY_COMPLEX)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<dcmplx> y{{1.0, 2.0}, {3.0, -1.0}};
    spida::ReportComplex1D<double, double> r{"cplx", x, y};
    r.setDirPath(m_dir);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "cplx.json");
    EXPECT_EQ(j["type"].get<std::string>(), "xy_complex");
    EXPECT_DOUBLE_EQ(j["yr"][0].get<double>(), 1.0);
    EXPECT_DOUBLE_EQ(j["yi"][0].get<double>(), 2.0);
    EXPECT_DOUBLE_EQ(j["yr"][1].get<double>(), 3.0);
    EXPECT_DOUBLE_EQ(j["yi"][1].get<double>(), -1.0);
}

TEST_F(ReportTest, REPORTCOMPLEX1D_POWER_MODE_WRITES_XY_WITH_POWER_META)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<dcmplx> y{{3.0, 4.0}, {1.0, 0.0}}; // |3+4i|^2=25, |1|^2=1
    spida::ReportComplex1D<double, double> r{"power", x, y};
    r.setDirPath(m_dir);
    r.setPower(true);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "power.json");
    EXPECT_EQ(j["type"].get<std::string>(), "xy");
    EXPECT_EQ(j["meta"]["field"].get<std::string>(), "power");
    EXPECT_NEAR(j["y"][0].get<double>(), 25.0, 1e-10);
    EXPECT_NEAR(j["y"][1].get<double>(), 1.0, 1e-10);
}

TEST_F(ReportTest, REPORTCOMPLEX1D_GETPOWER_DEFAULT_FALSE)
{
    std::vector<double> x{};
    std::vector<dcmplx> y{};
    spida::ReportComplex1D<double, double> r{"p", x, y};
    EXPECT_FALSE(r.getPower());
}

TEST_F(ReportTest, REPORTCOMPLEX1D_SET_POWER_TRUE_ROUNDTRIP)
{
    std::vector<double> x{};
    std::vector<dcmplx> y{};
    spida::ReportComplex1D<double, double> r{"p", x, y};
    r.setPower(true);
    EXPECT_TRUE(r.getPower());
}

// ============================================================
//  REPORT — Report2D
// ============================================================

TEST_F(ReportTest, REPORT2D_WRITES_XYZ_JSON)
{
    // 2x3 grid: x has 2 elements, y has 3, z has 6
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0, 2.0};
    std::vector<double> z{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    spida::Report2D<double, double, double> r{"surf", x, y, z};
    r.setDirPath(m_dir);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "surf.json");
    EXPECT_EQ(j["type"].get<std::string>(), "xyz");
    EXPECT_EQ(j["x"].size(), 2u);
    EXPECT_EQ(j["y"].size(), 3u);
    // z is a 2D array: 2 rows × 3 cols
    ASSERT_EQ(j["z"].size(), 2u);
    EXPECT_EQ(j["z"][0].size(), 3u);
}

TEST_F(ReportTest, REPORT2D_STRIDEX_HALVES_X_OUTPUT)
{
    // 4x2 grid
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{0.0, 1.0};
    std::vector<double> z(8, 1.0);
    spida::Report2D<double, double, double> r{"stride2d", x, y, z};
    r.setDirPath(m_dir);
    r.setStrideX(2);

    std::ofstream ofs;
    r.report(ofs);

    nlohmann::json j = readJson(m_dir / "stride2d.json");
    // stride 2 on x=4 elements → 2 output x values
    EXPECT_EQ(j["x"].size(), 2u);
    EXPECT_EQ(j["z"].size(), 2u);
}

TEST_F(ReportTest, REPORT2D_DEFAULT_STRIDES_ARE_ONE)
{
    std::vector<double> x{};
    std::vector<double> y{};
    std::vector<double> z{};
    spida::Report2D<double, double, double> r{"s", x, y, z};
    EXPECT_EQ(r.getStrideX(), 1u);
    EXPECT_EQ(r.getStrideY(), 1u);
}

// ============================================================
//  REPORT — Track
// ============================================================

TEST_F(ReportTest, TRACK_MAX_WRITES_XY_WITH_MAX_META)
{
    std::vector<double> data{1.0, 5.0, 3.0};
    spida::Track<double> t{"maxtrack", spida::TrackType::Max, data};
    t.setDirPath(m_dir);

    t.updateTracker(0.0);
    t.updateTracker(1.0);

    std::ofstream ofs;
    t.report(ofs);

    nlohmann::json j = readJson(m_dir / "maxtrack.json");
    EXPECT_EQ(j["type"].get<std::string>(), "xy");
    EXPECT_EQ(j["meta"]["track"].get<std::string>(), "max");
    ASSERT_EQ(j["x"].size(), 2u);
    EXPECT_DOUBLE_EQ(j["y"][0].get<double>(), 5.0);
    EXPECT_DOUBLE_EQ(j["y"][1].get<double>(), 5.0);
}

TEST_F(ReportTest, TRACK_MIN_WRITES_XY_WITH_MIN_META)
{
    std::vector<double> data{4.0, 1.0, 3.0};
    spida::Track<double> t{"mintrack", spida::TrackType::Min, data};
    t.setDirPath(m_dir);

    t.updateTracker(0.5);

    std::ofstream ofs;
    t.report(ofs);

    nlohmann::json j = readJson(m_dir / "mintrack.json");
    EXPECT_EQ(j["meta"]["track"].get<std::string>(), "min");
    EXPECT_DOUBLE_EQ(j["y"][0].get<double>(), 1.0);
}

TEST_F(ReportTest, TRACK_ACCUMULATES_TIME_STEPS)
{
    std::vector<double> data{2.0, 7.0};
    spida::Track<double> t{"accum", spida::TrackType::Max, data};
    t.setDirPath(m_dir);

    t.updateTracker(0.0);
    t.updateTracker(1.0);
    t.updateTracker(2.0);

    std::ofstream ofs;
    t.report(ofs);

    nlohmann::json j = readJson(m_dir / "accum.json");
    EXPECT_EQ(j["x"].size(), 3u);
    EXPECT_EQ(j["y"].size(), 3u);
}

// ============================================================
//  REPORT — TrackComplex
// ============================================================

TEST_F(ReportTest, TRACKCOMPLEX_MAX_WRITES_MAX_POWER_META)
{
    std::vector<dcmplx> data{{3.0, 4.0}, {1.0, 0.0}}; // |3+4i|=5 wins
    spida::TrackComplex<double> tc{"cmax", spida::TrackType::Max, data};
    tc.setDirPath(m_dir);

    tc.updateTracker(0.0);

    std::ofstream ofs;
    tc.report(ofs);

    nlohmann::json j = readJson(m_dir / "cmax.json");
    EXPECT_EQ(j["meta"]["track"].get<std::string>(), "max_power");
    // norm(3+4i) = 25 (std::norm returns |z|^2)
    EXPECT_NEAR(j["y"][0].get<double>(), 25.0, 1e-10);
}

TEST_F(ReportTest, TRACKCOMPLEX_MIN_WRITES_MIN_POWER_META)
{
    std::vector<dcmplx> data{{3.0, 4.0}, {1.0, 0.0}}; // |1|=1 is min
    spida::TrackComplex<double> tc{"cmin", spida::TrackType::Min, data};
    tc.setDirPath(m_dir);

    tc.updateTracker(0.0);

    std::ofstream ofs;
    tc.report(ofs);

    nlohmann::json j = readJson(m_dir / "cmin.json");
    EXPECT_EQ(j["meta"]["track"].get<std::string>(), "min_power");
    EXPECT_NEAR(j["y"][0].get<double>(), 1.0, 1e-10);
}

// ============================================================
//  STATS — Counter
// ============================================================

TEST(STATS_COUNTER_TEST, INITIAL_COUNT_IS_ZERO_BY_DEFAULT)
{
    detail::Counter c;
    c.addCounter("hits");
    std::ostringstream oss;
    c.report(oss);
    // Report format: "  key:                           <value right-justified in 16 chars>\n"
    EXPECT_NE(oss.str().find("hits"), std::string::npos);
    EXPECT_NE(oss.str().find(" 0 "), std::string::npos);
}

TEST(STATS_COUNTER_TEST, INITIAL_COUNT_WITH_NONZERO_START)
{
    detail::Counter c;
    c.addCounter("evts", 10u);
    std::ostringstream oss;
    c.report(oss);
    EXPECT_NE(oss.str().find("10"), std::string::npos);
}

TEST(STATS_COUNTER_TEST, INCREMENT_BY_DEFAULT_AMOUNT)
{
    detail::Counter c;
    c.addCounter("x");
    c.increment("x");
    c.increment("x");
    std::ostringstream oss;
    c.report(oss);
    EXPECT_NE(oss.str().find(" 2 "), std::string::npos);
}

TEST(STATS_COUNTER_TEST, INCREMENT_BY_CUSTOM_AMOUNT)
{
    detail::Counter c;
    c.addCounter("x", 0);
    c.increment("x", 5u);
    std::ostringstream oss;
    c.report(oss);
    EXPECT_NE(oss.str().find(" 5 "), std::string::npos);
}

TEST(STATS_COUNTER_TEST, INCREMENT_UNKNOWN_KEY_ADDS_IT)
{
    detail::Counter c;
    c.increment("new_key", 3u);
    std::ostringstream oss;
    c.report(oss);
    EXPECT_NE(oss.str().find("new_key"), std::string::npos);
}

// ============================================================
//  STATS — Tracker
// ============================================================

TEST(STATS_TRACKER_TEST, ADDED_TRACKER_APPEARS_IN_REPORT)
{
    detail::Tracker t;
    t.addTracker("energy", 0.0);
    std::ostringstream oss;
    t.report(oss);
    EXPECT_NE(oss.str().find("energy"), std::string::npos);
}

TEST(STATS_TRACKER_TEST, UPDATE_CHANGES_STORED_VALUE)
{
    detail::Tracker t;
    t.addTracker("val", 0.0);
    t.updateTracker("val", 42.0);
    std::ostringstream oss;
    t.report(oss);
    EXPECT_NE(oss.str().find("4.200e+01"), std::string::npos);
}

TEST(STATS_TRACKER_TEST, MULTIPLE_UPDATES_KEEP_LAST_VALUE)
{
    detail::Tracker t;
    t.addTracker("v");
    t.updateTracker("v", 1.0);
    t.updateTracker("v", 2.0);
    t.updateTracker("v", 99.5);
    std::ostringstream oss;
    t.report(oss);
    EXPECT_NE(oss.str().find("9.950e+01"), std::string::npos);
}

// ============================================================
//  STATS — Timer
// ============================================================

TEST(STATS_TIMER_TEST, ADDED_TIMER_APPEARS_IN_REPORT)
{
    detail::Timer t;
    t.addTimer("init");
    std::ostringstream oss;
    t.report(oss);
    EXPECT_NE(oss.str().find("init"), std::string::npos);
}

// DISABLED: Timer::endTimer has an inverted map_it != end() condition (should be ==),
// so it always returns false for an added timer.  The test currently asserts the
// wrong (buggy) behavior; re-enable and update the expectation once the bug is fixed.
TEST(STATS_TIMER_TEST, DISABLED_START_THEN_END_TIMER_RETURNS_FALSE_DUE_TO_LOGIC_BUG)
{
    // Documents the actual (buggy) runtime behavior.
    detail::Timer t;
    t.addTimer("work");
    t.startTimer("work");
    bool result = t.endTimer("work");
    EXPECT_FALSE(result);
}

TEST(STATS_TIMER_TEST, END_TIMER_ON_UNKNOWN_KEY_RETURNS_FALSE)
{
    detail::Timer t;
    bool result = t.endTimer("nonexistent");
    EXPECT_FALSE(result);
}

TEST(STATS_TIMER_TEST, START_WITHOUT_ADD_DOES_NOT_CRASH)
{
    detail::Timer t;
    EXPECT_NO_THROW(t.startTimer("implicit"));
}

// ============================================================
//  STATS — StatCenter
// ============================================================

TEST(STATS_STATCENTER_TEST, DEFAULT_CONSTRUCTION_NO_THROW)
{
    EXPECT_NO_THROW(detail::StatCenter sc{});
}

TEST(STATS_STATCENTER_TEST, ADD_AND_INCREMENT_COUNTER_NO_THROW)
{
    detail::StatCenter sc;
    EXPECT_NO_THROW(sc.addCounter("steps"));
    EXPECT_NO_THROW(sc.incrementCounter("steps"));
}

TEST(STATS_STATCENTER_TEST, STAT_UPDATE_WRITES_TO_STREAM)
{
    detail::StatCenter sc{"TESTCENTER", 1u};
    sc.addCounter("n");
    sc.incrementCounter("n");
    std::ostringstream oss;
    sc.statUpdate(oss);
    EXPECT_NE(oss.str().find("TESTCENTER"), std::string::npos);
}

TEST(STATS_STATCENTER_TEST, STAT_UPDATE_OBEYS_FREQUENCY)
{
    detail::StatCenter sc{"FREQ", 3u};
    std::ostringstream oss;
    sc.statUpdate(oss); // call 1
    sc.statUpdate(oss); // call 2
    EXPECT_EQ(oss.str().find("FREQ"), std::string::npos)
        << "Should not report before frequency threshold";
    sc.statUpdate(oss); // call 3 — should fire
    EXPECT_NE(oss.str().find("FREQ"), std::string::npos);
}

TEST(STATS_STATCENTER_TEST, SET_LOG_FREQUENCY_ZERO_THROWS)
{
    detail::StatCenter sc;
    EXPECT_THROW(sc.setLogFrequency(0u), std::invalid_argument);
}

TEST(STATS_STATCENTER_TEST, SET_HEADER_APPEARS_IN_REPORT)
{
    detail::StatCenter sc;
    sc.setHeader("MYHEADER");
    std::ostringstream oss;
    sc.report(oss);
    EXPECT_NE(oss.str().find("MYHEADER"), std::string::npos);
}

// ============================================================
//  INDEXING — nearestIndex
// ============================================================

TEST(INDEXING_TEST, EXACT_MATCH_RETURNS_CORRECT_INDEX)
{
    std::vector<double> v{0.0, 1.0, 2.0, 3.0};
    EXPECT_EQ(detail::nearestIndex(v, 2.0), 2u);
}

TEST(INDEXING_TEST, VALUE_AT_MIN_RETURNS_ZERO)
{
    // val == sorted_vec[0] hits the it == cbegin() early-return path
    std::vector<double> v{0.0, 1.0, 2.0};
    EXPECT_EQ(detail::nearestIndex(v, 0.0), 0u);
}

TEST(INDEXING_TEST, VALUE_AT_MAX_RETURNS_LAST_INDEX)
{
    std::vector<double> v{0.0, 1.0, 2.0};
    EXPECT_EQ(detail::nearestIndex(v, 2.0), 2u);
}

TEST(INDEXING_TEST, VALUE_LESS_THAN_MIN_THROWS)
{
    std::vector<double> v{0.0, 1.0, 2.0};
    EXPECT_THROW(detail::nearestIndex(v, -0.5), std::out_of_range);
}

TEST(INDEXING_TEST, VALUE_GREATER_THAN_MAX_THROWS)
{
    std::vector<double> v{0.0, 1.0, 2.0};
    EXPECT_THROW(detail::nearestIndex(v, 3.5), std::out_of_range);
}

TEST(INDEXING_TEST, NEAREST_BETWEEN_TWO_POINTS)
{
    // val=1.3 is closer to index 1 (dist 0.3) than index 2 (dist 0.7)
    // val=1.7 is closer to index 2 (dist 0.3) than index 1 (dist 0.7)
    std::vector<double> v{0.0, 1.0, 2.0, 3.0};
    EXPECT_EQ(detail::nearestIndex(v, 1.3), 1u);
    EXPECT_EQ(detail::nearestIndex(v, 1.7), 2u);
}

// ============================================================
//  THREADS — ThreadManager
// ============================================================

TEST(THREADS_TEST, CONSTRUCT_ONE_THREAD)
{
    detail::ThreadManager tm{1u};
    EXPECT_EQ(tm.getNumThreads(), 1u);
}

TEST(THREADS_TEST, CONSTRUCT_FOUR_THREADS)
{
    detail::ThreadManager tm{4u};
    EXPECT_EQ(tm.getNumThreads(), 4u);
}

TEST(THREADS_TEST, CONSTRUCT_ZERO_THROWS)
{
    EXPECT_THROW(detail::ThreadManager{0u}, detail::Exception);
}

TEST(THREADS_TEST, SET_NUM_THREADS_UPDATES_COUNT)
{
    detail::ThreadManager tm{2u};
    tm.setNumThreads(8u);
    EXPECT_EQ(tm.getNumThreads(), 8u);
}

TEST(THREADS_TEST, SET_NUM_THREADS_ZERO_THROWS)
{
    detail::ThreadManager tm{1u};
    EXPECT_THROW(tm.setNumThreads(0u), detail::Exception);
}

TEST(THREADS_TEST, GET_BOUNDS_ONE_THREAD_COVERS_WHOLE_RANGE)
{
    // Single-thread path: inner loop body (i = 1..0) never executes;
    // bounds = {0, size} covers the entire range.
    detail::ThreadManager tm{1u};
    const std::vector<unsigned> bounds = tm.getBounds(50u);
    ASSERT_EQ(bounds.size(), 2u);
    EXPECT_EQ(bounds.front(), 0u);
    EXPECT_EQ(bounds.back(), 50u);
}

TEST(THREADS_TEST, GET_BOUNDS_FOUR_THREADS_MONOTONE)
{
    // Verify front==0, back==size, and each successive bound is non-decreasing.
    detail::ThreadManager tm{4u};
    const std::vector<unsigned> bounds = tm.getBounds(100u);
    ASSERT_EQ(bounds.size(), 5u);
    EXPECT_EQ(bounds.front(), 0u);
    EXPECT_EQ(bounds.back(), 100u);
    for (std::size_t i = 1; i < bounds.size(); ++i)
        EXPECT_LE(bounds[i - 1], bounds[i]) << "bounds not monotone at index " << i;
}
