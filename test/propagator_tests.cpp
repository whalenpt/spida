
#include <gtest/gtest.h>
#include <filesystem>
#include <memory>
#include <vector>
#include "spida/propagator/propagator.h"
#include "spida/propagator/reporthandler.h"
#include "spida/helper/constants.h"
#include <pwutils/report.hpp>

using dcmplx = spida::dcmplx;
namespace fs = std::filesystem;

// ---- minimal concrete propagator for testing ----

class TestPropagatorCV : public spida::PropagatorCV {
public:
    explicit TestPropagatorCV(const fs::path& path, unsigned sz = 4)
        : PropagatorCV(path), m_field(sz, dcmplx(1.0)) {}
    void updateFields(double) override { m_update_count++; }
    std::vector<dcmplx>& propagator() override { return m_field; }
    int updateCount() const { return m_update_count; }
private:
    std::vector<dcmplx> m_field;
    int m_update_count{0};
};

// ---- fixture with temp directory ----

class PropagatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        m_dir = fs::temp_directory_path() / "spida_propagator_test";
        fs::create_directories(m_dir);
    }
    void TearDown() override {
        fs::remove_all(m_dir);
    }
    fs::path m_dir;

    // Helper: make a 1D dat report backed by local data
    std::unique_ptr<pw::Report1D<double,double>> make1DReport(
        const std::string& name,
        const std::vector<double>& x,
        const std::vector<double>& y)
    {
        return std::make_unique<pw::Report1D<double,double>>(name, x, y);
    }
};

// ============================================================
//  ReportHandler tests
// ============================================================

TEST(REPORT_HANDLER_TEST, EMPTY_HAS_NO_DATA)
{
    spida::ReportHandler rh;
    EXPECT_FALSE(rh.hasData1D());
    EXPECT_FALSE(rh.hasData2D());
    EXPECT_FALSE(rh.hasDataTrack());
}

TEST(REPORT_HANDLER_TEST, ADD_1D_REPORT_SETS_HAS_DATA)
{
    spida::ReportHandler rh;
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    rh.addReport(std::make_unique<pw::Report1D<double,double>>("test", x, y));
    EXPECT_TRUE(rh.hasData1D());
    EXPECT_FALSE(rh.hasData2D());
    EXPECT_FALSE(rh.hasDataTrack());
}

TEST(REPORT_HANDLER_TEST, ADD_TRACK_REPORT_SETS_HAS_DATA)
{
    spida::ReportHandler rh;
    std::vector<double> data{1.0, 2.0, 3.0};
    rh.addReport(std::make_unique<pw::Track<double>>("track", pw::TrackType::Max, data));
    EXPECT_FALSE(rh.hasData1D());
    EXPECT_TRUE(rh.hasDataTrack());
}

TEST(REPORT_HANDLER_TEST, SET_ITEM_WITH_NO_REPORTS_NO_CRASH)
{
    spida::ReportHandler rh;
    EXPECT_NO_THROW(rh.setItem("t", 1.0));
}

// ============================================================
//  PropagatorCV configuration tests (no file I/O)
// ============================================================

TEST(PROPAGATOR_CV_TEST, DIR_PATH_STORED)
{
    fs::path p = fs::temp_directory_path() / "spida_dirpath_test";
    TestPropagatorCV prop(p);
    EXPECT_EQ(prop.dirPath(), p);
}

TEST(PROPAGATOR_CV_TEST, SET_DIR_PATH)
{
    fs::path p1 = fs::temp_directory_path() / "spida_dir1";
    fs::path p2 = fs::temp_directory_path() / "spida_dir2";
    TestPropagatorCV prop(p1);
    prop.setDirPath(p2);
    EXPECT_EQ(prop.dirPath(), p2);
}

TEST(PROPAGATOR_CV_TEST, LOG_PROGRESS_DEFAULT_FALSE)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_FALSE(prop.logProgress());
}

TEST(PROPAGATOR_CV_TEST, SET_LOG_PROGRESS_TRUE)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    prop.setLogProgress(true);
    EXPECT_TRUE(prop.logProgress());
}

TEST(PROPAGATOR_CV_TEST, SET_STEPS_PER_OUTPUT_1D_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setStepsPerOutput1D(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, SET_STEPS_PER_OUTPUT_2D_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setStepsPerOutput2D(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, SET_STEPS_PER_OUTPUT_TRACK_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setStepsPerOutputTrack(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, SET_MAX_REPORTS_1D_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setMaxReports1D(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, SET_MAX_REPORTS_2D_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setMaxReports2D(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, REPORT_HANDLER_ACCESSIBLE)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    spida::ReportHandler& rh = prop.reportHandler();
    EXPECT_FALSE(rh.hasData1D());
}

// ============================================================
//  stepUpdate behavior tests (with file I/O via fixture)
// ============================================================

TEST_F(PropagatorTest, STEP_UPDATE_NO_REPORTS_ALWAYS_TRUE)
{
    TestPropagatorCV prop(m_dir);
    for (int i = 0; i < 10; i++)
        EXPECT_TRUE(prop.stepUpdate(static_cast<double>(i)));
}

TEST_F(PropagatorTest, STEP_UPDATE_NO_REPORTS_SKIPS_UPDATE_FIELDS)
{
    TestPropagatorCV prop(m_dir);
    prop.stepUpdate(0.0);
    prop.stepUpdate(1.0);
    EXPECT_EQ(prop.updateCount(), 0);
}

TEST_F(PropagatorTest, STEP_UPDATE_TRIGGERS_UPDATE_FIELDS_WITH_1D_REPORT)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));

    // steps_per_out1D defaults to 1, so every step triggers updateFields
    prop.stepUpdate(0.0);
    EXPECT_EQ(prop.updateCount(), 1);
    prop.stepUpdate(1.0);
    EXPECT_EQ(prop.updateCount(), 2);
}

TEST_F(PropagatorTest, STEP_UPDATE_RESPECTS_STEPS_PER_OUTPUT_1D)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setStepsPerOutput1D(3);

    prop.stepUpdate(0.0);  // step 1 — no report
    prop.stepUpdate(1.0);  // step 2 — no report
    EXPECT_EQ(prop.updateCount(), 0);
    prop.stepUpdate(2.0);  // step 3 — report fires
    EXPECT_EQ(prop.updateCount(), 1);
    prop.stepUpdate(3.0);  // step 4 — no report
    EXPECT_EQ(prop.updateCount(), 1);
}

TEST_F(PropagatorTest, STEP_UPDATE_RETURNS_FALSE_WHEN_1D_MAX_REACHED)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setMaxReports1D(1);

    // First step fires a report and immediately reaches max
    bool result = prop.stepUpdate(0.0);
    EXPECT_FALSE(result);
}

TEST_F(PropagatorTest, STEP_UPDATE_RETURNS_TRUE_BEFORE_1D_MAX_REACHED)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setMaxReports1D(3);
    prop.setStepsPerOutput1D(1);

    EXPECT_TRUE(prop.stepUpdate(0.0));   // report 1 → count=1, max=3 → true
    EXPECT_TRUE(prop.stepUpdate(1.0));   // report 2 → count=2, max=3 → true
    EXPECT_FALSE(prop.stepUpdate(2.0));  // report 3 → count=3, max=3 → false
}

TEST_F(PropagatorTest, STEP_UPDATE_SET_STEPS_PER_OUTPUT_SETS_ALL)
{
    // setStepsPerOutput sets 1D, 2D, and track together — verify by checking
    // that report and track fire on the same step number
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    std::vector<double> data{1.0, 2.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.addReport(std::make_unique<pw::Track<double>>("track", pw::TrackType::Max, data));
    prop.setStepsPerOutput(2);

    prop.stepUpdate(0.0);  // step 1 — no output
    EXPECT_EQ(prop.updateCount(), 0);
    prop.stepUpdate(1.0);  // step 2 — both 1D and track fire
    EXPECT_EQ(prop.updateCount(), 1);
}

TEST_F(PropagatorTest, REPORT_1D_CREATES_FILE)
{
    std::vector<double> x{0.0, 0.5, 1.0};
    std::vector<double> y{1.0, 0.5, 0.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("field", x, y));
    prop.report1D(0.0);

    // pw::Report1D writes to <dir>/<name>_<repNum>.json
    fs::path expected = m_dir / "field_0.json";
    EXPECT_TRUE(fs::exists(expected));
}

TEST_F(PropagatorTest, REPORT_TRACK_CREATES_FILE)
{
    std::vector<double> data{1.0, 2.0, 3.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(std::make_unique<pw::Track<double>>("maxtrack", pw::TrackType::Max, data));
    prop.reportTrack(0.0);

    fs::path expected = m_dir / "maxtrack.json";
    EXPECT_TRUE(fs::exists(expected));
}

TEST_F(PropagatorTest, PROPAGATOR_FIELD_ACCESSIBLE)
{
    TestPropagatorCV prop(m_dir, 8);
    std::vector<dcmplx>& field = prop.propagator();
    EXPECT_EQ(field.size(), 8u);
    field[0] = dcmplx(2.0, 0.0);
    EXPECT_EQ(prop.propagator()[0], dcmplx(2.0, 0.0));
}
