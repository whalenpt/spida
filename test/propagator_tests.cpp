
#include "propagator/reporthandler.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <filesystem>
#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include <utils/report.hpp>

using dcmplx = spida::dcmplx;
namespace fs = std::filesystem;

// ---- minimal concrete propagator for testing ----

class TestPropagatorCV : public spida::PropagatorCV {
public:
    explicit TestPropagatorCV(const fs::path& path, unsigned sz = 4)
        : PropagatorCV(path), m_field(sz, dcmplx(1.0))
    {
    }

    void updateFields(double) override
    {
        m_update_count++;
    }

    std::vector<dcmplx>& propagator() override
    {
        return m_field;
    }

    int updateCount() const
    {
        return m_update_count;
    }

private:
    std::vector<dcmplx> m_field;
    int m_update_count{0};
};

// ---- fixture with temp directory ----

class PropagatorTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        m_dir = fs::temp_directory_path() / "spida_propagator_test";
        fs::create_directories(m_dir);
    }

    void TearDown() override
    {
        fs::remove_all(m_dir);
    }

    fs::path m_dir;

    // Helper: make a 1D dat report backed by local data
    std::unique_ptr<spida::Report1D<double, double>> make1DReport(const std::string& name,
                                                                  const std::vector<double>& x,
                                                                  const std::vector<double>& y)
    {
        return std::make_unique<spida::Report1D<double, double>>(name, x, y);
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
    rh.addReport(std::make_unique<spida::Report1D<double, double>>("test", x, y));
    EXPECT_TRUE(rh.hasData1D());
    EXPECT_FALSE(rh.hasData2D());
    EXPECT_FALSE(rh.hasDataTrack());
}

TEST(REPORT_HANDLER_TEST, ADD_TRACK_REPORT_SETS_HAS_DATA)
{
    spida::ReportHandler rh;
    std::vector<double> data{1.0, 2.0, 3.0};
    rh.addReport(std::make_unique<spida::Track<double>>("track", spida::TrackType::Max, data));
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

TEST(PROPAGATOR_CV_TEST, REPORT_HANDLER_NO_DATA_1D_BY_DEFAULT)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_FALSE(prop.hasData1D());
}

// ============================================================
//  stepUpdate behavior tests (with file I/O via fixture)
// ============================================================

TEST_F(PropagatorTest, STEP_UPDATE_NO_REPORTS_ALWAYS_TRUE)
{
    TestPropagatorCV prop(m_dir);
    for (int i = 0; i < 10; i++)
        EXPECT_TRUE(prop.stepUpdate(static_cast<double>(i)));
    EXPECT_EQ(prop.updateCount(), 0);
}

TEST_F(PropagatorTest, STEP_UPDATE_NO_REPORTS_SKIPS_UPDATE_FIELDS)
{
    TestPropagatorCV prop(m_dir);
    (void) prop.stepUpdate(0.0);
    (void) prop.stepUpdate(1.0);
    EXPECT_EQ(prop.updateCount(), 0);
}

TEST_F(PropagatorTest, STEP_UPDATE_TRIGGERS_UPDATE_FIELDS_WITH_1D_REPORT)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));

    // steps_per_out1D defaults to 1, so every step triggers updateFields
    (void) prop.stepUpdate(0.0);
    EXPECT_EQ(prop.updateCount(), 1);
    (void) prop.stepUpdate(1.0);
    EXPECT_EQ(prop.updateCount(), 2);
}

TEST_F(PropagatorTest, STEP_UPDATE_RESPECTS_STEPS_PER_OUTPUT_1D)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setStepsPerOutput1D(3);

    (void) prop.stepUpdate(0.0); // step 1 — no report
    (void) prop.stepUpdate(1.0); // step 2 — no report
    EXPECT_EQ(prop.updateCount(), 0);
    (void) prop.stepUpdate(2.0); // step 3 — report fires
    EXPECT_EQ(prop.updateCount(), 1);
    (void) prop.stepUpdate(3.0); // step 4 — no report
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

    EXPECT_TRUE(prop.stepUpdate(0.0));  // report 1 → count=1, max=3 → true
    EXPECT_TRUE(prop.stepUpdate(1.0));  // report 2 → count=2, max=3 → true
    EXPECT_FALSE(prop.stepUpdate(2.0)); // report 3 → count=3, max=3 → false
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
    prop.addReport(std::make_unique<spida::Track<double>>("track", spida::TrackType::Max, data));
    prop.setStepsPerOutput(2);

    (void) prop.stepUpdate(0.0); // step 1 — no output
    EXPECT_EQ(prop.updateCount(), 0);
    (void) prop.stepUpdate(1.0); // step 2 — both 1D and track fire
    EXPECT_EQ(prop.updateCount(), 1);
}

TEST_F(PropagatorTest, REPORT_1D_CREATES_FILE)
{
    std::vector<double> x{0.0, 0.5, 1.0};
    std::vector<double> y{1.0, 0.5, 0.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("field", x, y));
    prop.report1D(0.0);

    // spida::Report1D writes to <dir>/<name>_<repNum>.json
    fs::path expected = m_dir / "field_0.json";
    EXPECT_TRUE(fs::exists(expected));
}

TEST_F(PropagatorTest, REPORT_TRACK_CREATES_FILE)
{
    std::vector<double> data{1.0, 2.0, 3.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(std::make_unique<spida::Track<double>>("maxtrack", spida::TrackType::Max, data));
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

// ============================================================
//  2D report surface — propagator.cpp:29-32, 136, 172-182
//  reporthandler.cpp:19-22, 34, 59-77, 99-103
// ============================================================

/// @brief Minimal 2D data used throughout the 2D-report tests.
/// x(2) × y(2) → z(4); satisfies Report2D invariant z.size()/y.size() == x.size().
namespace {
const std::vector<double> k2Dx{0.0, 1.0};
const std::vector<double> k2Dy{0.0, 1.0};
const std::vector<double> k2Dz{0.0, 0.5, 0.5, 1.0};
} // namespace

// ---- propagator.cpp:29-32 / reporthandler.cpp:19-22 ----

TEST_F(PropagatorTest, ADD_2D_REPORT_SETS_HAS_DATA)
{
    // Covers addReport(ReportData2D) in propagator.cpp and reporthandler.cpp.
    // Observable check: with stepsPerOutput2D=1 the first stepUpdate triggers updateFields.
    // No 1D or track reports added, so hasData1D() must remain false.
    TestPropagatorCV prop(m_dir);
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("surf", k2Dx, k2Dy, k2Dz));
    prop.setStepsPerOutput2D(1);
    EXPECT_FALSE(prop.hasData1D());
    (void) prop.stepUpdate(0.0);
    EXPECT_EQ(prop.updateCount(), 1);
}

// ---- propagator.cpp:136 (ready2D true branch) ----

TEST_F(PropagatorTest, STEP_UPDATE_TRIGGERS_UPDATE_FIELDS_WITH_2D_REPORT)
{
    // stepsPerOutput2D defaults to 1; step 1 satisfies 1 % 1 == 0.
    TestPropagatorCV prop(m_dir);
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("surf", k2Dx, k2Dy, k2Dz));
    prop.setStepsPerOutput2D(1);

    (void) prop.stepUpdate(0.0);
    EXPECT_EQ(prop.updateCount(), 1);
    (void) prop.stepUpdate(1.0);
    EXPECT_EQ(prop.updateCount(), 2);
}

TEST_F(PropagatorTest, STEP_UPDATE_RESPECTS_STEPS_PER_OUTPUT_2D)
{
    // stepsPerOutput2D=3: steps 1 and 2 (1%3!=0, 2%3!=0) do not fire; step 3 fires.
    TestPropagatorCV prop(m_dir);
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("surf", k2Dx, k2Dy, k2Dz));
    prop.setStepsPerOutput2D(3);

    (void) prop.stepUpdate(0.0); // step 1 — no report
    (void) prop.stepUpdate(1.0); // step 2 — no report
    EXPECT_EQ(prop.updateCount(), 0);
    (void) prop.stepUpdate(2.0); // step 3 — report fires
    EXPECT_EQ(prop.updateCount(), 1);
    (void) prop.stepUpdate(3.0); // step 4 — no report
    EXPECT_EQ(prop.updateCount(), 1);
}

TEST_F(PropagatorTest, STEP_UPDATE_RETURNS_FALSE_WHEN_2D_MAX_REACHED)
{
    // setMaxReports2D(1): after the first report2D fires, report_count2D==1 >= max==1.
    TestPropagatorCV prop(m_dir);
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("surf", k2Dx, k2Dy, k2Dz));
    prop.setMaxReports2D(1);

    bool result = prop.stepUpdate(0.0);
    EXPECT_FALSE(result);
}

// ---- propagator.cpp:172-182 (report2D body) ----

TEST_F(PropagatorTest, REPORT_2D_CREATES_FILE)
{
    // Direct call to report2D(t) must produce <name>_<repnum>.json in m_dir.
    TestPropagatorCV prop(m_dir);
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("surf2d", k2Dx, k2Dy, k2Dz));
    prop.report2D(0.0);

    const fs::path expected = m_dir / "surf2d_0.json";
    EXPECT_TRUE(fs::exists(expected));
}

TEST_F(PropagatorTest, REPORT_DISPATCHES_ALL_THREE_TYPES)
{
    // report(t) must delegate to report1D, report2D, and reportTrack when each
    // type has at least one registered report.
    std::vector<double> x1d{0.0, 1.0};
    std::vector<double> y1d{1.0, 2.0};
    std::vector<double> trackData{1.0, 2.0};

    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("r1d", x1d, y1d));
    prop.addReport(
        std::make_unique<spida::Report2D<double, double, double>>("r2d", k2Dx, k2Dy, k2Dz));
    prop.addReport(
        std::make_unique<spida::Track<double>>("rtrack", spida::TrackType::Max, trackData));
    prop.report(0.0);

    EXPECT_TRUE(fs::exists(m_dir / "r1d_0.json"));
    EXPECT_TRUE(fs::exists(m_dir / "r2d_0.json"));
    EXPECT_TRUE(fs::exists(m_dir / "rtrack.json"));
}

TEST_F(PropagatorTest, REPORT_2D_NOOP_WHEN_NO_2D_REPORTS_ADDED)
{
    // report2D with no registered 2D reports must not throw and must not
    // create any files — the early-return guard on hasData2D() protects this.
    TestPropagatorCV prop(m_dir);
    EXPECT_NO_THROW(prop.report2D(0.0));

    bool anyFile = false;
    for (auto const& entry : fs::directory_iterator(m_dir))
        anyFile = true;
    EXPECT_FALSE(anyFile);
}

// ---- reporthandler.cpp:19-22, 34, 59-77 (handler-level 2D path) ----

TEST(REPORT_HANDLER_TEST, ADD_2D_REPORT_SETS_HAS_DATA)
{
    // Directly exercises ReportHandler::addReport(ReportData2D) and hasData2D().
    spida::ReportHandler rh;
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    std::vector<double> z{0.0, 0.5, 0.5, 1.0};
    rh.addReport(std::make_unique<spida::Report2D<double, double, double>>("surf", x, y, z));
    EXPECT_TRUE(rh.hasData2D());
    EXPECT_FALSE(rh.hasData1D());
    EXPECT_FALSE(rh.hasDataTrack());
}

TEST(REPORT_HANDLER_TEST, REPORT_2D_CREATES_FILE)
{
    // Directly exercises ReportHandler::report2D(dir, rep_num).
    // Also hits the setItem 2D-defs loop (reporthandler.cpp:34) via setItem
    // called implicitly by report2D's parent propagator path; here we call
    // report2D directly so we hit lines 59-77 in isolation.
    const fs::path dir = fs::temp_directory_path() / "spida_rh_2d_test";
    fs::create_directories(dir);

    spida::ReportHandler rh;
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    std::vector<double> z{0.0, 0.5, 0.5, 1.0};
    rh.addReport(std::make_unique<spida::Report2D<double, double, double>>("rh2d", x, y, z));
    rh.report2D(dir, 0u);

    const fs::path expected = dir / "rh2d_0.json";
    EXPECT_TRUE(fs::exists(expected));

    fs::remove_all(dir);
}

TEST(REPORT_HANDLER_TEST, SET_ITEM_PROPAGATES_TO_2D_DEFS)
{
    // setItem must iterate over m_defs_2D (reporthandler.cpp:34) without crashing.
    spida::ReportHandler rh;
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    std::vector<double> z{0.0, 0.5, 0.5, 1.0};
    rh.addReport(std::make_unique<spida::Report2D<double, double, double>>("si2d", x, y, z));
    EXPECT_NO_THROW(rh.setItem("t", 3.14));
}

// ---- reporthandler.cpp:99-103 (reportData body) ----

TEST(REPORT_HANDLER_TEST, REPORT_DATA_CALLS_BOTH_1D_AND_2D)
{
    // reportData(dir, rep_num) must delegate to report1D and report2D,
    // producing one file for each registered type.
    const fs::path dir = fs::temp_directory_path() / "spida_rh_reportdata_test";
    fs::create_directories(dir);

    spida::ReportHandler rh;

    std::vector<double> x1{0.0, 1.0};
    std::vector<double> y1{1.0, 2.0};
    rh.addReport(std::make_unique<spida::Report1D<double, double>>("rd1d", x1, y1));

    std::vector<double> x2{0.0, 1.0};
    std::vector<double> y2{0.0, 1.0};
    std::vector<double> z2{0.0, 0.5, 0.5, 1.0};
    rh.addReport(std::make_unique<spida::Report2D<double, double, double>>("rd2d", x2, y2, z2));

    rh.reportData(dir, 0u);

    EXPECT_TRUE(fs::exists(dir / "rd1d_0.json"));
    EXPECT_TRUE(fs::exists(dir / "rd2d_0.json"));

    fs::remove_all(dir);
}

// ============================================================
//  P2: setLogFrequency(0) throw (propagator.cpp:95-99)
// ============================================================

TEST(PROPAGATOR_CV_TEST, SET_LOG_FREQUENCY_ZERO_THROWS)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_THROW(prop.setLogFrequency(0), std::invalid_argument);
}

TEST(PROPAGATOR_CV_TEST, SET_LOG_FREQUENCY_NONZERO_ACCEPTED)
{
    TestPropagatorCV prop(fs::temp_directory_path());
    EXPECT_NO_THROW(prop.setLogFrequency(5));
}

// ============================================================
//  P2: stepUpdate log-progress branch (propagator.cpp:140)
// ============================================================

TEST_F(PropagatorTest, STEP_UPDATE_CALLS_REPORT_STATS_WHEN_LOG_ENABLED)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setLogProgress(true);
    prop.setLogFrequency(1);

    EXPECT_NO_THROW((void) prop.stepUpdate(0.0));
    EXPECT_EQ(prop.updateCount(), 1);
}

// ============================================================
//  P2: report1D / reportTrack / report2D early-return guards
//  (propagator.cpp:162, 174, 187) — noop when no reports added
// ============================================================

TEST_F(PropagatorTest, REPORT_1D_NOOP_WHEN_NO_1D_REPORTS_ADDED)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_NO_THROW(prop.report1D(0.0));
    EXPECT_TRUE(fs::is_empty(m_dir));
}

TEST_F(PropagatorTest, REPORT_TRACK_NOOP_WHEN_NO_TRACK_REPORTS_ADDED)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_NO_THROW(prop.reportTrack(0.0));
    EXPECT_TRUE(fs::is_empty(m_dir));
}

TEST_F(PropagatorTest, REPORT_2D_NOOP_WHEN_NO_2D_REPORTS_ADDED_DIRECT_CALL)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_NO_THROW(prop.report2D(0.0));
    EXPECT_TRUE(fs::is_empty(m_dir));
}

// ============================================================
//  Cooperative cancellation (StopReason::CancelRequested)
// ============================================================

TEST_F(PropagatorTest, CANCEL_NOT_REQUESTED_BY_DEFAULT)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_FALSE(prop.cancelRequested());
    EXPECT_EQ(prop.stopReason(), spida::StopReason::None);
}

TEST_F(PropagatorTest, REQUEST_CANCEL_STOPS_NEXT_STEP_UPDATE)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_TRUE(prop.stepUpdate(0.0));
    prop.requestCancel();
    EXPECT_TRUE(prop.cancelRequested());
    EXPECT_FALSE(prop.stepUpdate(1.0));
    EXPECT_EQ(prop.stopReason(), spida::StopReason::CancelRequested);
}

TEST_F(PropagatorTest, MAX_REPORTS_REACHED_SETS_STOP_REASON)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setMaxReports1D(1);

    EXPECT_FALSE(prop.stepUpdate(0.0));
    EXPECT_EQ(prop.stopReason(), spida::StopReason::MaxReportsReached);
}

// ============================================================
//  Divergence detection (StopReason::Diverged)
// ============================================================

class DivergingPropagatorCV : public spida::PropagatorCV {
public:
    explicit DivergingPropagatorCV(const fs::path& path, unsigned sz = 4)
        : PropagatorCV(path), m_field(sz, dcmplx(1.0))
    {
    }

    void updateFields(double) override {}
    std::vector<dcmplx>& propagator() override { return m_field; }

    bool checkDiverged(double) override { return m_diverge_now; }
    void setDivergeNow(bool val) { m_diverge_now = val; }

private:
    std::vector<dcmplx> m_field;
    bool m_diverge_now{false};
};

TEST_F(PropagatorTest, CHECK_DIVERGED_DEFAULT_FALSE)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_FALSE(prop.checkDiverged(0.0));
}

TEST_F(PropagatorTest, DIVERGED_PROPAGATOR_STOPS_STEP_UPDATE)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{0.0, 1.0};
    DivergingPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y)); // steps_per_out1D=1 so updateFields runs every step

    EXPECT_TRUE(prop.stepUpdate(0.0));
    prop.setDivergeNow(true);
    EXPECT_FALSE(prop.stepUpdate(1.0));
    EXPECT_EQ(prop.stopReason(), spida::StopReason::Diverged);
}

// ============================================================
//  Progress observer (ProgressSnapshot)
// ============================================================

TEST_F(PropagatorTest, PROGRESS_OBSERVER_NOT_CALLED_WHEN_UNSET)
{
    TestPropagatorCV prop(m_dir);
    EXPECT_NO_THROW((void) prop.stepUpdate(0.0));
}

TEST_F(PropagatorTest, PROGRESS_OBSERVER_RECEIVES_SNAPSHOT_EACH_STEP)
{
    TestPropagatorCV prop(m_dir);
    prop.setFinalTime(5.0);
    std::vector<spida::ProgressSnapshot> seen;
    prop.setProgressObserver([&seen](const spida::ProgressSnapshot& s) { seen.push_back(s); });

    (void) prop.stepUpdate(1.0);
    (void) prop.stepUpdate(2.0);

    ASSERT_EQ(seen.size(), 2u);
    EXPECT_DOUBLE_EQ(seen[0].t, 1.0);
    ASSERT_TRUE(seen[0].tf.has_value());
    EXPECT_DOUBLE_EQ(*seen[0].tf, 5.0);
    EXPECT_EQ(seen[0].stepsTaken, 1u);
    EXPECT_FALSE(seen[0].currentStepSize.has_value()); // single-arg stepUpdate leaves it unset
    EXPECT_EQ(seen[1].stepsTaken, 2u);
}

// ============================================================
//  ReportHandler sink (Phase D)
// ============================================================

TEST_F(PropagatorTest, REPORT_SINK_RECEIVES_1D_PAYLOAD)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{2.0, 3.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));

    std::vector<std::string> names;
    prop.setReportSink([&names](std::string_view name, const nlohmann::json& j) {
        names.emplace_back(name);
        EXPECT_EQ(j.at("type"), "xy");
    });

    prop.report1D(0.0);
    ASSERT_EQ(names.size(), 1u);
    EXPECT_EQ(names[0], "u");
    // default behavior unchanged: file still written alongside the sink call
    EXPECT_FALSE(fs::is_empty(m_dir));
}

TEST_F(PropagatorTest, SET_WRITE_REPORT_FILES_FALSE_SKIPS_DISK_WRITE)
{
    std::vector<double> x{0.0, 1.0};
    std::vector<double> y{2.0, 3.0};
    TestPropagatorCV prop(m_dir);
    prop.addReport(make1DReport("u", x, y));
    prop.setWriteReportFiles(false);

    int sink_calls = 0;
    prop.setReportSink([&sink_calls](std::string_view, const nlohmann::json&) { sink_calls++; });

    prop.report1D(0.0);
    EXPECT_EQ(sink_calls, 1);
    EXPECT_TRUE(fs::is_empty(m_dir));
}
