#pragma once

#include "pwutils/pwdefs.h"
#include "pwutils/pwmath.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

namespace nlohmann {
template <typename T> struct adl_serializer<std::complex<T>> {
    static void to_json(json& j, const std::complex<T>& c)
    {
        j = json::array({c.real(), c.imag()});
    }
};
} // namespace nlohmann

namespace pw {

// ---------- File-system helpers ----------

inline void clearDirectory(const std::filesystem::path& dir_path)
{
    if (std::filesystem::exists(dir_path))
        std::filesystem::remove_all(dir_path);
}

inline void createDirectory(const std::filesystem::path& dir_path, bool overwrite = true)
{
    if (overwrite)
        clearDirectory(dir_path);
    if (!std::filesystem::exists(dir_path))
        std::filesystem::create_directory(dir_path);
}

[[nodiscard]] inline std::filesystem::path filePath(const std::filesystem::path& dir_path,
                                                    const std::string& nm,
                                                    int repNum,
                                                    const std::string& extension)
{
    return dir_path / (nm + "_" + std::to_string(repNum) + "." + extension);
}

[[nodiscard]] inline std::filesystem::path
filePath(const std::filesystem::path& dir_path, const std::string& nm, const std::string& extension)
{
    return dir_path / (nm + "." + extension);
}

// ---------- Base classes ----------

enum class TrackType { Max, Min };

class ReportBase {
public:
    explicit ReportBase(const std::string& nm) : m_name(nm) {}

    virtual ~ReportBase() = default;

    void report(std::ofstream& os) const
    {
        if (os.is_open())
            os.close();
        pw::createDirectory(m_dirpath, false);
        os.open(path());
        if (!os.is_open())
            throw std::runtime_error("Failed to open data file for output stream");
        reportImplement(os);
        os.close();
    }

    void report(std::ofstream& os, unsigned rep_num) const
    {
        if (os.is_open())
            os.close();
        pw::createDirectory(m_dirpath, false);
        os.open(path(static_cast<int>(rep_num)));
        if (!os.is_open())
            throw std::runtime_error("Failed to open data file for output stream");
        reportImplement(os);
        os.close();
    }

    friend std::ofstream& operator<<(std::ofstream& os, const ReportBase& r)
    {
        r.report(os);
        return os;
    }

    friend std::ofstream& operator<<(std::ofstream& os, const std::unique_ptr<ReportBase>& r)
    {
        r->report(os);
        return os;
    }

    void setName(std::string_view nm)
    {
        m_name = nm;
    }

    void setItem(const std::string& key, double val)
    {
        std::ostringstream oss;
        oss << std::setprecision(12) << val;
        m_metadata_map[key] = oss.str();
    }

    void setItem(const std::string& key, float val)
    {
        std::ostringstream oss;
        oss << std::setprecision(6) << val;
        m_metadata_map[key] = oss.str();
    }

    void setItem(const std::string& key, int val)
    {
        std::ostringstream oss;
        oss << val;
        m_metadata_map[key] = oss.str();
    }

    void setItem(const std::string& key, const std::string& nm)
    {
        m_metadata_map[key] = nm;
    }

    void removeItem(const std::string& key)
    {
        m_metadata_map.erase(key);
    }

    void setFileExtension(std::string_view ext)
    {
        m_extension = ext;
    }

    void setDirPath(const std::filesystem::path& dirpath)
    {
        pw::createDirectory(dirpath, false);
        m_dirpath = dirpath;
    }

    [[nodiscard]] const std::string& getName() const
    {
        return m_name;
    }

    [[nodiscard]] const metadataMap& getMetadata() const
    {
        return m_metadata_map;
    }

    [[nodiscard]] std::filesystem::path dirpath() const
    {
        return m_dirpath;
    }

    [[nodiscard]] std::filesystem::path path() const
    {
        return pw::filePath(m_dirpath, m_name, m_extension);
    }

    [[nodiscard]] std::filesystem::path path(int rep_num) const
    {
        return pw::filePath(m_dirpath, m_name, rep_num, m_extension);
    }

private:
    std::string m_name;
    metadataMap m_metadata_map;
    std::string m_extension{"json"};
    std::filesystem::path m_dirpath{"outfolder"};
    virtual void reportImplement(std::ofstream& os) const = 0;
};

class ReportData1D : public ReportBase {
    using ReportBase::ReportBase;
};

class ReportData2D : public ReportBase {
public:
    using ReportBase::ReportBase;
    ~ReportData2D() override = default;

    [[nodiscard]] unsigned getStrideX() const
    {
        return m_strideX;
    }

    [[nodiscard]] unsigned getStrideY() const
    {
        return m_strideY;
    }

    void setStrideX(unsigned s)
    {
        assert(s >= 1);
        m_strideX = s;
    }

    void setStrideY(unsigned s)
    {
        assert(s >= 1);
        m_strideY = s;
    }

private:
    unsigned m_strideX{1};
    unsigned m_strideY{1};
};

class TrackData : public ReportBase {
public:
    TrackData(const std::string& nm, TrackType ttype) : ReportBase(nm), m_ttype(ttype) {}

    ~TrackData() override = default;
    virtual void updateTracker(double t) = 0;

    [[nodiscard]] TrackType getTrackType() const
    {
        return m_ttype;
    }

    void setTrackType(TrackType t)
    {
        m_ttype = t;
    }

private:
    TrackType m_ttype;
};

// ---------- JSON metadata helper ----------

inline nlohmann::json buildMeta(const metadataMap& meta)
{
    nlohmann::json j = nlohmann::json::object();
    for (auto const& [k, v] : meta) {
        try {
            std::size_t pos;
            double d = std::stod(v, &pos);
            if (pos == v.size())
                j[k] = d;
            else
                j[k] = v;
        }
        catch (...) {
            j[k] = v;
        }
    }
    return j;
}

// ---------- Report1D: (x, y) real data ----------

template <typename Tx, typename Ty> class Report1D : public ReportData1D {
public:
    Report1D(const std::string& nm, const std::vector<Tx>& x, const std::vector<Ty>& y)
        : ReportData1D(nm), m_x(x), m_y(y)
    {
    }

    ~Report1D() override = default;

    Report1D(const std::string&, std::vector<Tx>&&, std::vector<Ty>&&) = delete;
    Report1D(const std::string&, std::vector<Tx>&&, const std::vector<Ty>&) = delete;
    Report1D(const std::string&, const std::vector<Tx>&, std::vector<Ty>&&) = delete;

private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;

    void reportImplement(std::ofstream& os) const override
    {
        nlohmann::json j;
        j["type"] = "xy";
        j["meta"] = buildMeta(getMetadata());
        j["x"] = m_x;
        j["y"] = m_y;
        os << j.dump(2);
    }
};

template <typename Tx, typename Ty>
Report1D(const std::string&, const std::vector<Tx>&, const std::vector<Ty>&) -> Report1D<Tx, Ty>;

// ---------- ReportComplex1D: (x, complex y) ----------

template <typename Tx, typename Ty> class ReportComplex1D : public ReportData1D {
public:
    ReportComplex1D(const std::string& nm,
                    const std::vector<Tx>& x,
                    const std::vector<std::complex<Ty>>& y)
        : ReportData1D(nm), m_x(x), m_y(y)
    {
    }

    ~ReportComplex1D() override = default;

    ReportComplex1D(const std::string&,
                    std::vector<Tx>&&,
                    std::vector<std::complex<Ty>>&&) = delete;
    ReportComplex1D(const std::string&,
                    std::vector<Tx>&&,
                    const std::vector<std::complex<Ty>>&) = delete;
    ReportComplex1D(const std::string&,
                    const std::vector<Tx>&,
                    std::vector<std::complex<Ty>>&&) = delete;

    void setPower(bool v)
    {
        m_power = v;
    }

    [[nodiscard]] bool getPower() const
    {
        return m_power;
    }

private:
    const std::vector<Tx>& m_x;
    const std::vector<std::complex<Ty>>& m_y;
    bool m_power{false};

    void reportPower(std::ofstream& os) const
    {
        std::vector<Ty> pow_y(m_y.size());
        for (std::size_t i = 0; i < m_y.size(); i++)
            pow_y[i] = static_cast<Ty>(std::norm(m_y[i]));
        nlohmann::json j;
        j["type"] = "xy";
        j["meta"] = buildMeta(getMetadata());
        j["meta"]["field"] = "power";
        j["x"] = m_x;
        j["y"] = pow_y;
        os << j.dump(2);
    }

    void reportComplex(std::ofstream& os) const
    {
        std::vector<Ty> yr(m_y.size()), yi(m_y.size());
        for (std::size_t i = 0; i < m_y.size(); i++) {
            yr[i] = m_y[i].real();
            yi[i] = m_y[i].imag();
        }
        nlohmann::json j;
        j["type"] = "xy_complex";
        j["meta"] = buildMeta(getMetadata());
        j["x"] = m_x;
        j["yr"] = yr;
        j["yi"] = yi;
        os << j.dump(2);
    }

    void reportImplement(std::ofstream& os) const override
    {
        if (m_power)
            reportPower(os);
        else
            reportComplex(os);
    }
};

// ---------- Report2D: (x, y, z) real surface data ----------

template <typename Tx, typename Ty, typename Tz> class Report2D : public ReportData2D {
public:
    Report2D(const std::string& nm,
             const std::vector<Tx>& x,
             const std::vector<Ty>& y,
             const std::vector<Tz>& z)
        : ReportData2D(nm), m_x(x), m_y(y), m_z(z)
    {
    }

    ~Report2D() override = default;

    Report2D(const std::string&, std::vector<Tx>&&, std::vector<Ty>&&, std::vector<Tz>&&) = delete;

private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;
    const std::vector<Tz>& m_z;

    void reportImplement(std::ofstream& os) const override
    {
        assert(m_y.size() == 0 || m_z.size() / m_y.size() == m_x.size());
        std::size_t sx = getStrideX(), sy = getStrideY();
        std::vector<Tx> xs;
        std::vector<Ty> ys;
        xs.reserve(pw::intceil(m_x.size(), sx));
        ys.reserve(pw::intceil(m_y.size(), sy));
        for (std::size_t i = 0; i < m_x.size(); i += sx)
            xs.push_back(m_x[i]);
        for (std::size_t j = 0; j < m_y.size(); j += sy)
            ys.push_back(m_y[j]);
        nlohmann::json zj = nlohmann::json::array();
        for (std::size_t i = 0; i < m_x.size(); i += sx) {
            nlohmann::json row = nlohmann::json::array();
            for (std::size_t j = 0; j < m_y.size(); j += sy)
                row.push_back(m_z[i * m_y.size() + j]);
            zj.push_back(std::move(row));
        }
        nlohmann::json j;
        j["type"] = "xyz";
        j["meta"] = buildMeta(getMetadata());
        j["x"] = xs;
        j["y"] = ys;
        j["z"] = zj;
        os << j.dump(2);
    }
};

template <typename Tx, typename Ty, typename Tz>
Report2D(const std::string&, const std::vector<Tx>&, const std::vector<Ty>&, const std::vector<Tz>&)
    -> Report2D<Tx, Ty, Tz>;

// ---------- ReportComplex2D: (x, y, complex z) ----------

template <typename Tx, typename Ty, typename Tz> class ReportComplex2D : public ReportData2D {
public:
    ReportComplex2D(const std::string& nm,
                    const std::vector<Tx>& x,
                    const std::vector<Ty>& y,
                    const std::vector<std::complex<Tz>>& z)
        : ReportData2D(nm), m_x(x), m_y(y), m_z(z)
    {
    }

    ~ReportComplex2D() override = default;

    ReportComplex2D(const std::string&,
                    std::vector<Tx>&&,
                    std::vector<Ty>&&,
                    std::vector<std::complex<Tz>>&&) = delete;

    void setPower(bool v)
    {
        m_power = v;
    }

    [[nodiscard]] bool getPower() const
    {
        return m_power;
    }

private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;
    const std::vector<std::complex<Tz>>& m_z;
    bool m_power{false};

    void buildAxes(std::size_t sx, std::size_t sy, std::vector<Tx>& xs, std::vector<Ty>& ys) const
    {
        xs.reserve(pw::intceil(m_x.size(), sx));
        ys.reserve(pw::intceil(m_y.size(), sy));
        for (std::size_t i = 0; i < m_x.size(); i += sx)
            xs.push_back(m_x[i]);
        for (std::size_t j = 0; j < m_y.size(); j += sy)
            ys.push_back(m_y[j]);
    }

    void reportPower(std::ofstream& os) const
    {
        std::size_t sx = getStrideX(), sy = getStrideY();
        std::vector<Tx> xs;
        std::vector<Ty> ys;
        buildAxes(sx, sy, xs, ys);
        nlohmann::json zj = nlohmann::json::array();
        for (std::size_t i = 0; i < m_x.size(); i += sx) {
            nlohmann::json row = nlohmann::json::array();
            for (std::size_t j = 0; j < m_y.size(); j += sy)
                row.push_back(static_cast<Tz>(std::norm(m_z[i * m_y.size() + j])));
            zj.push_back(std::move(row));
        }
        nlohmann::json j;
        j["type"] = "xyz";
        j["meta"] = buildMeta(getMetadata());
        j["meta"]["field"] = "power";
        j["x"] = xs;
        j["y"] = ys;
        j["z"] = zj;
        os << j.dump(2);
    }

    void reportComplex(std::ofstream& os) const
    {
        std::size_t sx = getStrideX(), sy = getStrideY();
        std::vector<Tx> xs;
        std::vector<Ty> ys;
        buildAxes(sx, sy, xs, ys);
        nlohmann::json zrj = nlohmann::json::array();
        nlohmann::json zij = nlohmann::json::array();
        for (std::size_t i = 0; i < m_x.size(); i += sx) {
            nlohmann::json rrow = nlohmann::json::array();
            nlohmann::json irow = nlohmann::json::array();
            for (std::size_t j = 0; j < m_y.size(); j += sy) {
                auto val = m_z[i * m_y.size() + j];
                rrow.push_back(val.real());
                irow.push_back(val.imag());
            }
            zrj.push_back(std::move(rrow));
            zij.push_back(std::move(irow));
        }
        nlohmann::json j;
        j["type"] = "xyz_complex";
        j["meta"] = buildMeta(getMetadata());
        j["x"] = xs;
        j["y"] = ys;
        j["zr"] = zrj;
        j["zi"] = zij;
        os << j.dump(2);
    }

    void reportImplement(std::ofstream& os) const override
    {
        assert(m_y.size() == 0 || m_z.size() / m_y.size() == m_x.size());
        if (m_power)
            reportPower(os);
        else
            reportComplex(os);
    }
};

template <typename Tx, typename Ty, typename Tz>
ReportComplex2D(const std::string&,
                const std::vector<Tx>&,
                const std::vector<Ty>&,
                const std::vector<std::complex<Tz>>&) -> ReportComplex2D<Tx, Ty, Tz>;

// ---------- Track: accumulates (t, max/min) pairs ----------

template <typename T> class Track : public TrackData {
public:
    Track(const std::string& nm, TrackType ttype, const std::vector<T>& data)
        : TrackData(nm, ttype), m_data(data)
    {
    }

    ~Track() override = default;

    Track(const std::string&, TrackType, std::vector<T>&&) = delete;

    void updateTracker(double t) override
    {
        m_t.push_back(t);
        if (getTrackType() == TrackType::Max)
            m_vals.push_back(pw::max(m_data));
        else
            m_vals.push_back(pw::min(m_data));
    }

private:
    const std::vector<T>& m_data;
    std::vector<double> m_t;
    std::vector<T> m_vals;

    void reportImplement(std::ofstream& os) const override
    {
        nlohmann::json j;
        j["type"] = "xy";
        j["meta"] = buildMeta(getMetadata());
        j["meta"]["track"] = (getTrackType() == TrackType::Max) ? "max" : "min";
        j["x"] = m_t;
        j["y"] = m_vals;
        os << j.dump(2);
    }
};

// ---------- TrackComplex: accumulates power of max/min complex value ----------

template <typename T> class TrackComplex : public TrackData {
public:
    TrackComplex(const std::string& nm, TrackType ttype, const std::vector<std::complex<T>>& data)
        : TrackData(nm, ttype), m_data(data)
    {
    }

    ~TrackComplex() override = default;

    TrackComplex(const std::string&, TrackType, std::vector<std::complex<T>>&&) = delete;

    void updateTracker(double t) override
    {
        m_t.push_back(t);
        if (getTrackType() == TrackType::Max)
            m_vals.push_back(static_cast<T>(std::norm(pw::max(m_data))));
        else
            m_vals.push_back(static_cast<T>(std::norm(pw::min(m_data))));
    }

private:
    const std::vector<std::complex<T>>& m_data;
    std::vector<double> m_t;
    std::vector<T> m_vals;

    void reportImplement(std::ofstream& os) const override
    {
        nlohmann::json j;
        j["type"] = "xy";
        j["meta"] = buildMeta(getMetadata());
        j["meta"]["track"] = (getTrackType() == TrackType::Max) ? "max_power" : "min_power";
        j["x"] = m_t;
        j["y"] = m_vals;
        os << j.dump(2);
    }
};

} // namespace pw
