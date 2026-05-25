#pragma once

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "pwutils/basereport.h"
#include "pwutils/pwmath.hpp"

namespace pw {

// Build a JSON object from a metadataMap, converting numeric strings to numbers.
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
        } catch (...) {
            j[k] = v;
        }
    }
    return j;
}

// ---------- Report1D: (x, y) real data ----------

template<typename Tx, typename Ty>
class Report1D : public ReportData1D {
public:
    Report1D(const std::string& nm, const std::vector<Tx>& x, const std::vector<Ty>& y)
        : ReportData1D(nm), m_x(x), m_y(y) {}
    ~Report1D() override = default;
private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;
    void reportImplement(std::ofstream& os) const override {
        nlohmann::json j{
            {"type", "xy"},
            {"meta", buildMeta(getMetadata())},
            {"x", m_x},
            {"y", m_y}
        };
        os << j.dump(2);
    }
};

// Deduction guide: allow Report1D{"name", x, y} without explicit template args
template<typename Tx, typename Ty>
Report1D(const std::string&, const std::vector<Tx>&, const std::vector<Ty>&) -> Report1D<Tx, Ty>;

// ---------- ReportComplex1D: (x, complex y) ----------

template<typename Tx, typename Ty>
class ReportComplex1D : public ReportData1D {
public:
    ReportComplex1D(const std::string& nm,
                    const std::vector<Tx>& x,
                    const std::vector<std::complex<Ty>>& y)
        : ReportData1D(nm), m_x(x), m_y(y) {}
    ~ReportComplex1D() override = default;
    void setPower(bool v) {m_power = v;}
    bool getPower() const {return m_power;}
private:
    const std::vector<Tx>& m_x;
    const std::vector<std::complex<Ty>>& m_y;
    bool m_power{false};
    void reportImplement(std::ofstream& os) const override {
        nlohmann::json meta = buildMeta(getMetadata());
        if (m_power) {
            std::vector<Ty> pow_y(m_y.size());
            for (std::size_t i = 0; i < m_y.size(); i++)
                pow_y[i] = static_cast<Ty>(std::norm(m_y[i]));
            meta["field"] = "power";
            nlohmann::json j{
                {"type", "xy"},
                {"meta", meta},
                {"x", m_x},
                {"y", pow_y}
            };
            os << j.dump(2);
        } else {
            std::vector<Ty> yr(m_y.size()), yi(m_y.size());
            for (std::size_t i = 0; i < m_y.size(); i++) {
                yr[i] = m_y[i].real();
                yi[i] = m_y[i].imag();
            }
            nlohmann::json j{
                {"type", "xy_complex"},
                {"meta", meta},
                {"x", m_x},
                {"yr", yr},
                {"yi", yi}
            };
            os << j.dump(2);
        }
    }
};

// ---------- Report2D: (x, y, z) real surface data ----------

template<typename Tx, typename Ty, typename Tz>
class Report2D : public ReportData2D {
public:
    Report2D(const std::string& nm,
             const std::vector<Tx>& x,
             const std::vector<Ty>& y,
             const std::vector<Tz>& z)
        : ReportData2D(nm), m_x(x), m_y(y), m_z(z) {}
    ~Report2D() override = default;
private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;
    const std::vector<Tz>& m_z;
    void reportImplement(std::ofstream& os) const override {
        assert(m_x.size() * m_y.size() == m_z.size());
        unsigned sx = getStrideX();
        unsigned sy = getStrideY();
        unsigned nx = pw::intceil(m_x.size(), sx);
        unsigned ny = pw::intceil(m_y.size(), sy);
        std::vector<Tx> xs; xs.reserve(nx);
        std::vector<Ty> ys; ys.reserve(ny);
        for (std::size_t i = 0; i < m_x.size(); i += sx) xs.push_back(m_x[i]);
        for (std::size_t j = 0; j < m_y.size(); j += sy) ys.push_back(m_y[j]);
        nlohmann::json zj = nlohmann::json::array();
        for (std::size_t i = 0; i < m_x.size(); i += sx) {
            nlohmann::json row = nlohmann::json::array();
            for (std::size_t j = 0; j < m_y.size(); j += sy)
                row.push_back(m_z[i * m_y.size() + j]);
            zj.push_back(std::move(row));
        }
        nlohmann::json j{
            {"type", "xyz"},
            {"meta", buildMeta(getMetadata())},
            {"x", xs},
            {"y", ys},
            {"z", zj}
        };
        os << j.dump(2);
    }
};

template<typename Tx, typename Ty, typename Tz>
Report2D(const std::string&, const std::vector<Tx>&, const std::vector<Ty>&, const std::vector<Tz>&)
    -> Report2D<Tx, Ty, Tz>;

// ---------- ReportComplex2D: (x, y, complex z) ----------

template<typename Tx, typename Ty, typename Tz>
class ReportComplex2D : public ReportData2D {
public:
    ReportComplex2D(const std::string& nm,
                    const std::vector<Tx>& x,
                    const std::vector<Ty>& y,
                    const std::vector<std::complex<Tz>>& z)
        : ReportData2D(nm), m_x(x), m_y(y), m_z(z) {}
    ~ReportComplex2D() override = default;
    void setPower(bool v) {m_power = v;}
    bool getPower() const {return m_power;}
private:
    const std::vector<Tx>& m_x;
    const std::vector<Ty>& m_y;
    const std::vector<std::complex<Tz>>& m_z;
    bool m_power{false};
    void reportImplement(std::ofstream& os) const override {
        assert(m_x.size() * m_y.size() == m_z.size());
        unsigned sx = getStrideX();
        unsigned sy = getStrideY();
        unsigned nx = pw::intceil(m_x.size(), sx);
        unsigned ny = pw::intceil(m_y.size(), sy);
        std::vector<Tx> xs; xs.reserve(nx);
        std::vector<Ty> ys; ys.reserve(ny);
        for (std::size_t i = 0; i < m_x.size(); i += sx) xs.push_back(m_x[i]);
        for (std::size_t j = 0; j < m_y.size(); j += sy) ys.push_back(m_y[j]);

        nlohmann::json meta = buildMeta(getMetadata());
        if (m_power) {
            meta["field"] = "power";
            nlohmann::json zj = nlohmann::json::array();
            for (std::size_t i = 0; i < m_x.size(); i += sx) {
                nlohmann::json row = nlohmann::json::array();
                for (std::size_t j = 0; j < m_y.size(); j += sy)
                    row.push_back(static_cast<Tz>(std::norm(m_z[i * m_y.size() + j])));
                zj.push_back(std::move(row));
            }
            nlohmann::json j{
                {"type", "xyz"},
                {"meta", meta},
                {"x", xs},
                {"y", ys},
                {"z", zj}
            };
            os << j.dump(2);
        } else {
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
            nlohmann::json j{
                {"type", "xyz_complex"},
                {"meta", meta},
                {"x", xs},
                {"y", ys},
                {"zr", zrj},
                {"zi", zij}
            };
            os << j.dump(2);
        }
    }
};

template<typename Tx, typename Ty, typename Tz>
ReportComplex2D(const std::string&, const std::vector<Tx>&, const std::vector<Ty>&,
                const std::vector<std::complex<Tz>>&) -> ReportComplex2D<Tx, Ty, Tz>;

// ---------- Track: accumulates (t, max/min) pairs ----------

template<typename T>
class Track : public TrackData {
public:
    Track(const std::string& nm, TrackType ttype, const std::vector<T>& data)
        : TrackData(nm, ttype), m_data(data) {}
    ~Track() override = default;
    void updateTracker(double t) override {
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
    void reportImplement(std::ofstream& os) const override {
        nlohmann::json meta = buildMeta(getMetadata());
        meta["track"] = (getTrackType() == TrackType::Max) ? "max" : "min";
        nlohmann::json j{
            {"type", "xy"},
            {"meta", meta},
            {"x", m_t},
            {"y", m_vals}
        };
        os << j.dump(2);
    }
};

// ---------- TrackComplex: accumulates power of max/min complex value ----------

template<typename T>
class TrackComplex : public TrackData {
public:
    TrackComplex(const std::string& nm, TrackType ttype,
                 const std::vector<std::complex<T>>& data)
        : TrackData(nm, ttype), m_data(data) {}
    ~TrackComplex() override = default;
    void updateTracker(double t) override {
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
    void reportImplement(std::ofstream& os) const override {
        nlohmann::json meta = buildMeta(getMetadata());
        meta["track"] = (getTrackType() == TrackType::Max) ? "max_power" : "min_power";
        nlohmann::json j{
            {"type", "xy"},
            {"meta", meta},
            {"x", m_t},
            {"y", m_vals}
        };
        os << j.dump(2);
    }
};

}
