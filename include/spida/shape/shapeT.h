// shapeT.h
#pragma once

#include "spida/grid/gridT.h"
#include "spida/helper/constants.h"
#include "spida/shape/shape.h"

#include <complex>
#include <stdexcept>
#include <string>

namespace spida {

// Interface class
// ShapeT in form of A0*shape((t-t_offset)/tp)*exp(-iC*(t-t_offset)^2 +
// i*slow_phase)*exp(-i*fast_phase*(t-t_offset))
//               ->  A0*shape(...)*exp(slowPhaseFactor)*exp(fastPhaseFactor)
class ShapeT : public Shape {
public:
    ShapeT(const GridT& grid, double A, double tp) : Shape(grid), m_t(grid.getT()), m_A(A), m_tp(tp)
    {
        if (tp == 0.0)
            throw std::domain_error("ShapeT: pulse width tp must be non-zero");
    }

    ~ShapeT() override = default;

    void setAmplitude(double v)
    {
        m_A = v;
    }

    // Precondition: v != 0.0  (throws std::domain_error otherwise)
    void setWidth(double v);

    void setChirp(double v)
    {
        m_chirp = v;
    }

    void setOffset(double v)
    {
        m_offset = v;
    }

    void setSlowPhase(double v)
    {
        m_slow_phase = v;
    }

    void setFastPhase(double v)
    {
        m_omega0 = v;
    }

    double amplitude() const override
    {
        return m_A;
    }

    double width() const
    {
        return m_tp;
    }

    double offset() const
    {
        return m_offset;
    }

    double chirp() const
    {
        return m_chirp;
    }

    double fastPhase() const
    {
        return m_omega0;
    }

    double slowPhase() const
    {
        return m_slow_phase;
    }

    [[nodiscard]] std::vector<dcmplx> shapeCV() const;
    [[nodiscard]] std::vector<double> shapeRV() const;
    [[nodiscard]] std::vector<dcmplx> envelope() const;

    const std::vector<double>& getT() const
    {
        return m_t;
    }

private:
    std::vector<double> m_t;
    double m_A;
    double m_tp;
    double m_offset{0.0};
    double m_chirp{0.0};
    double m_slow_phase{0.0};
    double m_omega0{0.0};

    dcmplx fastPhaseFactor(double t) const;
    dcmplx slowPhaseFactor(double t) const;
    virtual double compute(double t) const = 0;
};

class GaussT : public ShapeT {
public:
    using ShapeT::ShapeT;
    ~GaussT() override = default;

private:
    double compute(double t) const override;
};

class SechT : public ShapeT {
public:
    using ShapeT::ShapeT;
    ~SechT() override = default;

private:
    double compute(double t) const override;
};

class AiryT : public ShapeT {
public:
    AiryT(const GridT& grid, double A, double tp, double apod) : ShapeT(grid, A, tp), m_apod(apod)
    {
    }

    ~AiryT() override = default;

    void setApodization(double val)
    {
        m_apod = val;
    }

    double apodization() const
    {
        return m_apod;
    }

private:
    double compute(double t) const override;
    double m_apod{1.0};
};

class SuperGaussT : public ShapeT {
public:
    SuperGaussT(const GridT& grid, double A, double tp, double m) : ShapeT(grid, A, tp), m_M(m) {}

    ~SuperGaussT() override = default;

    void setM(double val)
    {
        m_M = val;
    }

    double M() const
    {
        return m_M;
    }

private:
    double compute(double t) const override;
    double m_M{1.0};
};

class BesselT : public ShapeT {
public:
    BesselT(const GridT& grid, double A, double tp, double apod) : ShapeT(grid, A, tp), m_apod(apod)
    {
    }

    ~BesselT() override = default;

    void setApodization(double val)
    {
        m_apod = val;
    }

    double apodization() const
    {
        return m_apod;
    }

private:
    double compute(double t) const override;
    double m_apod{1.0};
    // First zero of J_0, computed once at program start (not per instance).
    static const double s_j1;
};

} // namespace spida
