#include "spida/helper/interp.h"

#include <algorithm>
#include <format>

#include <pwutils/pwexcept.h>

namespace spida {

void checkXInterp(const std::vector<double>& xinterp, const std::vector<double>& xdata)
{
    if (xinterp.empty())
        throw pw::Exception("eval(std::vector<double> xinterp) xinterp must be a non-empty vector");
    if (xinterp[0] < xdata[0])
        throw pw::Exception(std::format(
            "eval(std::vector<double> xinterp) failed: xinterp[0] is less than the class "
            "m_xvec[0], increase xinterp[0] to at least {}",
            xdata[0]));
    if (xinterp.back() > xdata.back())
        throw pw::Exception(std::format(
            "eval(std::vector<double> xinterp) failed: xinterp[size-1] is greater than the "
            "class m_xvec[m_size-1], decrease xinterp[size-1] to below or equal to {}",
            xdata.back()));
}

void checkXInterp(double xinterp, const std::vector<double>& xdata)
{
    if (xinterp < xdata[0])
        throw pw::Exception(
            std::format("eval(double xinterp) failed: xinterp is less than the class m_xvec[0], "
                        "increase xinterp[0] to at least {}",
                        xdata[0]));
    if (xinterp > xdata.back())
        throw pw::Exception(std::format(
            "eval(double xinterp) failed: xinterp is greater than the class m_xvec[m_size-1], "
            "decrease xinterp to below or equal to {}",
            xdata.back()));
}

void checkData(const std::vector<double>& xdata, const std::vector<double>& ydata)
{
    if (xdata.size() < 2)
        throw pw::Exception("Interp constructor error: xvec size must be at least 2");
    if (xdata.size() != ydata.size())
        throw pw::Exception(
            "Interp constructor error: the xvec and yvec must be of the same vector size");
}

LinearInterp::LinearInterp(const std::vector<double>& xvec, const std::vector<double>& yvec)
    : m_xvec(xvec), m_yvec(yvec)
{
    checkData(xvec, yvec);
}

void LinearInterp::computeInterpY(const std::vector<double>& xinterp,
                                  std::vector<double>& yinterp) const
{
    yinterp.clear();
    yinterp.reserve(xinterp.size());
    for (auto xval : xinterp) {
        auto it = std::lower_bound(this->m_xvec.begin(), this->m_xvec.end(), xval);
        auto idx = static_cast<std::size_t>(std::distance(this->m_xvec.begin(), it));
        if (this->m_xvec[idx] == xval) {
            yinterp.push_back(this->m_yvec[idx]);
        }
        else {
            double xa = this->m_xvec[idx - 1], xb = this->m_xvec[idx];
            double ya = this->m_yvec[idx - 1], yb = this->m_yvec[idx];
            yinterp.push_back(ya + ((yb - ya) / (xb - xa)) * (xval - xa));
        }
    }
}

std::vector<double> LinearInterp::eval(const std::vector<double>& xinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    std::vector<double> yinterp;
    this->computeInterpY(xinterp, yinterp);
    return yinterp;
}

void LinearInterp::eval(const std::vector<double>& xinterp, std::vector<double>& yinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    this->computeInterpY(xinterp, yinterp);
}

double LinearInterp::eval(double xinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    auto it = std::lower_bound(this->m_xvec.begin(), this->m_xvec.end(), xinterp);
    auto idx = static_cast<std::size_t>(std::distance(this->m_xvec.begin(), it));
    if (this->m_xvec[idx] == xinterp)
        return this->m_yvec[idx];
    double xa = this->m_xvec[idx - 1], xb = this->m_xvec[idx];
    double ya = this->m_yvec[idx - 1], yb = this->m_yvec[idx];
    return ya + ((yb - ya) / (xb - xa)) * (xinterp - xa);
}

void tridisolve(const std::vector<double>& a,
                const std::vector<double>& b,
                const std::vector<double>& c,
                const std::vector<double>& d,
                std::vector<double>& x)
{
    if (a.size() < 2)
        throw pw::Exception("tridisolve(a,b,c,d,x) error: subdiagonal 'a' size must be at least 2");
    if (b.size() < 3)
        throw pw::Exception("tridisolve(a,b,c,d,x) error: diagonal 'b' size must be at least 3");
    if (c.size() < 2)
        throw pw::Exception(
            "tridisolve(a,b,c,d,x) error: superdiagonal 'c' size must be at least 2");
    if (d.size() < 3)
        throw pw::Exception("tridisolve(a,b,c,d,x) error: R.H.S. 'd' size must be at least 3");
    if (a.size() != b.size() - 1)
        throw pw::Exception(
            std::format("tridisolve(a,b,c,d,x) error: subdiagonal 'a' size of {} must be equal to "
                        "the size of the diagonal 'b' {} minus one",
                        a.size(),
                        b.size()));
    if (c.size() != b.size() - 1)
        throw pw::Exception(std::format(
            "tridisolve(a,b,c,d,x) error: superdiagonal 'c' size of {} must be equal to "
            "the size of the diagonal 'b' {} minus one",
            c.size(),
            b.size()));
    if (b.size() != d.size())
        throw pw::Exception(
            "tridisolve(a,b,c,d,x) error: diagonal 'b' size must equal R.H.S. 'd' size");
    auto sz = b.size();
    x.assign(sz, 0.0);
    std::vector<double> bhat(sz);
    std::vector<double> dhat(sz);
    bhat[0] = b[0];
    dhat[0] = d[0];
    for (std::size_t i = 1; i < sz; i++) {
        double mu = a[i - 1] / bhat[i - 1];
        bhat[i] = b[i] - mu * c[i - 1];
        dhat[i] = d[i] - mu * dhat[i - 1];
    }
    x[sz - 1] = dhat[sz - 1] / bhat[sz - 1];
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(sz) - 2; i >= 0; i--)
        x[i] = (dhat[i] - c[i] * x[i + 1]) / bhat[i];
}

SplineInterp::SplineInterp(const std::vector<double>& xvec, const std::vector<double>& yvec)
    : m_xvec(xvec), m_yvec(yvec), m_dk(xvec.size()), m_ck(xvec.size() - 1), m_bk(xvec.size() - 1)
{
    checkData(xvec, yvec);
    if (xvec.size() < 3)
        throw pw::Exception(
            "SplineInterp constructor error: xvec size must be at least 3 for cubic spline");
    this->initializeCoefficients();
}

void SplineInterp::buildTridiSystem(const std::vector<double>& h,
                                    const std::vector<double>& delta,
                                    std::vector<double>& a,
                                    std::vector<double>& b,
                                    std::vector<double>& c,
                                    std::vector<double>& r) const
{
    auto sz = this->m_xvec.size();
    for (std::size_t i = 0; i < sz - 2; i++)
        a[i] = h[i + 1];
    a[sz - 2] = h[sz - 2] + h[sz - 3];
    b[0] = h[1];
    for (std::size_t i = 1; i < sz - 1; i++)
        b[i] = 2 * (h[i - 1] + h[i]);
    b[sz - 1] = h[sz - 3];
    c[0] = h[0] + h[1];
    for (std::size_t i = 1; i < sz - 1; i++)
        c[i] = h[i - 1];
    r[0] = ((h[0] + 2 * c[0]) * h[1] * delta[0] + h[0] * h[0] * delta[1]) / c[0];
    for (std::size_t i = 1; i < sz - 1; i++)
        r[i] = 3 * (h[i - 1] * delta[i] + h[i] * delta[i - 1]);
    r[sz - 1] = (h[sz - 2] * h[sz - 2] * delta[sz - 3] +
                 (2 * a[sz - 2] + h[sz - 2]) * h[sz - 3] * delta[sz - 2]) /
                a[sz - 2];
}

void SplineInterp::initializeCoefficients()
{
    auto sz = this->m_xvec.size();
    std::vector<double> delta(sz - 1), h(sz - 1);
    std::vector<double> a(sz - 1), b(sz), c(sz - 1), r(sz);
    for (std::size_t i = 0; i < sz - 1; i++) {
        h[i] = this->m_xvec[i + 1] - this->m_xvec[i];
        delta[i] = (this->m_yvec[i + 1] - this->m_yvec[i]) / h[i];
    }
    this->buildTridiSystem(h, delta, a, b, c, r);
    tridisolve(a, b, c, r, this->m_dk);
    for (std::size_t i = 0; i < sz - 1; i++) {
        this->m_ck[i] = (3 * delta[i] - 2 * this->m_dk[i] - this->m_dk[i + 1]) / h[i];
        this->m_bk[i] = (this->m_dk[i] - 2 * delta[i] + this->m_dk[i + 1]) / (h[i] * h[i]);
    }
}

void SplineInterp::computeInterpY(const std::vector<double>& xinterp,
                                  std::vector<double>& yinterp) const
{
    yinterp.clear();
    yinterp.reserve(xinterp.size());
    for (auto xval : xinterp) {
        auto it = std::lower_bound(this->m_xvec.begin(), this->m_xvec.end(), xval);
        auto idx = static_cast<std::size_t>(std::distance(this->m_xvec.begin(), it));
        if (this->m_xvec[idx] == xval) {
            yinterp.push_back(this->m_yvec[idx]);
        }
        else {
            double s = xval - this->m_xvec[idx - 1];
            yinterp.push_back(this->m_yvec[idx - 1] + s * this->m_dk[idx - 1] +
                              s * s * this->m_ck[idx - 1] + s * s * s * this->m_bk[idx - 1]);
        }
    }
}

std::vector<double> SplineInterp::eval(const std::vector<double>& xinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    std::vector<double> yinterp;
    this->computeInterpY(xinterp, yinterp);
    return yinterp;
}

void SplineInterp::eval(const std::vector<double>& xinterp, std::vector<double>& yinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    this->computeInterpY(xinterp, yinterp);
}

double SplineInterp::eval(double xinterp) const
{
    checkXInterp(xinterp, this->m_xvec);
    auto it = std::lower_bound(this->m_xvec.begin(), this->m_xvec.end(), xinterp);
    auto idx = static_cast<std::size_t>(std::distance(this->m_xvec.begin(), it));
    if (this->m_xvec[idx] == xinterp)
        return this->m_yvec[idx];
    double s = xinterp - this->m_xvec[idx - 1];
    return this->m_yvec[idx - 1] + s * this->m_dk[idx - 1] + s * s * this->m_ck[idx - 1] +
           s * s * s * this->m_bk[idx - 1];
}

} // namespace spida
