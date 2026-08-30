#include "utils/math.hpp"

#include "utils/strings.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace detail {

int intceil(int x, int y)
{
    int q = x / y + static_cast<int>((x % y != 0));
    return q;
}

unsigned int intceil(unsigned int x, unsigned int y)
{
    unsigned int q = x / y + static_cast<int>((x % y != 0));
    return q;
}

std::size_t intceil(std::size_t x, std::size_t y)
{
    return x / y + static_cast<std::size_t>(x % y != 0);
}

int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

bool isInteger(const std::string& s) noexcept
{
    size_t pos{0};
    [[maybe_unused]] int val;
    try {
        val = std::stoi(s, &pos);
    }
    catch (...) {
        return false;
    }
    if (pos != s.size()) // no characters other than integers
        return false;
    return true;
}

bool isIntegers(const std::vector<std::string>& row) noexcept
{
    if (row.empty())
        return false;
    return std::all_of(row.cbegin(), row.cend(), isInteger);
}

bool isIntegers(const std::string& s) noexcept
{
    auto str = detail::eatWhiteSpace(s);
    if (str.empty())
        return false;
    auto str_vec = detail::parseString(str, ' ');
    return isIntegers(str_vec);
}

bool isDouble(const std::string& s) noexcept
{
    size_t pos{0};
    [[maybe_unused]] double val;
    try {
        val = std::stod(s, &pos);
    }
    catch (...) {
        return false;
    }
    if (pos != s.size()) // use the whole string
        return false;
    return true;
}

bool isDoubles(const std::vector<std::string>& row) noexcept
{
    if (row.empty())
        return false;
    return std::all_of(row.cbegin(), row.cend(), isDouble);
}

bool isDoubles(const std::string& s) noexcept
{
    auto str = detail::eatWhiteSpace(s);
    if (str.empty())
        return false;
    auto str_vec = detail::parseString(str, ' ');
    return isDoubles(str_vec);
}

} // namespace detail
