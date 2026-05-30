#pragma once

#include "spida/grid/grid.h"

namespace spida {

// Base interface class
class Shape {
public:
    explicit Shape([[maybe_unused]] const Grid&) noexcept {}

    virtual ~Shape() = default;
    virtual double amplitude() const = 0;
};

} // namespace spida
