#pragma once

#include "spida/grid/grid.h"

namespace spida{

// Base interface class
class Shape 
{
    public:
        explicit Shape(const Grid&) {}; 
        virtual ~Shape() = default;
        virtual double amplitude() const = 0;
};

}