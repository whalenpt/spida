
#ifndef SPIDA_SHAPE_H_
#define SPIDA_SHAPE_H_

namespace spida{

// Base interface class
class Shape 
{
    public:
        Shape(const Grid& grid) {}
        virtual ~Shape() {};
        virtual double amplitude() const = 0;
};


}

#endif





