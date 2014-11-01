class varPoint;
class varRegion;

#include "varPoint.h"
#include "varRegion.h"
#include "inline_variant.hpp"
#include <cstring>

varPoint::varPoint() : Point(), dimension(0), coords(0) {}
varPoint::varPoint(const coord_type* coords, uint32_t dimension) : Point(), dimension(dimension) {
    this->coords = new coord_type[dimension];
    memcpy(this->coords, coords, dimension*sizeof(coord_type));
}
varPoint::varPoint(const varPoint& p) : Point(), dimension(p.dimension) {
    coords = new coord_type[dimension];
    memcpy(coords, p.coords, dimension*sizeof(coord_type));
}
varPoint::~varPoint() {
    delete [] coords;
}
varPoint& varPoint::operator=(const varPoint& p) {
    makeDimension(p.dimension);
    memcpy(coords, p.coords, dimension*sizeof(coord_type));
    return *this;
}
bool varPoint::operator==(const varPoint& p) const {
    for (auto i = 0u; i < dimension; i++) {
        if (match(coords[i], [](REAL_TYPE) { return true; }, [](DISCRETE_TYPE) { return false; })) {
            auto thisx = boost::get<REAL_TYPE>(coords[i]);
            auto rx = boost::get<REAL_TYPE>(p.coords[i]);
            if (std::fabs(thisx - rx) > std::numeric_limits<double>::epsilon()) return false;
        } else {
            auto thisx = boost::get<DISCRETE_TYPE>(coords[i]);
            auto rx = boost::get<DISCRETE_TYPE>(p.coords[i]);
            if (thisx != rx) return false;
        }
    }
    return true;
}
varPoint* varPoint::clone() {
    return new varPoint(*this);
}
uint32_t varPoint::getByteArraySize() {
    return sizeof(uint32_t) + dimension*sizeof(coord_type);
}
void varPoint::loadFromByteArray(const byte* ptr) {
    uint32_t dim; memcpy(&dim, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    makeDimension(dim);
    memcpy(coords, ptr, dimension*sizeof(coord_type));
}
void varPoint::storeToByteArray(byte** data, uint32_t& len) {
    len = getByteArraySize();
    *data = new byte[len];
    byte* ptr = *data;
    memcpy(ptr, &dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, coords, dimension*sizeof(coord_type));
}
bool varPoint::intersectsShape(const SpatialIndex::IShape& s) const {
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return pr->containsPoint(*this);
    assert(0);
}
bool varPoint::containsShape(const SpatialIndex::IShape&) const {
    return false;
}
bool varPoint::touchesShape(const SpatialIndex::IShape& s) const {
    const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
    if (ppt != 0) {
        if (*this == *ppt) return true;
        return false;
    }
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return pr->touchesPoint(*this);
    assert(0);
}
void varPoint::getCenter(varPoint& out) const {
    out = *this;
}
uint32_t varPoint::getDimension() const {
    return dimension;
}
void varPoint::getMBR(varRegion& out) const {
    out = varRegion(coords, coords, dimension);
}
double varPoint::getArea() const {
    return 0.0;
}
double varPoint::getMinimumDistance(const SpatialIndex::IShape& in) const {
    const varPoint* ppt = dynamic_cast<const varPoint*>(&in);
    if (ppt != 0) return getMinimumDistance(*ppt);
    const varRegion* pr = dynamic_cast<const varRegion*>(&in);
    if (pr != 0) return pr->getMinimumDistance(*this);
    assert(0);
}
double varPoint::getMinimumDistance(const varPoint& p) const {
    auto ret = 0.0;
    for (auto i = 0u; i < dimension; i++) {
        ret += std::pow(match(coords[i], [&](REAL_TYPE px) {return px - boost::get<REAL_TYPE>(p.coords[i]);}, [&](DISCRETE_TYPE px) {return (REAL_TYPE)(px - boost::get<DISCRETE_TYPE>(p.coords[i]));}), 2.0);
    }
    return std::sqrt(ret);
}
/*coord_type varPoint::getCoordinate(uint32_t index) const {
    return coords[index];
}*/
void varPoint::makeInfinite(uint32_t dimension) {
    makeDimension(dimension);
    for (auto i = 0u; i < dimension; i++) {
        coords[i] = std::numeric_limits<DISCRETE_TYPE>::max();
    }
}
void varPoint::makeDimension(uint32_t dim) {
    if (dimension != dim) {
        delete[] coords;
        dimension = dim;
        coords = new coord_type[dimension];
    }
}
