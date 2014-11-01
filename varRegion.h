#pragma once
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include <boost/variant.hpp>

#define DISCRETE_TYPE long long
#define REAL_TYPE double
typedef boost::variant<DISCRETE_TYPE, REAL_TYPE> coord_type;

class varRegion : public SpatialIndex::Region {
public:
    coord_type* low;
    coord_type* high;
    uint32_t dimension;

    varRegion();
    varRegion(const coord_type* low, const coord_type* high, uint32_t dimension);
    varRegion(const varPoint& low, const varPoint& high);
    void init(const coord_type* low, const coord_type* high, uint32_t dimension);
    ~varRegion();
    varRegion& operator=(const varRegion& r);
    bool operator==(const varRegion& r);
    varRegion* clone();
    uint32_t getByteArraySize();
    void loadFromByteArray(const byte* ptr);
    void storeToByteArray(byte** data, uint32_t& len);
    bool intersectsShape(const SpatialIndex::IShape& s) const;
    bool containsShape(const SpatialIndex::IShape &s) const;
    bool touchesShape(const SpatialIndex::IShape& s) const;
    void getCenter(varPoint& out) const;
    uint32_t getDimension();
    void getMBR(varRegion& out) const;
    double getArea() const;
    double getMinimumDistance(const SpatialIndex::IShape& s) const;
    bool intersectsRegion(const varRegion& r) const;
    bool containsRegion(const varRegion& r) const;
    // THIS FUNCTION HAS BEEN CORRECTED FROM THE ORIGINAL;
    bool touchesRegion(const varRegion& r) const;
    double getMinimumDistance(const varRegion& r) const;
    bool containsPoint(const varPoint& p) const;
    bool touchesPoint(const varPoint& p) const;
    double getMinimumDistance(const varPoint& p) const;
    varRegion getIntersectingRegion(const varRegion& r) const;
    double getIntersectingArea(const varRegion& r) const;
    double getMargin() const;
    void combineRegion(const varRegion& r);
    void combinePoint(const varPoint& p);
    void getCombinedRegion(varRegion& out, const varRegion& in) const;
    /*coord_type getLow(uint32_t index) const;
    coord_type getHigh(uint32_t index) const*/;
    void makeInfinite(uint32_t dimension);
    void makeDimension(uint32_t dimension);
};
