#pragma once
#include "common.h"
#include "libspatialindex/include/spatialindex/SpatialIndex.h"

class varPoint : public SpatialIndex::Point {
public:
    uint32_t dimension;
    coord_type* coords;
    varPoint();
    varPoint(const coord_type* coords, uint32_t dimension);
    varPoint(const varPoint& p);
    ~varPoint();
    varPoint& operator=(const varPoint& p);
    bool operator==(const varPoint& p) const;
    varPoint* clone();
    uint32_t getByteArraySize();
    void loadFromByteArray(const byte* ptr);
    void storeToByteArray(byte** data, uint32_t& len);
    bool intersectsShape(const SpatialIndex::IShape& s) const;
    bool containsShape(const SpatialIndex::IShape&) const;
    bool touchesShape(const SpatialIndex::IShape& s) const;
    void getCenter(varPoint& out) const;
    uint32_t getDimension() const;
    void getMBR(varRegion& out) const;
    double getArea() const;
    double getMinimumDistance(const SpatialIndex::IShape& in) const ;
    double getMinimumDistance(const varPoint& p) const;
    //coord_type getCoordinate(uint32_t index) const;
    void makeInfinite(uint32_t dimension);
    void makeDimension(uint32_t dim);
};
