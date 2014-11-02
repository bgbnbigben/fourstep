class varPoint;
class varRegion;

#include "varPoint.h"
#include "varRegion.h"
#include <cstring>
#include <boost/variant.hpp>
#include "inline_variant.hpp"

varRegion::varRegion() : SpatialIndex::Region(), low(nullptr), high(nullptr), dimension(0) {}
varRegion::varRegion(const coord_type* low, const coord_type* high, uint32_t dimension) : SpatialIndex::Region() {
    init(low, high, dimension);
}
varRegion::varRegion(const varPoint& low, const varPoint& high) : SpatialIndex::Region() {
    init(low.coords, high.coords, low.dimension);
}
void varRegion::init(const coord_type* low_, const coord_type* high_, uint32_t dimension_) {
    dimension = dimension_;
    low = new coord_type[dimension];
    high = new coord_type[dimension];
    memcpy(this->low, low_, dimension*sizeof(coord_type));
    memcpy(this->high, high_, dimension*sizeof(coord_type));
}
varRegion::~varRegion() {
    delete[] low;
    delete[] high;
}
varRegion& varRegion::operator=(const varRegion& r) {
    makeDimension(r.dimension);
    memcpy(low, r.low, dimension*sizeof(coord_type));
    memcpy(high, r.high, dimension*sizeof(coord_type));
    return *this;
}
bool varRegion::operator==(const varRegion& r) {
    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [&](REAL_TYPE p) {
                auto otherLow = boost::get<REAL_TYPE>(r.low[i]);
                auto thisHigh = boost::get<REAL_TYPE>(high[i]);
                auto otherHigh = boost::get<REAL_TYPE>(r.high[i]);
                return (std::fabs(p - otherLow) > std::numeric_limits<double>::epsilon() || std::fabs(thisHigh - otherHigh) > std::numeric_limits<double>::epsilon()); },
            [&] (DISCRETE_TYPE p) {
                return p != boost::get<DISCRETE_TYPE>(r.low[i]) || boost::get<DISCRETE_TYPE>(high[i]) != boost::get<DISCRETE_TYPE>(r.high[i]);
        })) return false;
    }
    return true;
}
varRegion* varRegion::clone() {
    return new varRegion(*this);
}
uint32_t varRegion::getByteArraySize() {
    return sizeof(uint32_t) + 2*dimension*sizeof(coord_type);
}
void varRegion::loadFromByteArray(const byte* ptr) {
    uint32_t dim;
    memcpy(&dim, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    makeDimension(dim);
    memcpy(low, ptr, dimension*sizeof(coord_type));
    ptr += dimension*sizeof(coord_type);
    memcpy(high, ptr, dimension*sizeof(coord_type));
    ptr += dimension*sizeof(coord_type);
}
void varRegion::storeToByteArray(byte** data, uint32_t& len) {
    len = getByteArraySize();
    *data = new byte[len];
    byte* ptr = *data;
    memcpy(ptr, &dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, low, dimension*sizeof(coord_type));
    ptr += dimension*sizeof(coord_type);
    memcpy(ptr, high, dimension*sizeof(coord_type));
}
bool varRegion::intersectsShape(const SpatialIndex::IShape& s) const {
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return intersectsRegion(*pr);
    const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
    if (ppt != 0) return containsPoint(*ppt);
    assert(0);
}
bool varRegion::containsShape(const SpatialIndex::IShape &s) const {
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return containsRegion(*pr);
    const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
    if (ppt != 0) return containsPoint(*ppt);
    assert(0);
}
bool varRegion::touchesShape(const SpatialIndex::IShape& s) const {
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return touchesRegion(*pr);
    const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
    if (ppt != 0) return touchesPoint(*ppt);
    assert(0);
}
void varRegion::getCenter(varPoint& out) const {
    out.makeDimension(dimension);
    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
            out.coords[i] = (boost::get<REAL_TYPE>(low[i]) + boost::get<REAL_TYPE>(high[i])) / 2.0;
        }, [&](DISCRETE_TYPE) {
            out.coords[i] = (boost::get<DISCRETE_TYPE>(low[i]) + boost::get<DISCRETE_TYPE>(high[i])) / 2.0;
        });
    }
}
uint32_t varRegion::getDimension() {
    return dimension;
}
void varRegion::getMBR(varRegion& out) const {
    out = *this;
}
double varRegion::getArea() const {
    auto area = 1.0;
    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
            area *= boost::get<REAL_TYPE>(high[i]) - boost::get<REAL_TYPE>(low[i]);
        }, [&](DISCRETE_TYPE) {
            area *= boost::get<DISCRETE_TYPE>(high[i]) - boost::get<DISCRETE_TYPE>(low[i]);
        });
    }
    return area;
}
double varRegion::getMinimumDistance(const SpatialIndex::IShape& s) const {
    const varRegion* pr = dynamic_cast<const varRegion*>(&s);
    if (pr != 0) return getMinimumDistance(*pr);
    const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
    if (ppt != 0) return getMinimumDistance(*ppt);
    assert(0);
}
bool varRegion::intersectsRegion(const varRegion& r) const {
    if (dimension != r.dimension) {
        assert(false);
    }
    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [](REAL_TYPE) {return true;}, [](DISCRETE_TYPE) {return false;})) {
            auto thislow = boost::get<REAL_TYPE>(low[i]);
            auto thishigh = boost::get<REAL_TYPE>(high[i]);
            auto rlow = boost::get<REAL_TYPE>(r.low[i]);
            auto rhigh = boost::get<REAL_TYPE>(r.high[i]);
            if (thislow > rhigh || thishigh < rlow) return false;
        } else {
            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
            auto rlow = boost::get<DISCRETE_TYPE>(r.low[i]);
            auto rhigh = boost::get<DISCRETE_TYPE>(r.high[i]);
            if (thislow > rhigh || thishigh < rlow) return false;
        }
    }
    return true;
}
bool varRegion::containsRegion(const varRegion& r) const {
    if (dimension != r.dimension) {
        assert(false);
    }

    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [](REAL_TYPE) {return true;}, [](DISCRETE_TYPE) {return false;})) {
            auto thislow = boost::get<REAL_TYPE>(low[i]);
            auto thishigh = boost::get<REAL_TYPE>(high[i]);
            auto rlow = boost::get<REAL_TYPE>(r.low[i]);
            auto rhigh = boost::get<REAL_TYPE>(r.high[i]);
            if (thislow > rlow || thishigh < rhigh) return false;
        } else {
            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
            auto rlow = boost::get<DISCRETE_TYPE>(r.low[i]);
            auto rhigh = boost::get<DISCRETE_TYPE>(r.high[i]);
            if (thislow > rlow || thishigh < rhigh) return false;
        }
    }
    return true;
}
// THIS FUNCTION HAS BEEN CORRECTED FROM THE ORIGINAL
bool varRegion::touchesRegion(const varRegion& r) const {
    if (dimension != r.dimension) {
        assert (false);
    }
    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [](REAL_TYPE) {return true;}, [](DISCRETE_TYPE) {return false;})) {
            auto thislow = boost::get<REAL_TYPE>(low[i]);
            auto thishigh = boost::get<REAL_TYPE>(high[i]);
            auto rlow = boost::get<REAL_TYPE>(r.low[i]);
            auto rhigh = boost::get<REAL_TYPE>(r.high[i]);
            if ((std::fabs(thislow - rlow) >= std::numeric_limits<double>::epsilon()) &&
                (std::fabs(thishigh - rhigh) >= std::numeric_limits<double>::epsilon()))
                return false;
        } else {
            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
            auto rlow = boost::get<DISCRETE_TYPE>(r.low[i]);
            auto rhigh = boost::get<DISCRETE_TYPE>(r.high[i]);
            if (thislow != rlow || thishigh != rhigh)
                return false;
        }
    }

    return true;
}
double varRegion::getMinimumDistance(const varRegion& r) const {
    if (dimension != r.dimension) assert(false);
    auto ret = 0.0;

    for (auto i = 0u; i < dimension; i++) {
        auto x = 0.0;
        if (match(low[i], [](REAL_TYPE) {return true;}, [](DISCRETE_TYPE) {return false;})) {
            auto thislow = boost::get<REAL_TYPE>(low[i]);
            auto thishigh = boost::get<REAL_TYPE>(high[i]);
            auto rlow = boost::get<REAL_TYPE>(r.low[i]);
            auto rhigh = boost::get<REAL_TYPE>(r.high[i]);

            if (rhigh < thislow)
                x = std::abs(rhigh - thislow);
            else if (thishigh < rlow)
                x = std::abs(thishigh - rlow);
        } else {
            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
            auto rlow = boost::get<DISCRETE_TYPE>(r.low[i]);
            auto rhigh = boost::get<DISCRETE_TYPE>(r.high[i]);

            if (rhigh < thislow)
                x = std::abs(rhigh - thislow);
            else if (thishigh < rlow)
                x = std::abs(thishigh - rlow);
        }
        ret += x*x;
    }
    return std::sqrt(ret);
}
bool varRegion::containsPoint(const varPoint& p) const {
    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i],
                    [&](REAL_TYPE){
                        auto thislow = boost::get<REAL_TYPE>(low[i]);
                        auto thishigh = boost::get<REAL_TYPE>(high[i]);
                        if (thislow > boost::get<REAL_TYPE>(p.coords[i]) || thishigh < boost::get<REAL_TYPE>(p.coords[i])) return true;
                        return false;
                    },
                    [&](DISCRETE_TYPE){
                        auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
                        auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
                        if (thislow > boost::get<DISCRETE_TYPE>(p.coords[i]) || thishigh < boost::get<DISCRETE_TYPE>(p.coords[i])) return true;
                        return false;

                    })) {
            return false;
        }
    }
    return true;
}
bool varRegion::touchesPoint(const varPoint& p) const {
    for (auto i = 0u; i < dimension; i++) {
         if (match(low[i],
                    [&](REAL_TYPE){
                        auto thislow = boost::get<REAL_TYPE>(low[i]);
                        auto thishigh = boost::get<REAL_TYPE>(high[i]);
                        if (std::fabs(thislow - boost::get<REAL_TYPE>(p.coords[i])) <= std::numeric_limits<double>::epsilon() || std::fabs(thishigh - boost::get<REAL_TYPE>(p.coords[i])) <= std::numeric_limits<double>::epsilon()) return true;
                        return false;
                    },
                    [&](DISCRETE_TYPE){
                        auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
                        auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
                        if (thislow == boost::get<DISCRETE_TYPE>(p.coords[i]) || thishigh == boost::get<DISCRETE_TYPE>(p.coords[i])) return true;
                        return false;

                    })) {
            return true;
        }
    }
    return false;
}
double varRegion::getMinimumDistance(const varPoint& p) const {
    auto ret = 0.0;
    for (auto i = 0u; i < dimension; i++) {
        ret += match(low[i], [&](REAL_TYPE l) {
            auto h = boost::get<REAL_TYPE>(high[i]);
            auto pi = boost::get<REAL_TYPE>(p.coords[i]);
            if (pi < l) return std::pow(l - pi, 2.0);
            else if (pi > h) return std::pow(pi - h, 2.0);
            return 0.0;
        }, [&](DISCRETE_TYPE l) {
            auto h = boost::get<DISCRETE_TYPE>(high[i]);
            auto pi = boost::get<DISCRETE_TYPE>(p.coords[i]);
            if (pi < l) return std::pow(l - pi, 2.0);
            else if (pi > h) return std::pow(pi - h, 2.0);
            return 0.0;
        });
    }
    return std::sqrt(ret);
}
varRegion varRegion::getIntersectingRegion(const varRegion& r) const {
    varRegion ret; ret.makeInfinite(dimension);
    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [&](REAL_TYPE) {
                return boost::get<REAL_TYPE>(low[i]) > boost::get<REAL_TYPE>(r.high[i]) || boost::get<REAL_TYPE>(high[i]) < boost::get<REAL_TYPE>(r.low[i]);
            }, [&](DISCRETE_TYPE) {
                return boost::get<DISCRETE_TYPE>(low[i]) > boost::get<DISCRETE_TYPE>(r.high[i]) || boost::get<DISCRETE_TYPE>(high[i]) < boost::get<DISCRETE_TYPE>(r.low[i]);
            })) {
            return ret;
        }
    }

    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
            ret.low[i] = std::max(boost::get<REAL_TYPE>(low[i]), boost::get<REAL_TYPE>(r.low[i]));
            ret.high[i] = std::max(boost::get<REAL_TYPE>(high[i]), boost::get<REAL_TYPE>(r.high[i]));
        }, [&](DISCRETE_TYPE) {
            ret.low[i] = std::max(boost::get<DISCRETE_TYPE>(low[i]), boost::get<DISCRETE_TYPE>(r.low[i]));
            ret.high[i] = std::max(boost::get<DISCRETE_TYPE>(high[i]), boost::get<DISCRETE_TYPE>(r.high[i]));
        });
    }
    return ret;
}
double varRegion::getIntersectingArea(const varRegion& r) const {
    auto ret = 1.0;
    double f1, f2;

    for (auto i = 0u; i < dimension; i++) {
        if (match(low[i], [&](REAL_TYPE) {
            auto l = boost::get<REAL_TYPE>(low[i]);
            auto h = boost::get<REAL_TYPE>(high[i]);
            auto rl = boost::get<REAL_TYPE>(r.low[i]);
            auto rh = boost::get<REAL_TYPE>(r.high[i]);
            if (l > rh || h < rl) return true;
            f1 = std::max(l, rl);
            f2 = std::min(h, rh);
            return false;
        }, [&](DISCRETE_TYPE) {
            auto l = boost::get<DISCRETE_TYPE>(low[i]);
            auto h = boost::get<DISCRETE_TYPE>(high[i]);
            auto rl = boost::get<DISCRETE_TYPE>(r.low[i]);
            auto rh = boost::get<DISCRETE_TYPE>(r.high[i]);
            if (l > rh || h < rl) return true;
            f1 = std::max(l, rl);
            f2 = std::min(h, rh);
            return false;
        })) {
            return 0.0;
        }
        ret *= (f2 - f1);
    }
    return ret;
}
double varRegion::getMargin() const {
    auto mul = std::pow(2.0, static_cast<double>(dimension) - 1.0);
    auto margin = 0.0;
    for (auto i = 0u; i < dimension; i++) {
        margin += match(low[i], [&](REAL_TYPE) {
                return boost::get<REAL_TYPE>(high[i]) - boost::get<REAL_TYPE>(low[i]);
            }, [&](DISCRETE_TYPE) {
                return (REAL_TYPE)(boost::get<DISCRETE_TYPE>(high[i]) - boost::get<DISCRETE_TYPE>(low[i]));
            }) * mul;
    }
    return margin;
}
void varRegion::combineRegion(const varRegion& r) {
    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
                low[i] = std::min(boost::get<REAL_TYPE>(low[i]), boost::get<REAL_TYPE>(r.low[i]));
                high[i] = std::min(boost::get<REAL_TYPE>(high[i]), boost::get<REAL_TYPE>(r.high[i]));
            }, [&](DISCRETE_TYPE) {
                low[i] = std::min(boost::get<DISCRETE_TYPE>(low[i]), boost::get<DISCRETE_TYPE>(r.low[i]));
                high[i] = std::min(boost::get<DISCRETE_TYPE>(high[i]), boost::get<DISCRETE_TYPE>(r.high[i]));
            });
    }
}
void varRegion::combinePoint(const varPoint& p) {
    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
                low[i] = std::min(boost::get<REAL_TYPE>(low[i]), boost::get<REAL_TYPE>(p.coords[i]));
                high[i] = std::min(boost::get<REAL_TYPE>(high[i]), boost::get<REAL_TYPE>(p.coords[i]));
            }, [&](DISCRETE_TYPE) {
                low[i] = std::min(boost::get<DISCRETE_TYPE>(low[i]), boost::get<DISCRETE_TYPE>(p.coords[i]));
                high[i] = std::min(boost::get<DISCRETE_TYPE>(high[i]), boost::get<DISCRETE_TYPE>(p.coords[i]));
            });
    }
}
void varRegion::getCombinedRegion(varRegion& out, const varRegion& in) const {
    out = *this;
    out.combineRegion(in);
}
/*coord_type varRegion::getLow(uint32_t index) const {
    return low[index];
}
coord_type varRegion::getHigh(uint32_t index) const {
    return high[index];
}*/
void varRegion::makeInfinite(uint32_t dimension) {
    makeDimension(dimension);
    for (auto i = 0u; i < dimension; i++) {
        match(low[i], [&](REAL_TYPE) {
                low[i] = std::numeric_limits<REAL_TYPE>::max();
                high[i] = -std::numeric_limits<REAL_TYPE>::max();
            },
            [&](DISCRETE_TYPE) {
                low[i] = std::numeric_limits<DISCRETE_TYPE>::max();
                high[i] = -std::numeric_limits<DISCRETE_TYPE>::max();
            });
    }
}
void varRegion::makeDimension(uint32_t dimension) {
    if (this->dimension != dimension) {
        std::vector<bool> continuous(dimension);
        std::transform(low, low+dimension, continuous.begin(),
                [](const coord_type& p) { return p.type() == typeid(REAL_TYPE); });
        delete[] low;
        delete[] high;
        low = high = nullptr;
        this->dimension = dimension;
        low = new coord_type[dimension];
        high = new coord_type[dimension];
        // TODO use continuous maybe
    }
}

