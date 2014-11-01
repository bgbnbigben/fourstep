#include <boost/variant.hpp>
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include "inline_variant.hpp"
#include <iostream>
#include <cstring>

#define DISCRETE_TYPE long long
#define REAL_TYPE double

template <class T>
class Point { 
    public:
        const T val;
        const T left;
        const T right;

        Point() = delete;
        Point(const T v, const T l, const T r) : val(v), left(l), right(r) {};
        Point(const Point& p) : val(p.val), left(p.left), right(p.right) {};
        Point(const Point&& p) : val(p.val), left(p.left), right(p.right) {};

        Point<T>& operator=(const Point<T>& p) {
            // Urgh.
            auto t = const_cast<T*>(&this->val);
            *t = p.val;
            t = const_cast<T*>(&this->left);
            *t = p.left;
            t = const_cast<T*>(&this->right);
            *t = p.right;

            return *this;
        }

        Point<T>& operator=(const Point<T>&& p) {
            auto t = const_cast<T*>(&this->val);
            *t = p.val;
            t = const_cast<T*>(&this->left);
            *t = p.left;
            t = const_cast<T*>(&this->right);
            *t = p.right;

            return *this;
        }

        T operator()() const { return val; }
};

typedef boost::variant<Point<DISCRETE_TYPE>, Point<REAL_TYPE>> point_type;
typedef std::vector<point_type> points_vector;
typedef boost::variant<DISCRETE_TYPE, REAL_TYPE> coord_type;

class varPoint : /*public SpatialIndex::IObject, public virtual SpatialIndex::IShape*/ public SpatialIndex::Point {
public:
    uint32_t dimension;
    coord_type* coords;
    friend class varRegion;
    varPoint() : Point(), dimension(0), coords(0) {}
    varPoint(const coord_type* coords, uint32_t dimension) : Point(), dimension(dimension) {
        this->coords = new coord_type[dimension];
        memcpy(this->coords, coords, dimension*sizeof(coord_type));
    }
    varPoint(const varPoint& p) : Point(), dimension(p.dimension) {
        coords = new coord_type[dimension];
        memcpy(coords, p.coords, dimension*sizeof(coord_type));
    }
    ~varPoint() {
        delete [] coords;
    }
    varPoint& operator=(const varPoint& p) {
        makeDimension(p.dimension);
        memcpy(coords, p.coords, dimension*sizeof(coord_type));
        return *this;
    }
    bool operator==(const varPoint& p) const {
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
    varPoint* clone() {
        return new varPoint(*this);
    }
    uint32_t getByteArraySize() {
        return sizeof(uint32_t) + dimension*sizeof(coord_type);
    }
    void loadFromByteArray(const byte* ptr) {
        uint32_t dim; memcpy(&dim, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);
        makeDimension(dim);
        memcpy(coords, ptr, dimension*sizeof(coord_type));
    }
    void storeToByteArray(byte** data, uint32_t& len) {
        len = getByteArraySize();
        *data = new byte[len];
        byte* ptr = *data;
        memcpy(ptr, &dimension, sizeof(uint32_t));
        ptr += sizeof(uint32_t);
        memcpy(ptr, coords, dimension*sizeof(coord_type));
    }
    bool intersectsShape(const SpatialIndex::IShape& s) const {
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return pr->containsPoint(*this);
        assert(0);
    }
    bool containsShape(const SpatialIndex::IShape&) const {
        return false;
    }
    bool touchesShape(const SpatialIndex::IShape& s) const {
        const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
        if (ppt != 0) {
            if (*this == *ppt) return true;
            return false;
        }
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return pr->touchesPoint(*this);
        assert(0);
    }
    void getCenter(varPoint& out) const {
        out = *this;
    }
    uint32_t getDimension() const {
        return dimension;
    }
    void getMBR(varRegion& out) const {
        out = varRegion(coords, coords, dimension);
    }
    double getArea() const {
        return 0.0;
    }
    double getMinimumDistance(const SpatialIndex::IShape& in) const {
        const varPoint* ppt = dynamic_cast<const varPoint*>(&in);
        if (ppt != 0) return getMinimumDistance(*ppt);
        const varRegion* pr = dynamic_cast<const varRegion*>(&in);
        if (pr != 0) return pr->getMinimumDistance(*this);
        assert(0);
    }
    double getMinimumDistance(const varPoint& p) const {
        auto ret = 0.0;
        for (auto i = 0u; i < dimension; i++) {
            ret += std::pow(match(coords[i], [&](REAL_TYPE px) {return px - boost::get<REAL_TYPE>(p.coords[i]);}, [&](DISCRETE_TYPE px) {return (REAL_TYPE)(px - boost::get<DISCRETE_TYPE>(p.coords[i]));}), 2.0);
        }
        return std::sqrt(ret);
    }
    coord_type getCoordinate(uint32_t index) const {
        return coords[index];
    }
    void makeInfinite(uint32_t dimension) {
        makeDimension(dimension);
        for (auto i = 0u; i < dimension; i++) {
            coords[i] = std::numeric_limits<DISCRETE_TYPE>::max();
        }
    }
    void makeDimension(uint32_t dim) {
        if (dimension != dim) {
            delete[] coords;
            dimension = dim;
            coords = new coord_type[dimension];
        }
    }
};

class varRegion : /*public SpatialIndex::IObject, public virtual SpatialIndex::IShape*/public SpatialIndex::Region {
public:
    coord_type* low;
    coord_type* high;
    uint32_t dimension;

    varRegion() : SpatialIndex::Region(), low(nullptr), high(nullptr), dimension(0) {}
    varRegion(const coord_type* low, const coord_type* high, uint32_t dimension) : SpatialIndex::Region() {
        init(low, high, dimension);
    }
    varRegion(const varPoint& low, const varPoint& high) : SpatialIndex::Region() {
        init(low.coords, high.coords, low.dimension);
    }
    void init(const coord_type* low, const coord_type* high, uint32_t dimension) {
        low = new coord_type[dimension];
        high = new coord_type[dimension];
        memcpy(this->low, low, dimension*sizeof(coord_type));
        memcpy(this->high, high, dimension*sizeof(coord_type));
    }
    ~varRegion() {
        delete[] low;
        delete[] high;
    }
    varRegion& operator=(const varRegion& r) {
        makeDimension(r.dimension);
        memcpy(low, r.low, dimension*sizeof(coord_type));
        memcpy(high, r.high, dimension*sizeof(coord_type));
        return *this;
    }
    bool operator==(const varRegion& r) {
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
    varRegion* clone() {
        return new varRegion(*this);
    }
    uint32_t getByteArraySize() {
        return sizeof(uint32_t) + 2*dimension*sizeof(coord_type);
    }
    void loadFromByteArray(const byte* ptr) {
        uint32_t dim;
        memcpy(&dim, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);
        makeDimension(dim);
        memcpy(low, ptr, dimension*sizeof(coord_type));
        ptr += dimension*sizeof(coord_type);
        memcpy(high, ptr, dimension*sizeof(coord_type));
        ptr += dimension*sizeof(coord_type);
    }
    void storeToByteArray(byte** data, uint32_t& len) {
        len = getByteArraySize();
        *data = new byte[len];
        byte* ptr = *data;
        memcpy(ptr, &dimension, sizeof(uint32_t));
        ptr += sizeof(uint32_t);
        memcpy(ptr, low, dimension*sizeof(coord_type));
        ptr += dimension*sizeof(coord_type);
        memcpy(ptr, high, dimension*sizeof(coord_type));
    }
    bool intersectsShape(const SpatialIndex::IShape& s) const {
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return intersectsRegion(*pr);
        const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
        if (ppt != 0) return containsPoint(*ppt);
        assert(0);
    }
    bool containsShape(const SpatialIndex::IShape &s) const {
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return containsRegion(*pr);
        const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
        if (ppt != 0) return containsPoint(*ppt);
    }
    bool touchesShape(const SpatialIndex::IShape& s) const {
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return touchesRegion(*pr);
        const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
        if (ppt != 0) return touchesPoint(*ppt);
    }
    void getCenter(varPoint& out) const {
        out.makeDimension(dimension);
        for (auto i = 0u; i < dimension; i++) {
            match(low[i], [&](REAL_TYPE) {
                out.coords[i] = (boost::get<REAL_TYPE>(low[i]) + boost::get<REAL_TYPE>(high[i])) / 2.0;
            }, [&](DISCRETE_TYPE) {
                out.coords[i] = (boost::get<DISCRETE_TYPE>(low[i]) + boost::get<DISCRETE_TYPE>(high[i])) / 2.0;
            });
        }
    }
    uint32_t getDimension() {
        return dimension;
    }
    void getMBR(varRegion& out) const {
        out = *this;
    }
    double getArea() const {
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
    double getMinimumDistance(const SpatialIndex::IShape& s) const {
        const varRegion* pr = dynamic_cast<const varRegion*>(&s);
        if (pr != 0) return getMinimumDistance(*pr);
        const varPoint* ppt = dynamic_cast<const varPoint*>(&s);
        if (ppt != 0) return getMinimumDistance(*ppt);
    }
    bool intersectsRegion(const varRegion& r) const {
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
    bool containsRegion(const varRegion& r) const {
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
    bool touchesRegion(const varRegion& r) const {
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
    double getMinimumDistance(const varRegion& r) const {
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
    bool containsPoint(const varPoint& p) const {
        for (auto i = 0u; i < dimension; i++) {
            if (match(low[i],
                        [&](REAL_TYPE){
                            auto thislow = boost::get<REAL_TYPE>(low[i]);
                            auto thishigh = boost::get<REAL_TYPE>(high[i]);
                            if (thislow > boost::get<REAL_TYPE>(p.getCoordinate(i)) || thishigh < boost::get<REAL_TYPE>(p.getCoordinate(i))) return true;
                            return false;
                        },
                        [&](DISCRETE_TYPE){
                            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
                            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
                            if (thislow > boost::get<DISCRETE_TYPE>(p.getCoordinate(i)) || thishigh < boost::get<DISCRETE_TYPE>(p.getCoordinate(i))) return true;
                            return false;

                        })) {
                return false;
            }
        }
        return true;
    }
    bool touchesPoint(const varPoint& p) const {
        for (auto i = 0u; i < dimension; i++) {
             if (match(low[i],
                        [&](REAL_TYPE){
                            auto thislow = boost::get<REAL_TYPE>(low[i]);
                            auto thishigh = boost::get<REAL_TYPE>(high[i]);
                            if (std::fabs(thislow - boost::get<REAL_TYPE>(p.getCoordinate(i))) <= std::numeric_limits<double>::epsilon() || std::fabs(thishigh - boost::get<REAL_TYPE>(p.getCoordinate(i))) <= std::numeric_limits<double>::epsilon()) return true;
                            return false;
                        },
                        [&](DISCRETE_TYPE){
                            auto thislow = boost::get<DISCRETE_TYPE>(low[i]);
                            auto thishigh = boost::get<DISCRETE_TYPE>(high[i]);
                            if (thislow == boost::get<DISCRETE_TYPE>(p.getCoordinate(i)) || thishigh == boost::get<DISCRETE_TYPE>(p.getCoordinate(i))) return true;
                            return false;

                        })) {
                return true;
            }
        }
        return false;
    }
    double getMinimumDistance(const varPoint& p) const {
        auto ret = 0.0;
        for (auto i = 0u; i < dimension; i++) {
            ret += match(low[i], [&](REAL_TYPE l) {
                auto h = boost::get<REAL_TYPE>(high[i]);
                auto pi = boost::get<REAL_TYPE>(p.getCoordinate(i));
                if (pi < l) return std::pow(l - pi, 2.0);
                else if (pi > h) return std::pow(pi - h, 2.0);
            }, [&](DISCRETE_TYPE l) {
                auto h = boost::get<DISCRETE_TYPE>(high[i]);
                auto pi = boost::get<DISCRETE_TYPE>(p.getCoordinate(i));
                if (pi < l) return std::pow(l - pi, 2.0);
                else if (pi > h) return std::pow(pi - h, 2.0);
            });
        }
        return std::sqrt(ret);
    }
    varRegion getIntersectingRegion(const varRegion& r) const {
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
    }
    double getIntersectingArea(const varRegion& r) const {
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
    double getMargin() const {
        double mul = std::pow(2.0, static_cast<double>(dimension) - 1.0);
        double margin = 0.0;
        for (auto i = 0u; i < dimension; i++) {
            margin += match(low[i], [&](REAL_TYPE) {
                    return boost::get<REAL_TYPE>(high[i]) - boost::get<REAL_TYPE>(low[i]);
                }, [&](DISCRETE_TYPE) {
                    return (REAL_TYPE)(boost::get<DISCRETE_TYPE>(high[i]) - boost::get<DISCRETE_TYPE>(low[i]));
                }) * mul;
        return margin;
        }
    }
    void combineRegion(const varRegion& r) {
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
    void combinePoint(const varPoint& p) {
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
    void getCombinedRegion(varRegion& out, const varRegion& in) const {
        out = *this;
        out.combineRegion(in);
    }
    /*coord_type getLow(uint32_t index) const {
        return low[index];
    }
    coord_type getHigh(uint32_t index) const {
        return high[index];
    }*/
    void makeInfinite(uint32_t dimension) {
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
    void makeDimension(uint32_t dimension) {
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
    
};

void test() {
    assert(varRegion({-1, -1}, {1, 1}) == varRegion({-1, -1}, {1, 1}));
    assert(varRegion({-1, -1}, {1, 1}) == varRegion({-1, -1}, {1, 0}));
}

SpatialIndex::ISpatialIndex* tree;
SpatialIndex::id_type id;

REAL_TYPE extractReal(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return p(); }, [](const Point<DISCRETE_TYPE>& p) { return (REAL_TYPE)p(); });
}

DISCRETE_TYPE extractDiscrete(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return (DISCRETE_TYPE)p(); }, [](const Point<DISCRETE_TYPE>& p) { return p(); });
}

REAL_TYPE rosenbrock(const points_vector& x) {
    tree->intersectsWithQuery(point_region, visitor);
    if (visitor.found())
        return visitor.data();

    auto start = extractReal(x[0]);
    auto ret_val = (1.0 - start)*(1 - start);
    std::cerr << "Testing ("; 
    for (auto i = x.size() - 1; i > 0; i--) {
        std::cerr << extractReal(x[i]) << ", ";
        auto curr = extractReal(x[i]);
        auto next = extractReal(x[i-1]);
        ret_val += 100*(curr - next*next)*(curr - next*next);
    }
    std::cerr << extractReal(x[0]) << ")" << std::endl;

    tree->insertData(sizeof(REAL_TYPE), reinterpret_cast<const byte*>(&ret_val), point_region, id++);
    return ret_val;
}

// f [in], x [in]
std::tuple<points_vector, REAL_TYPE> mesh_search(std::function<REAL_TYPE(const points_vector&)> f, const points_vector x) {
    std::vector<bool> continuous(x.size());
    std::transform(x.begin(), x.end(), continuous.begin(),
            [](const point_type& p) { return p.type() == typeid(Point<double>); });
    for (auto i : continuous)
        std::cerr << (i ? "Continuous" : "Discrete") << std::endl;
    std::vector<boost::variant<DISCRETE_TYPE, REAL_TYPE>> left(x.size());
    // The mesh right window
    std::vector<boost::variant<DISCRETE_TYPE, REAL_TYPE>> right(x.size());
    for (auto i = 0u; i < x.size(); i++) {
        if (continuous[i]) {
            left[i] = boost::get<Point<REAL_TYPE>>(x[i]).left;
            right[i] = boost::get<Point<REAL_TYPE>>(x[i]).right;
        } else {
            left[i] = boost::get<Point<DISCRETE_TYPE>>(x[i]).left;
            right[i] = boost::get<Point<DISCRETE_TYPE>>(x[i]).right;
        }
    }

    REAL_TYPE bestF = std::numeric_limits<REAL_TYPE>::max();
    points_vector bestX(x);

    bool improvement = true;
    int constrictions = 0; // take 2^constrictions number of nodes per line.
    while (improvement || constrictions < 10) {
        std::cerr << "Had an improvement" << std::endl;
        improvement = false;
        points_vector test(bestX);

        long long numPoints = 1;
        long long nodesPerRow = pow(2, constrictions) + 1;
        std::cerr << "There are up to " << nodesPerRow << " npr" << std::endl;
        for (auto i = 0u; i < test.size(); i++) {
            if (continuous[i]) {
                numPoints *= nodesPerRow;
            } else {
                auto& p = boost::get<Point<DISCRETE_TYPE>>(bestX[i]);
                auto v = std::min(nodesPerRow, p.right - p.left + 1);
                numPoints *= v;
            }
        }

        // there are std::min(npr, right - left + 1) nodes in a continuous row.
        // This has been rounded down to an odd number

        for (auto i = 0u; i < test.size(); i++) {
            if (continuous[i]) {
                auto& p = boost::get<Point<REAL_TYPE>>(test[i]);
                if (i == 0) {
                    test[i] = std::move(Point<REAL_TYPE>(p.left - (p.right - p.left)/nodesPerRow, p.left, p.right));
                } else {
                    test[i] = std::move(Point<REAL_TYPE>(p.left, p.left, p.right));
                }
            } else {
                auto& p = boost::get<Point<DISCRETE_TYPE>>(test[i]);
                const auto best = boost::get<Point<DISCRETE_TYPE>>(bestX[i]);
                auto inRow = std::min(nodesPerRow, best.right - best.left + 1);
                auto left = std::min(inRow / 2, best.val - best.left);
                if (i == 0) {
                    test[i] = std::move(Point<DISCRETE_TYPE>(best.val - left - 1, p.left, p.right));
                } else {
                    test[i] = std::move(Point<DISCRETE_TYPE>(best.val - left, p.left, p.right));
                }
            }
        }

        std::cerr << "We have " << numPoints << " points on constriction " << constrictions << std::endl;
        for (auto pointNum = 0; pointNum < numPoints; pointNum++) {
            int i = 0;
            bool tick = false;
            do {
                if (continuous[i]) {
                    auto& current = boost::get<Point<REAL_TYPE>>(test[i]);
                    test[i] = std::move(Point<REAL_TYPE>(current.val + (current.right - current.left) / nodesPerRow, current.left, current.right));
                    if (current.val > current.right) {
                        test[i] = std::move(Point<REAL_TYPE>(current.left, current.left, current.right));
                        i++;
                        tick = true;
                    } else {
                        tick = false;
                    }
                } else {
                    auto& current = boost::get<Point<DISCRETE_TYPE>>(test[i]);
                    const auto best = boost::get<Point<DISCRETE_TYPE>>(bestX[i]);
                    auto inRow = std::min(nodesPerRow, best.right - best.left + 1);
                    auto left = std::min(inRow / 2, best.val - best.left);
                    auto right = inRow - left + 1;
                    test[i] = std::move(Point<DISCRETE_TYPE>(current.val + 1, current.left, current.right));
                    if (current.val > best.val + right) {
                        test[i] = std::move(Point<DISCRETE_TYPE>(best.val - left, current.left, current.right));
                        i++;
                        tick = true;
                    } else {
                        tick = false;
                    }
                }
            } while (tick);
            if (f(test) < bestF) {
                std::cerr << "Best found" << std::endl;
                bestF = f(test);
                bestX = test;
                improvement = true;
                constrictions = 0;
            }
        }
        constrictions++;
    }
    std::cerr << bestF << std::endl;
    return std::make_tuple(bestX, bestF);
}


int main() {
    points_vector points = {Point<double>(0.0, -1.0, 2.0), Point<long long>(0, -1, 2)};
    points_vector output;
    REAL_TYPE f;
    std::tie(output, f) = mesh_search(rosenbrock, {Point<double>(0.0, -1.0, 2.0), Point<double>(0.0, -1.0, 2.0)});
    std::cout << "Best f was " << f << std::endl;
    for (auto p : output)
        std::cout << extractReal(p) << " ";
    std::cout << std::endl;

    std::string basename = "base.rtree";
    SpatialIndex::IStorageManager* diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(basename, 4096);
    SpatialIndex::StorageManager::IBuffer* file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    SpatialIndex::id_type indexIdentifier;
    tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 100, 100, points.size(), SpatialIndex::RTree::RV_RSTAR, indexIdentifier);


    std::tie(output, f) = mesh_search(rosenbrock, points);
    std::cout << "Best f was " << f << std::endl;
    for (auto p : output)
        std::cout << extractReal(p) << " ";
    std::cout << std::endl;

    return 0;
}
