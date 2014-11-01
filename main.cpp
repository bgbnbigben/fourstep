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

class varPoint;
class varRegion;
#include "varPoint.h"
#include "varRegion.h"

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
