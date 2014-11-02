#include "common.h"
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include <iostream>
#include <random>
#include "point.h"

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

class varPoint;
class varRegion;
#include "varPoint.h"
#include "varRegion.h"
#include "particle.h"

SpatialIndex::ISpatialIndex* tree;
SpatialIndex::id_type id;

REAL_TYPE extractReal(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return p(); }, [](const Point<DISCRETE_TYPE>& p) { return (REAL_TYPE)p(); });
}

DISCRETE_TYPE extractDiscrete(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return (DISCRETE_TYPE)p(); }, [](const Point<DISCRETE_TYPE>& p) { return p(); });
}

class CheckerVisitor : public SpatialIndex::IVisitor {
    bool found_;
    double data_;
public:
    CheckerVisitor() : SpatialIndex::IVisitor(), found_(false), data_(0.0) {}
    void visitNode(const SpatialIndex::INode&) {}
    void visitData(const SpatialIndex::IData& d) {
        byte* data; uint32_t len; d.getData(len, &data);
        found_ = true;
        data_ = *(reinterpret_cast<REAL_TYPE*>(data));
        delete[] data;
    }
    void visitData(std::vector<const SpatialIndex::IData*>&) {}
    bool found() const {
        return found_;
    }
    double data() const {
        return data_;
    }
};

unsigned long long function_calls;
unsigned long long total_function_calls;

double cube_width = 0.05;

REAL_TYPE rosenbrock(const points_vector& x) {
    ++total_function_calls;
    std::vector<double> low(x.size());
    std::vector<double> high(x.size());
    std::vector<double> casted(x.size());
    for (auto i = 0u; i < x.size(); i++) {
        low[i] = match(x[i], [&](Point<REAL_TYPE> p) {
                return std::max(p.left, p() - cube_width);
            }, [&] (Point<DISCRETE_TYPE> p) {
                return std::max((double) p.left, (double)p() - cube_width);
            });
        high[i] = match(x[i], [&](Point<REAL_TYPE> p) {
                return std::min(p.right, p() + cube_width);
            }, [&] (Point<DISCRETE_TYPE> p) {
                return std::min((double) p.right, (double)p() + cube_width);
            });
        casted[i] = extractReal(x[i]);
    }
    SpatialIndex::Region point_region(&low[0], &high[0], x.size());
    SpatialIndex::Point point(&casted[0], x.size());
    CheckerVisitor visitor;
    tree->pointLocationQuery(point, visitor);
    if (visitor.found()) {
        return visitor.data();
    }

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
    ++function_calls;
    return ret_val;
}

void test() {
    coord_type low[2]; low[0] = -1ll, low[1] = -1ll;
    coord_type high[2]; high[0] = 1ll, high[1] = 0ll;
    coord_type higher[2]; higher[0] = 1ll, higher[1] = 1ll;
    assert(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(higher, 2u)));
    assert(!(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(high, 2u))));
    points_vector x = {Point<REAL_TYPE>(1.2, -1.0, 2.0), Point<REAL_TYPE>(1.0, -1.0, 2.0)};
    std::cout << "This point = " << rosenbrock(x) << std::endl;
    x = {Point<REAL_TYPE>(1.0, -1.0, 2.0), Point<REAL_TYPE>(1.0, -1.0, 2.0)};
    std::cout << "new point = " << rosenbrock(x) << std::endl;
}

std::tuple<points_vector, REAL_TYPE> mesh_search(std::function<REAL_TYPE(const points_vector&)> f, const points_vector& x) {
    std::vector<bool> continuous(x.size());
    std::transform(x.begin(), x.end(), continuous.begin(),
            [](const point_type& p) { return p.type() == typeid(Point<double>); });
    for (auto i : continuous)
        std::cerr << (i ? "Continuous" : "Discrete") << std::endl;
    std::vector<coord_type> left(x.size());
    // The mesh right window
    std::vector<coord_type> right(x.size());
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
    auto meshWidth = 1.0;
    while (improvement || constrictions < 20) {
        improvement ? (std::cerr << "Had an improvement" << std::endl) : (std::cerr << "Going around again" << std::endl);
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
                auto v = std::min((double)nodesPerRow, std::min(meshWidth + 1, (double)(p.right - p.left + 1)));
                numPoints *= v;
            }
        }

        // there are min(npr, width, right - left + 1) nodes in a continuous row.
        // This has been rounded down to an odd number

        for (auto i = 0u; i < test.size(); i++) {
            if (continuous[i]) {
                auto& p = boost::get<Point<REAL_TYPE>>(test[i]);
                auto left = (meshWidth < (p.right - p.left)) ? p() - meshWidth / 2 : p.left;
                auto right = (meshWidth < (p.right - p.left)) ? p() + meshWidth / 2 : p.right;
                if (i == 0) {
                    test[i] = std::move(Point<REAL_TYPE>(left - (right - left)/nodesPerRow, p.left, p.right));
                } else {
                    test[i] = std::move(Point<REAL_TYPE>(left, p.left, p.right));
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
    std::string basename = "base.rtree";
    SpatialIndex::IStorageManager* diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(basename, 4096);
    SpatialIndex::StorageManager::IBuffer* file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    SpatialIndex::id_type indexIdentifier;
    tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 100, 100, /*points.size()*/2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

    std::vector<points_vector> start_points = {{Point<double>(0.0, -1.0, 2.0), Point<long long>(0, -1, 2)},
                                                {Point<double>(0.1, -1.0, 2.0), Point<long long>(-1, -1, 2)}};
    points_vector output;
    REAL_TYPE f;

    std::tie(output, f) = particle_swarm(rosenbrock, start_points);
    std::cout << "Seeding mesh_search with (";
    for (auto p : output) {
        std::cout << extractReal(p) << " ";
    }
    std::cout << "\b) for a value of " << f << std::endl;
    //output = start_points[0];

    std::tie(output, f) = mesh_search(rosenbrock, output);
    std::cout << "Best f was " << f << std::endl;
    for (auto p : output)
        std::cout << extractReal(p) << " ";
    std::cout << std::endl;

    std::cout << "A total of " << function_calls << " calls out of " << total_function_calls << std::endl;

    //test();

    return 0;
}
