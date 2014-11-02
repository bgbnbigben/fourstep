#include "common.h"
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include <iostream>
#include <random>
#include "point.h"
#include "mesh_search.h"


//class varPoint;
//class varRegion;
//#include "varPoint.h"
//#include "varRegion.h"
#include "particle.h"

SpatialIndex::ISpatialIndex* tree;
SpatialIndex::id_type id;

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
    //assert(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(higher, 2u)));
    //assert(!(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(high, 2u))));
    points_vector x = {Point<REAL_TYPE>(1.2, -1.0, 2.0), Point<REAL_TYPE>(1.0, -1.0, 2.0)};
    std::cout << "This point = " << rosenbrock(x) << std::endl;
    x = {Point<REAL_TYPE>(1.0, -1.0, 2.0), Point<REAL_TYPE>(1.0, -1.0, 2.0)};
    std::cout << "new point = " << rosenbrock(x) << std::endl;
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
