#include "common.h"
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include <iostream>
#include <random>
#include "point.h"
#include "mesh_search.h"
//#include <mpi.h>


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


#define DIE_TAG 0xd1ed1e


int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);

    MPI::Status stat;
    auto numprocs = MPI::COMM_WORLD.Get_size();
    auto rank = MPI::COMM_WORLD.Get_rank();
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    std::vector<points_vector> start_points = {{Point<double>(0.0, -1.0, 2.0), Point<long long>(0, -1, 2)},
                                                {Point<double>(0.1, -1.0, 2.0), Point<long long>(-1, -1, 2)}};
    auto numDimensions = start_points[0].size();
    // Fuck it, bin process 0 as only a data marshall. Whatever.
    if (rank == 0) {
        auto globalBest = std::numeric_limits<double>::max();
        points_vector globalData(numDimensions);
        auto partitions = 100;
        auto generatePartitions = [&](auto numPartitions) {
            std::vector<points_vector> ret(partitions);
            std::vector<int> splits(numDimensions, -1);
            auto maxPartitions = std::lround(std::pow(numDimensions, 1./partitions));
            auto numPartitionsLeft = partitions;
            auto continuous = 0u;
            for (auto i = 0u; i < numDimensions; i++) {
                match(start_points[0][i], [&](Point<REAL_TYPE>) {
                        continuous++;
                    }, [&](Point<DISCRETE_TYPE> p) {
                        splits[i] = std::min((DISCRETE_TYPE)maxPartitions, p.right - p.left + 1);
                        numPartitionsLeft = std::lround(numPartitionsLeft / (double)splits[i]);
                    });
            }
            for (auto& split: splits) {
                if (split == -1)
                    split = std::lround(std::pow(numPartitionsLeft, 1./continuous));
            }

            auto actualPartitions = 1u;
            for (auto split: splits)
                actualPartitions *= split;

            std::cout << "Getting " << actualPartitions << " when we wanted " << numPartitions << std::endl;
            ret[0] = start_points[0];
            for (auto i = 1u; i < actualPartitions; i++) {
                ret[i] = ret[i-1];
                auto tick = false;
                auto j = 0;
                do {
                    match(ret[i][j], [&](Point<REAL_TYPE> p) {
                        REAL_TYPE left;
                        auto start_point = boost::get<Point<REAL_TYPE>>(start_points[0][j]);
                        auto jump = (start_point.right - start_point.left) / splits[j];
                        if (std::fabs(p.right - start_point.right) <= std::numeric_limits<REAL_TYPE>::epsilon()) {
                            tick = true;
                            left = p.left;
                        } else {
                            left = p.right;
                        }
                        ret[i][j] = Point<REAL_TYPE>(left + jump/2, left, left + jump);
                    }, [&](Point<DISCRETE_TYPE> p) {
                        DISCRETE_TYPE left;
                        if (p.right == boost::get<Point<DISCRETE_TYPE>>(start_points[0][j]).right) {
                            tick = true;
                            left = p.left;
                        } else {
                            left = p.right;
                        }
                        ret[i][j] = Point<DISCRETE_TYPE>(left, left, left + 1);
                    });
                    j++;
                } while (tick);
            }
            return ret;
        };

        std::vector<points_vector> partition_vector = generatePartitions(partitions);
        for (auto partition : partition_vector) {
            MPI::COMM_WORLD.Isend(&partition.front(), numDimensions * sizeof(partition[0]), MPI::CHAR, i, DATA_TAG);
            i = (i + 1) % numprocs; if (i == 0) ++i;
        }
        for (auto i = 1; i < numprocs; i++) {
            std::clog << "Sending DIE_TAG to process " << i << std::endl;
            MPI::COMM_WORLD.Isend(&partition_vector[0].front(), numDimensions*sizeof(partition_vector[0][0]), MPI::CHAR, i, DIE_TAG);
        }
        auto counter = 0;
        while (counter < partitions) {
            double candidateBest;
            unsigned long long candidateFunctionCalls, candidateTotalFunctionCalls;
            points_vector candidateData(numDimensions);
            MPI::COMM_WORLD.Recv(&candidateFunctionCalls, 1, MPI::UNSIGNED_LONG_LONG, MPI::ANY_SOURCE, MPI::ANY_TAG, stat)
            MPI::COMM_WORLD.Recv(&candidateTotalFunctionCalls, 1, MPI::UNSIGNED_LONG_LONG, stat.Get_source(), MPI::ANY_TAG, stat)
            MPI::COMM_WORLD.Recv(&candidateBest, 1, MPI::DOUBLE, stat.Get_source(), MPI::ANY_TAG, stat)
            MPI::COMM_WORLD.Recv(&candidateData.front(), numDimensions * sizeof(candidateData[0]), MPI::CHAR, stat.Get_source(), MPI::ANY_TAG, stat)
            if (candidateBest < globalBest) {
                globalBest = candidateBest;
                globalData = candidateData;
            }
            function_calls += candidateFunctionCalls;
            total_function_calls += candidateTotalFunctionCalls;
            counter++;
        }
    } else {
        std::string basename = std::to_string(rank) + ".rtree";
        SpatialIndex::IStorageManager* diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(basename, 4096);
        SpatialIndex::StorageManager::IBuffer* file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
        SpatialIndex::id_type indexIdentifier;
        tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 100, 100, /*points.size()*/2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

        while (1) {
            points_vector data(numDimensions);
            REAL_TYPE f
            MPI::COMM_WORLD.Recv(&data.front(), dim*sizeof(data[0]), MPI::CHAR, 0, MPI::ANY_TAG, stat); 
            if (stat.Get_tag() == DIE_TAG) break;

            std::tie(output, f) = particle_swarm(rosenbrock, data);
            std::tie(output, f) = mesh_search(rosenbrock, output);
            MPI::COMM_WORLD.Isend(function_calls, 1, MPI::UNSIGNED_LONG_LONG, 0, 0);
            MPI::COMM_WORLD.Isend(total_function_calls, 1, MPI::UNSIGNED_LONG_LONG, 0, 0);
            MPI::COMM_WORLD.Isend(f, 1, MPI::DOUBLE, 0, 0);
            MPI::COMM_WORLD.Isend(&output.front(), output.size() * sizeof(data[0]), MPI::CHAR, 0, DIE_TAG);
        }
    }

    MPI::Finalize();
    return 0;

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
