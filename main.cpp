#include "common.h"
#include "libspatialindex/include/spatialindex/SpatialIndex.h"
#include "clp/build/include/coin/CoinMpsIO.hpp"
#include <iostream>
#include <random>
#include "point.h"
#include "mesh_search.h"
#include <mpi.h>
#include "parser.h"
#include <csignal>


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
int rank;
volatile sig_atomic_t done = 0;

void term (int signum) {
    done = 1;
}

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
    SpatialIndex::Region point_region(&low.front(), &high.front(), x.size());
    SpatialIndex::Point point(&casted.front(), x.size());
    CheckerVisitor visitor;
    tree->pointLocationQuery(point, visitor);
    if (visitor.found()) {
        return visitor.data();
    }

    auto start = extractReal(x[0]);
    auto ret_val = (1.0 - start)*(1 - start);
    std::cout << "Testing ("; 
    for (auto i = x.size() - 1; i > 0; i--) {
        std::cout << extractReal(x[i]) << ", ";
        auto curr = extractReal(x[i]);
        auto next = extractReal(x[i-1]);
        ret_val += 100*(curr - next*next)*(curr - next*next);
    }
    std::cout << extractReal(x[0]) << ")" << std::endl;

    tree->insertData(sizeof(REAL_TYPE), reinterpret_cast<const byte*>(&ret_val), point_region, id++);
    ++function_calls;
    return ret_val;
}

void test() {
    std::string basename = "rtrees/test.rtree";
    SpatialIndex::IStorageManager* diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(basename, 4096);
    SpatialIndex::StorageManager::IBuffer* file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    SpatialIndex::id_type indexIdentifier;
    tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 100, 100, 2 /*numDimensions*/, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

    coord_type low[2]; low[0] = 0ll, low[1] = 0ll;
    coord_type high[2]; high[0] = 1ll, high[1] = 1ll;
    coord_type higher[2]; higher[0] = 0ll, higher[1] = 1ll;
    //assert(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(higher, 2u)));
    //assert(!(varRegion(varPoint(low, 2u), varPoint(higher, 2u)) == varRegion(varPoint(low, 2u), varPoint(high, 2u))));
    points_vector x = {Point<DISCRETE_TYPE>(0, 0, 1), Point<DISCRETE_TYPE>(0, 0, 1)};
    mesh_search(rosenbrock, x);
}

#define DATA_TAG 0xd474
#define DIE_TAG 0xd1ed1e

std::vector<double> lower;
std::vector<double> upper;
std::vector<double> objCoeffs;
std::vector<std::vector<double>> constraints;
std::vector<char> type;
double objOffset;

REAL_TYPE mpsFunction(const points_vector& x) {
try{
    ++total_function_calls;
    std::vector<double> low(x.size());
    std::vector<double> high(x.size());
    std::vector<double> casted(x.size());
    # pragma omp parallel for 
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
    SpatialIndex::Region point_region(&low.front(), &high.front(), x.size());
    SpatialIndex::Point point(&casted.front(), x.size());
    CheckerVisitor visitor;
    tree->pointLocationQuery(point, visitor);
    if (visitor.found()) {
        return visitor.data();
    }

    auto ret = objOffset;
    auto penalty = 1;
    # pragma omp parallel for 
    for (auto i = 0u; i < x.size(); i++) {
        ret += match(x[i], [&](Point<REAL_TYPE> p) {
            return p() * objCoeffs[i];
        }, [&](Point<DISCRETE_TYPE> p) {
            return p() * objCoeffs[i];
        });
    }

    # pragma omp parallel for 
    for (auto i = 0u; i < constraints.size(); i++) {
        auto constraint = 0.0;
        assert(constraints[i].size() == x.size());
        for (auto j = 0u; j < x.size(); j++) {
            constraint += match(x[j], [&](Point<REAL_TYPE> p) {
                return p() * constraints[i][j];
            }, [&](Point<DISCRETE_TYPE> p) {
                return p() * constraints[i][j];
            });
        }
        switch (type[i]) {
            case 'G':
                if (constraint - lower[i] < 0) {
                    ret += penalty * (lower[i] - constraint);
                }
                break;
            case 'L':
                if (constraint - upper[i] > 0) {
                    ret += penalty * (constraint - upper[i]);
                }
                break;
            case 'R':
                if (constraint - lower[i] < 0) {
                    ret += penalty * (lower[i] - constraint);
                }
                if (constraint - upper[i] > 0) {
                    ret += penalty * (upper[i] - constraint);
                }
                break;
            case 'E':
                if (constraint != upper[i]) {
                    ret += penalty * (upper[i] - constraint) * ((upper[i] > constraint) ? 1. : -1.);
                }
                break;
        }
    }

    tree->insertData(sizeof(REAL_TYPE), reinterpret_cast<const byte*>(&ret), point_region, (id++ % 1000));
    ++function_calls;

    return ret;
} catch (Tools::IllegalArgumentException& e) {
    std::cerr << e.what() << std::endl;
    throw;
}
}

int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);

    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);

    std::vector<points_vector> start_points(1);
    start_points[0] = parseGams(argv[1]);

    /*
    CoinMpsIO m; 
    m.readMps(argv[1], "");

    objOffset = m.objectiveOffset();
    objCoeffs.resize(m.getNumCols());
    lower.resize(m.getNumRows());
    upper.resize(m.getNumRows());
    type.resize(m.getNumRows());
    constraints.resize(m.getNumRows());
    std::cout << "There are " << m.getNumRows() << " constraints " << std::endl;
    for (auto i = 0u; i < constraints.size(); i++) {
        constraints[i].resize(m.getNumCols());
        type[i] = m.getRowSense()[i];
        switch (type[i]) {
            case 'G':
                lower[i] = m.getRowLower()[i];
                break;
            case 'E':
            case 'L':
                upper[i] = m.getRowUpper()[i];
                break;
            case 'R':
                lower[i] = m.getRowLower()[i];
                upper[i] = m.getRowUpper()[i];
                break;
        }
    }

    for (auto i = 0; i < m.getNumCols(); i++) {
        objCoeffs[i] = m.getObjCoefficients()[i];
    }

    auto matrix = m.getMatrixByRow();
    for (auto i = 0, j = 0; i < matrix->getNumElements(); i++) {
        if (i == 0)
            constraints[j][matrix->getIndices()[0]] = matrix->getElements()[0];
        else if (matrix->getIndices()[i] < matrix->getIndices()[i-1]) {
            constraints[++j][matrix->getIndices()[i]] = matrix->getElements()[i];
        } else {
            constraints[j][matrix->getIndices()[i]] = matrix->getElements()[i];
        }
    }

    std::vector<points_vector> start_points(1);
    for (auto i = 0; i < m.getNumCols(); i++) {
        if (m.integerColumns() != nullptr && m.integerColumns()[i]) {
            auto upper = m.getColUpper()[i];
            if (upper == std::numeric_limits<double>::max()) upper = 100000;
            start_points[0].push_back(Point<DISCRETE_TYPE>(m.getColLower()[i], m.getColLower()[i], upper));
        } else {
            auto upper = m.getColUpper()[i];
            if (upper == std::numeric_limits<double>::max()) upper = 100000.;
            start_points[0].push_back(Point<REAL_TYPE>(m.getColLower()[i], m.getColLower()[i], upper));
        }
    }
    */


    MPI::Status stat;
    auto numprocs = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    auto numDimensions = start_points[0].size();
    std::cout << "There are " << numDimensions << " dimensions" << std::endl;
    // Fuck it, bin process 0 as only a data marshall. Whatever.
    if (rank == 0) {
        auto globalBest = std::numeric_limits<double>::max();
        points_vector globalData;
        auto partitions = 64;
        auto generatePartitions = [&](auto numPartitions) {
            auto actualPartitions = numPartitions;
            std::vector<int> splits(numDimensions, 1);
            std::vector<int> splitIndex; splitIndex.reserve(numDimensions);
            std::cout << "cont are ";
            for (auto i = 0u; i < start_points[0].size(); i++) {
                match(start_points[0][i], [](Point<DISCRETE_TYPE>) {
                    }, [&](Point<REAL_TYPE>) {
                        std::cout << i << " ";
                        splitIndex.push_back(i);
                    });
            }
            std::cout << std::endl;
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(splitIndex.begin(), splitIndex.end(), g);
            auto found = false;
            for (auto i = 10; i > 1 && !found; i--) {
                // attempt to get i ^ x * (i-1) ^ (splitIndex.size() - x) == numPartitions. 
                // Solve for the largest possible i. Solved if x > 0
                for (auto j = std::log(numPartitions) / std::log(i); j > 0; j--) {
                    if (std::pow(i, j) * std::pow(i - 1, splitIndex.size() - j) <= numPartitions) {
                        actualPartitions = std::pow(i, j) * std::pow(i - 1, splitIndex.size() - j);
                        for (auto k = 0; k < splitIndex.size(); k++) {
                            if (k <= j)
                                splits[splitIndex[k]] = j;
                            else
                                splits[splitIndex[k]] = j-1;
                        }
                        found = true;
                        std::cout << "Solving i^x * (i-1) ^ (nd - x) == np with " << i << "^" << j << "*" << i-1 << "^" << splitIndex.size() - j << std::endl;
                        break;
                    }
                }
            }
            std::cout << "Getting " << actualPartitions << " when we wanted " << numPartitions << std::endl;
            std::vector<points_vector> ret(actualPartitions);
            ret[0] = start_points[0];
            for (auto i = 0u; i < actualPartitions; i++) {
                if (i > 0)
                    ret[i] = ret[i-1];
                auto tick = false;
                auto j = 0;
                do {
                    assert(start_points[0][j].which() == ret[i][j].which());
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
                } while (tick && j < ret[i].size());
            }
            return ret;
        };

        std::vector<points_vector> partition_vector = generatePartitions(partitions);
        auto proc_i = 1;
        for (auto& partition : partition_vector) {
            std::cout << "Sending to process " << proc_i << std::endl;
            MPI::COMM_WORLD.Isend(&partition.front(), numDimensions * sizeof(partition[0]), MPI::CHAR, proc_i, DATA_TAG);
            proc_i = (proc_i + 1) % numprocs; if (proc_i == 0) ++proc_i;
        }
        for (auto i = 1; i < numprocs; i++) {
            std::cout << "Sending DIE_TAG to process " << i << std::endl;
            MPI::COMM_WORLD.Isend(&partition_vector[0].front(), numDimensions*sizeof(partition_vector[0][0]), MPI::CHAR, i, DIE_TAG);
        }
        auto counter = 0;
        while (counter < partition_vector.size()) {
            double candidateBest;
            unsigned long long candidateFunctionCalls, candidateTotalFunctionCalls;
            // Don't care about the value; it'll be overwritten
            points_vector candidateData(numDimensions, Point<REAL_TYPE>(0., 0., 0.));
            std::cout << "Master blocking while trying to receive function calls back" << std::endl;
            MPI::COMM_WORLD.Recv(&candidateFunctionCalls, 1, MPI::UNSIGNED_LONG_LONG, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
            std::cout << "Got " << candidateFunctionCalls << " calls\n";
            MPI::COMM_WORLD.Recv(&candidateTotalFunctionCalls, 1, MPI::UNSIGNED_LONG_LONG, stat.Get_source(), MPI::ANY_TAG, stat);
            std::cout << "Out of " << candidateTotalFunctionCalls << " calls\n";
            MPI::COMM_WORLD.Recv(&candidateBest, 1, MPI::DOUBLE, stat.Get_source(), MPI::ANY_TAG, stat);
            std::cout << "Candiate best was " << candidateBest << std::endl;
            MPI::COMM_WORLD.Recv(&candidateData.front(), numDimensions * sizeof(candidateData[0]), MPI::CHAR, stat.Get_source(), MPI::ANY_TAG, stat);
            std::cout << "Received the point " << std::endl;
            if (candidateBest < globalBest) {
                globalBest = candidateBest;
                globalData = candidateData;
            }
            function_calls += candidateFunctionCalls;
            total_function_calls += candidateTotalFunctionCalls;
            counter++;
        }

        std::cout << "GLOBAL BEST: " << globalBest << "\nFUNC_CALLS: " << function_calls << "\n:TOTAL_CALLS: " << total_function_calls << std::endl;
    } else {
        try {
        std::string basename = "rtrees/" + std::string(argv[1]) + "-" + std::to_string(rank) + ".rtree";
        std::cerr << "basename is " << basename << std::endl;
        SpatialIndex::IStorageManager* diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(basename, 4096);
        std::cerr << "diskfile done" << std::endl;
        SpatialIndex::StorageManager::IBuffer* file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
        std::cerr << "rev buffer " << std::endl;
        SpatialIndex::id_type indexIdentifier;
        tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 1000, 1000, numDimensions, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
        } catch (Tools::IllegalArgumentException& e) {
            std::cerr << e.what() << std::endl;
            throw e;
        }
        auto numSwarmParticles = 30u;

        while (1) {
            REAL_TYPE f;
            // As above, who cares; it'll be overwritten
            points_vector start_partition(numDimensions, Point<REAL_TYPE>(0., 0., 0.));
            points_vector output(numDimensions, Point<REAL_TYPE>(0., 0., 0.));
            MPI::COMM_WORLD.Recv(&start_partition.front(), numDimensions*sizeof(start_partition[0]), MPI::CHAR, 0, MPI::ANY_TAG, stat); 
            std::cout << "Process " << rank << " received data of size " << numDimensions * sizeof(start_partition[0]) << std::endl;
            if (stat.Get_tag() == DIE_TAG) break;

            std::random_device rd;
            std::mt19937 gen(rd());
            std::vector<points_vector> data(numSwarmParticles);
            for (auto i = 0u; i < numSwarmParticles; i++) {
                data[i].resize(numDimensions, Point<REAL_TYPE>(0., 0., 0.));
                for (auto j = 0u; j < numDimensions; j++)  {
                    data[i][j] = start_partition[j];
                }
                for (auto j = 0u; j < numDimensions; j++) {
                    match(data[i][j], [&](Point<REAL_TYPE> p) {
                        std::uniform_real_distribution<REAL_TYPE> real(p.left, p.right);
                        data[i][j] = Point<REAL_TYPE>(real(gen), p.left, p.right);
                    }, [&](Point<DISCRETE_TYPE> p) {
                        std::uniform_int_distribution<DISCRETE_TYPE> disc(p.left, p.right);
                        data[i][j] = Point<DISCRETE_TYPE>(disc(gen), p.left, p.right);
                    });
                }
            }

            std::cout << "swarm " << std::endl;
            //std::tie(output, f) = particle_swarm(mpsFunction, data);
            std::tie(output, f) = particle_swarm(gamsFunc, data);
            std::cout << "mesh " << std::endl;
            std::tie(output, f) = mesh_search(gamsFunc, output);
            std::cout << "Sending " << f << " which is at point" << "\n";
            for (auto i = 0u; i < output.size(); i++) {
                match(output[i], [&](Point<REAL_TYPE> p) {
                    std::cout << p() << ",";
                }, [&](Point<DISCRETE_TYPE> p) {
                    std::cout << p() << ",";
                });
            }
            std::cout << std::endl;
            MPI::COMM_WORLD.Send(&function_calls, 1, MPI::UNSIGNED_LONG_LONG, 0, 0);
            MPI::COMM_WORLD.Send(&total_function_calls, 1, MPI::UNSIGNED_LONG_LONG, 0, 0);
            MPI::COMM_WORLD.Send(&f, 1, MPI::DOUBLE, 0, 0);
            MPI::COMM_WORLD.Send(&output.front(), output.size() * sizeof(output[0]), MPI::CHAR, 0, DIE_TAG);
            std::cout << "sent" << std::endl;
        }
    }

    MPI::Finalize();
    return 0;
}
