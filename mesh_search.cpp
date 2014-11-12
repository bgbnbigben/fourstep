#include "mesh_search.h"

namespace {
    template <typename T, typename U = T>
    auto min(T a, U b) -> decltype(a < b ? a : b) {
        return a < b ? a : b;
    }
}

std::tuple<points_vector, REAL_TYPE> mesh_search(std::function<REAL_TYPE(const points_vector&)> f, const points_vector& x) {
    std::vector<bool> continuous(x.size());
    std::transform(x.begin(), x.end(), continuous.begin(),
            [](const point_type& p) { return p.type() == typeid(Point<double>); });
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
    while ((improvement && constrictions < 10) || constrictions < 10) {
        improvement ? (std::cerr << "Had an improvement\n") : (std::cerr << "Going around again\n");
        std::cerr << "Constriction " << constrictions << std::endl;
        improvement = false;
        points_vector test(bestX);

        unsigned long long numPoints = 1;
        unsigned long long nodesPerRow = pow(2, constrictions) + 1;
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
                    test[i] = std::move(Point<REAL_TYPE>(left - (right - left)/(double)nodesPerRow, p.left, p.right));
                } else {
                    test[i] = std::move(Point<REAL_TYPE>(left, p.left, p.right));
                }
            } else {
                auto& p = boost::get<Point<DISCRETE_TYPE>>(test[i]);
                const auto best = boost::get<Point<DISCRETE_TYPE>>(bestX[i]);
                auto inRow = ::min(nodesPerRow, (best.right - best.left + 1));
                auto left = ::min(inRow / 2, best.val - best.left);
                if (i == 0) {
                    test[i] = std::move(Point<DISCRETE_TYPE>(best.val - left - 1, p.left, p.right));
                } else {
                    test[i] = std::move(Point<DISCRETE_TYPE>(best.val - left, p.left, p.right));
                }
            }
        }

        std::cerr << "We have " << numPoints << " points on constriction " << constrictions << std::endl;
        for (auto pointNum = 0u; pointNum < numPoints; pointNum++) {
            int i = 0;
            bool tick = false;
            do {
                if (continuous[i]) {
                    auto& current = boost::get<Point<REAL_TYPE>>(test[i]);
                    test[i] = std::move(Point<REAL_TYPE>(current.val + (current.right - current.left) / (double)nodesPerRow, current.left, current.right));
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
                    auto inRow = ::min(nodesPerRow, (best.right - best.left + 1));
                    auto left = ::min(inRow / 2, best.val - best.left);
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

