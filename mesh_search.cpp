#include "mesh_search.h"
#include <csignal>

namespace {
    template <typename T, typename U = T>
    auto min(T a, U b) -> decltype(a < b ? a : b) {
        return a < b ? a : b;
    }
}

extern sig_atomic_t done;

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
    long long constrictions = 0; // take 2*constrictions number of nodes per line. Or something.
    auto meshWidth = 1024.0;
    while (!done && ((improvement && constrictions < 20) || constrictions < 20)) {
        std::cerr << "Constriction " << constrictions << std::endl;
        improvement = false;
        points_vector test(bestX);
        meshWidth = pow(2.0, -constrictions);

        // Alright, so.
        /* 
         * x[i] is continuous.
         * Attempt to set x[i] + e_i * mesh_width. If this is OOB, attempt
         * x[i] - e_i * mesh_width, if this is OOB, we're kinda fucked and I
         * guess we just poll some minimal distance L/R? If this is the case I
         * think we repeat this procedure with mesh_width = x[i].right - x[i],
         * unless it's zero at which point we do the same thing left. this
         * should allow for some minimal increment, but we *have* to store this
         * value so that we an add the negative of this to the set v (positive
         * basis blah blah)
         *
         * x[i] is discrete.
         * Set x[i] + 1_i, if it's OOB set x[i] - 1_i. There's literally no
         * other place to go otherwise.
         *
         * Poll the new x, compared to the current centre and its best.
         *
         * If this point is the lowest on the mesh, take mesh_width / 2 for the
         * continuous variables
         *
         * For the discrete variables, increase the `hamming distnace' such
         * that we gain the next phi(hd) directions. In theory the current
         * discrete code works for this but it's a giant flaming pile of shit
         * that needs to be cleaned up.
         */
        std::vector<coord_type> units(bestX.size());
        for (auto i = 0u; i < bestX.size(); i++) {
            match(bestX[i], [&] (const Point<REAL_TYPE>& p) {
                    REAL_TYPE dist;
                    if (std::fabs(p.right - p() - meshWidth) < std::numeric_limits<REAL_TYPE>::epsilon()) {
                        if (std::fabs(p() - p.left - meshWidth) > std::numeric_limits<REAL_TYPE>::epsilon())
                            dist = -meshWidth;
                        else if (std::fabs(p.right - p() - p() - p.left) > std::numeric_limits<REAL_TYPE>::epsilon()) 
                            dist = (p.right - p());
                        else
                            dist = -(p() - p.left);
                    } else {
                        dist = meshWidth;
                    }
                    // YOLO always move right.
                    units[i] = dist;
                }, [&] (const Point<DISCRETE_TYPE>& p) {
                    DISCRETE_TYPE dist = 1;
                    if (p() == p.right)
                        dist = -1;
                    units[i] = dist;
                });
        }
        auto clamp = [](auto a, auto l, auto r) {
            return std::max(std::min(a, r), l);
        };

        // +1 for the negative sum.
        for (auto i = 0u; i < units.size() + 1; i++) {
            // make ALL THE COPIES
            auto currentX = bestX;
            if (i < units.size()) {
                assert(units[i].which() == currentX[i].which());
                match(currentX[i], [&] (const Point<REAL_TYPE>& p) {
                        currentX[i] = Point<REAL_TYPE>(clamp(p() + boost::get<REAL_TYPE>(units[i]), p.left, p.right), p.left, p.right);
                    }, [&] (const Point<DISCRETE_TYPE>& p) {
                        currentX[i] = Point<DISCRETE_TYPE>(clamp(p() + boost::get<DISCRETE_TYPE>(units[i]), p.left, p.right), p.left, p.right);
                    });
            } else {
                // Don't need to assert since if the above passes then this
                // will also pass.
                for (auto j = 0u; j < units.size(); j++) {
                    match(currentX[j], [&] (const Point<REAL_TYPE>& p) {
                            currentX[j] = Point<REAL_TYPE>(clamp(p() - boost::get<REAL_TYPE>(units[j]), p.left, p.right), p.left, p.right);
                        }, [&] (const Point<DISCRETE_TYPE>& p) {
                            currentX[j] = Point<DISCRETE_TYPE>(clamp(p() - boost::get<DISCRETE_TYPE>(units[j]), p.left, p.right), p.left, p.right);
                        });
                }
            }
            if (i < currentX.size())
                match(currentX[i], [&](Point<REAL_TYPE> p) {
                        assert(p() >= p.left && p() <= p.right && "P is out of bounds.....");
                    }, [&](Point<DISCRETE_TYPE> p) {
                        assert(p() >= p.left && p() <= p.right && "P is out of bounds.....");
                    });
            auto currentF = f(currentX);
            if (currentF < bestF) {
                std::cout << "Found " << bestF << std::endl;
                bestF = currentF;
                bestX = currentX;
                improvement = true;
                // GREEEEEEEEEEEEEEEED 
                break;
            }
        }
        if (improvement == false)
            constrictions++;
        else
            constrictions = 0;
    }
    if (done) std::cerr << "RAN OUT OF TIME" << std::endl;

    std::cerr << bestF << std::endl;
    return std::make_tuple(bestX, bestF);
}

