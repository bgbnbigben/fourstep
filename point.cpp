#include "point.h"

REAL_TYPE extractReal(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return p(); }, [](const Point<DISCRETE_TYPE>& p) { return (REAL_TYPE)p(); });
}

DISCRETE_TYPE extractDiscrete(const point_type& x) {
    return match(x, [](const Point<REAL_TYPE>& p) { return (DISCRETE_TYPE)p(); }, [](const Point<DISCRETE_TYPE>& p) { return p(); });
}
