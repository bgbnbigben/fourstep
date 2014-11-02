#pragma once
#include "common.h"
#include "point.h"
std::tuple<points_vector, REAL_TYPE> mesh_search(std::function<REAL_TYPE(const points_vector&)>, const points_vector&);
