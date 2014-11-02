#pragma once
#include <boost/variant.hpp>
#include "inline_variant.hpp"
#define DISCRETE_TYPE long long
#define REAL_TYPE double

typedef boost::variant<DISCRETE_TYPE, REAL_TYPE> coord_type;
