#pragma once
#include "common.h"
#include "point.h"

class Particle {
public:
    points_vector current;
    points_vector best;
    std::vector<double> velocity;
    REAL_TYPE value;
    Particle() {};
    Particle(points_vector);
    void move();
    void updateVelocity(double, double, double, const points_vector&);
};
std::tuple<points_vector, REAL_TYPE> particle_swarm(std::function<REAL_TYPE(const points_vector&)>, const std::vector<points_vector>&);
