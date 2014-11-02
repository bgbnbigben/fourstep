#include "particle.h"

Particle::Particle(points_vector start) : current(start), best(start), value(std::numeric_limits<REAL_TYPE>::max()) {
    std::random_device rd;
    std::mt19937 gen(rd());
    velocity.resize(start.size());
    for (auto i = 0u; i < start.size(); i++) {
        velocity[i] = match(start[i], [&](Point<REAL_TYPE> p) {
                return extractReal(current[i]) - std::uniform_real_distribution<double>(p.left, p.right)(gen);
            }, [&](Point<DISCRETE_TYPE> p) {
                return (double)extractDiscrete(current[i]) - std::uniform_real_distribution<double>(p.left, p.right)(gen);
            }) / 2.;
    }
}
void Particle::move() {
    for (auto i = 0u; i < current.size(); i++) {
        match(current[i], [&](const Point<REAL_TYPE>& p) {
                REAL_TYPE* t = const_cast<REAL_TYPE*>(&p.val);
                *t = std::min(std::max(p.val + velocity[i], p.left), p.right);
            }, [&](const Point<DISCRETE_TYPE>& p) {
                DISCRETE_TYPE* t = const_cast<DISCRETE_TYPE*>(&p.val);
                *t = std::lrint(std::min(std::max((double)p.val + velocity[i], (double)p.left), (double)p.right));
        });
    }
}
void Particle::updateVelocity(double w, double r1, double r2, const points_vector& bestX) {
    for (auto i = 0u; i < velocity.size(); i++) {
        match(bestX[i], [&](Point<REAL_TYPE>) {
                velocity[i] = w * velocity[i] + r1 * (boost::get<Point<REAL_TYPE>>(best[i])() - velocity[i]) + r2 * (boost::get<Point<REAL_TYPE>>(bestX[i])() - velocity[i]);
            }, [&](Point<DISCRETE_TYPE>) {
                velocity[i] = w * velocity[i] + r1 * (boost::get<Point<DISCRETE_TYPE>>(best[i])() - velocity[i]) + r2 * (boost::get<Point<DISCRETE_TYPE>>(bestX[i])() - velocity[i]);
            });
    }
}

std::tuple<points_vector, REAL_TYPE> particle_swarm(std::function<REAL_TYPE(const points_vector&)> f, const std::vector<points_vector>& x) {
    std::vector<Particle> particles(x.size());
    for (auto i = 0u; i < x.size(); i++) {
        particles[i] = Particle(x[i]);
    }
    REAL_TYPE bestF = std::numeric_limits<REAL_TYPE>::max();
    points_vector bestX(x[0]);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0, 1);
    auto C1 = 1.496;
    auto C2 = 1.496;
    auto w_max = 0.7298;
    auto w_min = 0.3;
    auto w = w_max;
    auto NUM_STEPS = 100;
    for (auto steps = 0; steps < NUM_STEPS; steps++) {
        auto decrease_stage = 3*NUM_STEPS / 4;
        if (steps < decrease_stage) {
            w = w_min + (w_max - w_min) * (decrease_stage - steps) / decrease_stage;
        } else {
            w = w_min;
        }
        for (auto& particle: particles) {
            auto newVal = f(particle.current);
            if (newVal < particle.value) {
                particle.best = particle.current;
            }
            particle.value = newVal;
            if (particle.value < bestF) {
                bestF = particle.value;
                bestX = particle.current;
            }
            particle.move();

            auto r1 = C1 * uniform(gen);
            auto r2 = C2 * uniform(gen);
            particle.updateVelocity(w, r1, r2, bestX);
        }
    }
    return std::make_tuple(bestX, bestF);
}
