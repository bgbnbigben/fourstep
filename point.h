#pragma once
#include "common.h"

template <class T>
class Point { 
    public:
        const T val;
        const T left;
        const T right;

        Point() = delete;
        Point(const T v, const T l, const T r) : val(v), left(l), right(r) {};
        Point(const Point& p) : val(p.val), left(p.left), right(p.right) {};
        Point(const Point&& p) : val(p.val), left(p.left), right(p.right) {};

        Point<T>& operator=(const Point<T>& p) {
            // Urgh.
            auto t = const_cast<T*>(&this->val);
            *t = p.val;
            t = const_cast<T*>(&this->left);
            *t = p.left;
            t = const_cast<T*>(&this->right);
            *t = p.right;

            return *this;
        }

        Point<T>& operator=(const Point<T>&& p) {
            auto t = const_cast<T*>(&this->val);
            *t = p.val;
            t = const_cast<T*>(&this->left);
            *t = p.left;
            t = const_cast<T*>(&this->right);
            *t = p.right;

            return *this;
        }

        T operator()() const { return val; }
};

typedef boost::variant<Point<DISCRETE_TYPE>, Point<REAL_TYPE>> point_type;
typedef std::vector<point_type> points_vector;

REAL_TYPE extractReal(const point_type&);
DISCRETE_TYPE extractDiscrete(const point_type&);
