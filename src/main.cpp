/* CelMech
 * Copyright (C) 2019, Szilard Ledan <szledan@gmail.com>
 * All rights reserved.
 *
 * APACHE LICENSE, VERSION 2.0
 * https://www.apache.org/licenses/LICENSE-2.0
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace cm {

struct Vec3 {
    Vec3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
    double x, y, z;

    double length() const { return std::sqrt(x * x + y * y + z * z); }
    Vec3& operator+=(const Vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    Vec3& operator-=(const Vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
};

Vec3 operator+(Vec3 lhs, const Vec3& rhs) { lhs += rhs; return lhs; }
Vec3 operator-(Vec3 lhs, const Vec3& rhs) { lhs -= rhs; return lhs; }
Vec3 operator*(const double c, const Vec3& rhs) { return Vec3(c * rhs.x, c * rhs.y, c * rhs.z); }
std::ostream& operator<<(std::ostream& os, const Vec3& obj)
{
    os << obj.x << "," << obj.y << "," << obj.z;
    return os;
}

struct Object {
    Object(double m_, Vec3 r_, Vec3 v_) : m(m_ ? m_ : 1.0), r(r_), v(v_) {}
    Vec3 r, v, a;
    double m;
};

std::ostream& operator<<(std::ostream& os, const Object& obj)
{
    os << obj.m << ",";
    os << obj.r << ",";
    os << obj.v;
    return os;
}


void visualize(const std::vector<cm::Object>& objects)
{
    const int SIZE_X = 40;
    const int SIZE_Y = 20;
    static bool isFirst = true;

    if (!isFirst) {
        std::cout << "\033[" + std::to_string(SIZE_Y) + "F";
    }

    std::vector<std::string> map(SIZE_Y, std::string(SIZE_X, '.'));

    for (int i = 0; i < objects.size(); ++i) {
        const cm::Object& o = objects[i];
        const int x = int(o.r.x);
        const int y = int(o.r.y);
        if (x > (-SIZE_X / 2) && x < (SIZE_X / 2) && y > (-SIZE_Y / 2) && y < (SIZE_Y / 2)) {
            map[y + SIZE_Y / 2][x + SIZE_X / 2] = 'X';
        }
    }

    for (int i = 0; i < map.size(); ++i) {
        std::cout << map[i] << std::endl;
    }

    isFirst = false;
}

} // namespace cm

int main(int argc, char* argv[])
{
    const double k = 1;
    const double kk = k * k;
    int stepNum = 10000;

    if (argc > 1) {
        stepNum = std::atoi(argv[1]);
    }

    std::vector<cm::Object> objects = {
        cm::Object(0.01, cm::Vec3(-10, 0, 0), cm::Vec3(0, 0.001, 0)),
        cm::Object(0.01, cm::Vec3(10, 0, 0), cm::Vec3(0, -0.001, 0)),
    };

    std::cout << "step_num: " << stepNum << std::endl;

    // Calculate
    for (int it = 0; it < stepNum; ++it) {
        for (int i = 0; i < objects.size(); ++i) {
            cm::Object& co = objects[i];
            co.a = cm::Vec3();
            for (int j = 0; j < objects.size(); ++j) {
                if (i != j) {
                    cm::Object& o = objects[j];
                    const double d = (o.r - co.r).length();
                    co.a += kk * ((co.m * o.m) / (d * d * d)) * ((o.r - co.r));
                }
            }
            co.v += co.a;
            co.r += co.v;
        };
        visualize(objects);
    }

    return 0;
}
