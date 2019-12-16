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
#include <unistd.h>

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

struct System {
    // double U = 1/2 * SUM(i; SUM(j, i!=j; k^2 * m_i * m_j / r_ij));
    // Vec3 r = SUM(i; m_i * r_i) / SUM(i; r_i);
    // Vec3 m = SUM(i; m_i);
    // m * r = a * t + b; Vec3 a, b; double t;
};

std::ostream& operator<<(std::ostream& os, const Object& obj)
{
    os << obj.m << ",";
    os << obj.r << ",";
    os << obj.v;
    return os;
}

void moveCursor(int r, int c)
{
    static int rw = 0;
    static int cl = 0;
    if (r - rw) {
        std::cout << "\033[" << std::abs(r - rw) << (r - rw > 0 ? "B" : "A");
        rw = r;
    }
    if (c - cl) {
        std::cout << "\033[" << std::abs(c - cl) << (c - cl > 0 ? "C" : "D");
        cl = c;
    }
}

void putChar(int r, int c, char ch)
{
    moveCursor(r, c);
    std::cout << ch << "\033[D";
}

void putStr(int r, int c, std::string str)
{
    moveCursor(r, c);
    std::cout << str;
    for (int i = 0; i < str.size(); ++i) {
        std::cout << "\033[D";
    }
}

const int SIZE_X = 40;
const int SIZE_Y = 30;

int create()
{
    std::cout << std::string(SIZE_Y, '\n') << std::endl;
    std::cout << "\033[" + std::to_string(SIZE_Y) + "F";
    return 0;
}

void visualize(const std::vector<cm::Object>& objects, char ch, int it)
{
    for (int i = 0; i < objects.size(); ++i) {
        const cm::Object& o = objects[i];
        const int x = int(o.r.x);
        const int y = int(o.r.y);
        if (x > (-SIZE_X / 2) && x < (SIZE_X / 2) && y > (-SIZE_Y / 2) && y < (SIZE_Y / 2)) {
            putChar(y + SIZE_Y / 2, x + SIZE_X / 2, ch);
        }
    }
    putStr(0, SIZE_X + 1, "it: " + std::to_string(it));
}

} // namespace cm

int main(int argc, char* argv[])
{
    const double G = 1; // Standard constant: G = 6.67430(15) * 10^−11 (m^3 * kg^−1 * s^−2) = 6.6743e-11
    const double GM_sun = 1; // Standard gravitational parameter : GM_sun = 1.32712440018(9) * 10^20 (m^3 * s^−2) =  1.32712440018e20
    const double AU = 1; // ...
    const double k = 1; // Gaussian gravitational constant : k = (G * M_sun)^0.5 * AU^−1.5 = 0.01720209895 rad/day
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
    cm::create();
    std::cout << "\e[?25l";
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
        if (it % 10 == 0) {
            visualize(objects, 'o', it);
        }
        usleep(100);
        if (it % 10 == 0) {
            visualize(objects, '.', it);
        }
    }
    visualize(objects, 'o', stepNum);
    cm::moveCursor(cm::SIZE_Y, 0);
    std::cout << "\e[?25h";

    return 0;
}
