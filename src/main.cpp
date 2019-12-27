/* CelMech
 * Copyright (C) 2019, Szilard Ledan <szledan@gmail.com>
 * All rights reserved.
 *
 * APACHE LICENSE, VERSION 2.0
 * https://www.apache.org/licenses/LICENSE-2.0
 */

#include <assert.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ncurses.h>
#include <string>
#include <vector>
#include <unistd.h>

namespace cm {

struct Vec3 {
    Vec3(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
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
    Object(const double m_, const Vec3& r_, const Vec3& v_) : m(m_ ? m_ : 1.0), r(r_), v(v_) {}
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

struct System {
    typedef typename std::vector<Object> Objects;

    System& addObject(Object obj)
    {
        _objects.push_back(obj);
        return *this;
    }

    System& addObject(const double m, const Vec3& r, const Vec3& v) { return addObject(Object(m, r, v)); }

    const size_t size() const { return _objects.size(); }
    const bool empty() const { return _objects.size() == 0; }

    const Objects& objects() const { return _objects; }
    Objects& objects() { return _objects; }

private:
    std::vector<Object> _objects;

    // double U = 1/2 * SUM(i; SUM(j, i!=j; k^2 * m_i * m_j / r_ij));
    // Vec3 r = SUM(i; m_i * r_i) / SUM(i; r_i);
    // Vec3 m = SUM(i; m_i);
    // m * r = a * t + b; Vec3 a, b; double t;
};

std::ostream& operator<<(std::ostream& os, const System& system)
{
    for (const auto& o : system.objects()) {
        os << o << std::endl;
    }
    return os;
}

struct Universe {
    double G = 1; // Standard constant: G = 6.67430(15) * 10^−11 (m^3 * kg^−1 * s^−2) = 6.6743e-11
    double GM_sun = 1; // Standard gravitational parameter : GM_sun = 1.32712440018(9) * 10^20 (m^3 * s^−2) =  1.32712440018e20
    double AU = 1; // ...
    double k = 1; // Gaussian gravitational constant : k = (G * M_sun)^0.5 * AU^−1.5 = 0.01720209895 rad/day
    double kk = k * k;
};

struct Simulator {
    Simulator& setSystem(System& system)
    {
        _system = system;
        return *this;
    }

    Simulator& setUniverse(Universe& universe)
    {
        _universe = universe;
        return *this;
    }

    const System& simulate()
    {
        System::Objects& objects = _system.objects();
        for (int i = 0; i < objects.size(); ++i) {
            cm::Object& co = objects[i];
            co.a = cm::Vec3();
            for (int j = 0; j < objects.size(); ++j) {
                if (i != j) {
                    cm::Object& o = objects[j];
                    const double d = (o.r - co.r).length();
                    co.a += _universe.kk * ((co.m * o.m) / (d * d * d)) * ((o.r - co.r));
                }
            }
            co.v += co.a;
            co.r += co.v;
        };

        return _system;
    }

private:
    System _system;
    Universe _universe;
};

template<typename T>
struct Visualizer {
    T& setSimulator(Simulator& simulator_)
    {
        simulator = simulator_;
        return static_cast<T&>(*this);
    }

    virtual void run(int argc, char* argv[]) = 0;

    Simulator simulator;
};

struct ConsoleVisualizer : public Visualizer<ConsoleVisualizer> {
    ConsoleVisualizer& setViewport(unsigned int width, unsigned int height)
    {
        _width = width;
        _height = height;
        return *this;
    }

    virtual void run(int argc, char* argv[])
    {
        createCanvas();
        char ch;
        while (ch != 27) {
            ch = getch();
            const System& system = simulator.simulate();;
            visualize(system, 'o');
            usleep(200);
        }
        finalize();
    }

    void finalize()
    {
        refresh();
        endwin();
//        putStr(0, _width + 1, "it: " + std::to_string(it));
//        moveCursor(_height, 0);
//        std::cout << "\e[?25h"; // Cursor ON
    }

    void visualize(const System& system, const char ch)
    {
        for (int i = 0; i < system.size(); ++i) {
            const Object& o = system.objects()[i];
            const int x = int(o.r.x);
            const int y = int(o.r.y);
            if (x > (-_width / 2) && x < (_width / 2) && y > (-_height / 2) && y < (_height / 2)) {
                putChar(y + _height / 2, x + _width / 2, ch);
            }

            mvprintw(_height - i, 0, "%d: %10f %10f", i, o.r.x, o.r.y);
        }
        move(0, 0);
    }

private:
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
        mvprintw(r, c, "%c", ch);
//        moveCursor(r, c);
//        std::cout << ch << "\033[D";
    }

    void putStr(int r, int c, std::string str)
    {
        moveCursor(r, c);
        std::cout << str;
        for (int i = 0; i < str.size(); ++i) {
            std::cout << "\033[D";
        }
    }

    void createCanvas()
    {
        initscr();
        noecho();
        scrollok(stdscr, TRUE);
        nodelay(stdscr, TRUE);
//        std::cout << "\e[?25l"; // Cursor OFF
//        std::cout << std::string(_height, '\n') << std::endl;
//        std::cout << "\033[" + std::to_string(_height) + "F";
    }

    int _width = 40;
    int _height = 30;
    System _prevSystem;
};

} // namespace cm

int main(int argc, char* argv[])
{
    cm::System twoObj = cm::System()
            .addObject(0.01, cm::Vec3(-10, 0, 0), cm::Vec3(0, 0.001, 0))
            .addObject(0.01, cm::Vec3(10, 0, 0), cm::Vec3(0, -0.001, 0));
    cm::Simulator simulator = cm::Simulator()
            .setSystem(twoObj);

    cm::ConsoleVisualizer()
            .setSimulator(simulator)
            .setViewport(40, 30)
            .run(argc, argv);

    return 0;
}
