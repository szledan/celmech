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
#include <map>
#include <ncurses.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace cm {

struct Coord {
    Coord(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
    double x, y, z;

    double length() const { return std::sqrt(x * x + y * y + z * z); }

    Coord& operator+=(const Coord& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    Coord& operator-=(const Coord& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
    Coord& operator*=(const double rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
    Coord& operator/=(const double rhs) { return this->operator*=(1.0 / rhs); }
};

Coord operator+(Coord lhs, const Coord& rhs) { lhs += rhs; return lhs; }
Coord operator-(Coord lhs, const Coord& rhs) { lhs -= rhs; return lhs; }
Coord operator*(Coord lhs, const double c) { lhs *= c; return lhs; }
Coord operator*(const double c, Coord rhs) { return rhs * c; }
Coord operator/(const double c, Coord rhs) { rhs /= c; return rhs; }
Coord operator/(Coord lhs, const double c) { lhs /= c; return lhs; }

std::ostream& operator<<(std::ostream& os, const Coord& obj)
{
    os << obj.x << "," << obj.y << "," << obj.z;
    return os;
}

typedef double Shape;
typedef double Mass;

struct Object {
    Object(const double m_ = 1.0, const Coord& r_ = 0.0, const Coord& v_ = 0.0) : m(m_ ? m_ : 1.0), r(r_), v(v_) {}
    Coord r, v, a;
    Shape R;
    Mass m;
};

std::ostream& operator<<(std::ostream& os, const Object& obj)
{
    os << obj.m << ",";
    os << obj.r << ",";
    os << obj.v;
    return os;
}

struct System {
    typedef typename std::map<int, Object> Objects;

    System& addObject(Object obj)
    {
        _objects[_uid++] = obj;
        return *this;
    }

    System& addObject(const double m, const Coord& r, const Coord& v) { return addObject(Object(m, r, v)); }

    const size_t size() const { return _objects.size(); }
    const bool empty() const { return _objects.size() == 0; }

    const Objects& objects() const { return _objects; }
    Objects& objects() { return const_cast<Objects&>(const_cast<const System*>(this)->objects()); }

private:
    Objects _objects;
    int _uid = 0;

    // double U = 1/2 * SUM(i; SUM(j, i!=j; k^2 * m_i * m_j / r_ij));
    // Vec3 r = SUM(i; m_i * r_i) / SUM(i; r_i);
    // Vec3 m = SUM(i; m_i);
    // m * r = a * t + b; Vec3 a, b; double t;
};

std::ostream& operator<<(std::ostream& os, const System& system)
{
    for (auto& o : system.objects()) {
        os << o.second << std::endl;
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

    void init()
    {
        _table.resize(_system.size() * (_system.size() - 1) / 2);
        System::Objects& objects = _system.objects();
        for (int i = 0, m = 0; i < objects.size(); ++i) {
            const cm::Object& co = objects[i];
            for (int j = i + 1; j < objects.size(); ++j, ++m) {
                const cm::Object& o = objects[j];
                _table[m].M = _universe.kk * co.m * o.m;
                _table[m].a = 0.0;
            }
        }
    }

    const System& simulate()
    {
        System::Objects& objects = _system.objects();

        for (int i = 0, m = 0; i < objects.size() - 1; ++i) {
            const cm::Object& co = objects[i];
            _table[m].a = 0.0;
            for (int j = i + 1; j < objects.size(); ++j, ++m) {
                const cm::Object& o = objects[j];
                const cm::Coord rr = o.r - co.r;
                const double d = rr.length();
                _table[m].a = _table[m].M / (d * d * d) * rr;
            }
        }

        for (int i = 0, m = 0; i < objects.size(); ++i) {
            const cm::Object& co = objects[i];
            for (int j = i + 1; j < objects.size(); ++j, ++m) {
                const cm::Object& o = objects[j];
                const cm::Coord rr = o.r - co.r;
                const double d = rr.length();
                _table[m].a = _table[m].M / (d * d * d) * rr;
            }
        }

        for (int i = 0; i < objects.size(); ++i) {
            cm::Object& co = objects[i];
            co.a = cm::Coord();
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
    struct Coef {
        double M;
        cm::Coord a;
    };

    std::vector<Coef> _table;

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
        _canvas.w = width;
        _canvas.h = height;
        return *this;
    }

    virtual void run(int argc, char* argv[])
    {
        createCanvas();
        simulator.init();
        char ch;
        while (ch != 27 && ch != 'q') {
            ch = getch();
            const System& system = simulator.simulate();
            visualize(system, 'o');
            usleep(200);
        }
        finalize();
    }

    void visualize(const System& system, const char ch)
    {
        for (int i = 0; i < system.size(); ++i) {
            const cm::Object& o = system.objects().at(i);
            const int x = int(o.r.x);
            const int y = int(o.r.y);
            if (x > (-_canvas.w / 2) && x < (_canvas.w / 2) && y > (-_canvas.h / 2) && y < (_canvas.h / 2)) {
                mvaddch(y + _canvas.h / 2, x + _canvas.w / 2, ch);
            }

            mvprintw(_canvas.h - system.size() + i, 0, "%d: %10f %10f", i, o.r.x, o.r.y);
        }

        move(_screen.h - 1, 0);
    }

private:
    void createCanvas()
    {
        initscr();
        noecho();
        scrollok(stdscr, TRUE);
        nodelay(stdscr, TRUE);

        getmaxyx(stdscr, _screen.h, _screen.w);
        _canvas = _screen;
        _canvas.h -= 1;
    }

    void finalize()
    {
        refresh();
        endwin();
    }

    struct Size {
        int w = 0;
        int h = 0;
    } _canvas, _screen;


    System _prevSystem;
};

} // namespace cm

int main(int argc, char* argv[])
{
    cm::System twoObj = cm::System()
            .addObject(0.01, cm::Coord(-10, 0, 0), cm::Coord(0, 0.001, 0))
            .addObject(0.01, cm::Coord(10, 0, 0), cm::Coord(0, -0.001, 0))
            .addObject(0.0001, cm::Coord(15, 0, 0), cm::Coord(0, 0.0005, 0))
            .addObject(0.0001, cm::Coord(-15, 0, 0), cm::Coord(0, -0.0005, 0));
    cm::Simulator simulator = cm::Simulator()
            .setSystem(twoObj);

    cm::ConsoleVisualizer()
            .setSimulator(simulator)
            .setViewport(40, 30)
            .run(argc, argv);

    return 0;
}
