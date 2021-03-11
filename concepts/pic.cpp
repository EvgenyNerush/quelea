// Для компиляции нужно включить поддержку концептов (доступна с GCC 8, осень 2019 года):
// g++ -fconcepts pic.cpp && ./a.out

// Почему бы вместо конфигов не писать каждый раз свой cpp файл?  Или вынести конфиг в config.hpp,
// если удастся сделать структуру pic.cpp общей для всех конфигураций? Это позволит без труда
// написать, например, конфиг с внешним полем (заданным через функцию) какого угодно вида, не
// занимаясь "велосипедостроением" по разбору конфиг-файла, созданием библиотек типичных огибающих
// и т. п.

// Можно сделать структуру pic.cpp довольно общей, но, скорее всего, при этом для того, чтобы она
// не зависила от сортов используемых частиц, придётся создать базовый класс - массив частиц одного
// сорта, и создать наследников этого класса (массив электронов со своим пушером, способом вывода в
// файл и т. п., массив ионов со своим способом размещения в памяти, пушером и т. д.), чего делать
// не хочется. Это просто усложнит код.

#include <iostream>
#include <functional>
#include <cmath>
#include "particles.hpp"

using namespace std;

double a0 = 2;
double dt = 0.1;
double t_end = 10;

Vec3<double> sine_field(double t, Vec3<double> r) {
    return Vec3<double>(a0 * sin(r.x - t), 0, 0);
}

Vec3<double> zero_f(double t, Vec3<double> r) {
    return zero<double>;
}

class Indexed_particle {
    public:
    double x, y, z, px, py, pz, g;
    long index;
};

auto pusherf = [](double t, Indexed_particle& p) {
    return Pusher3D3P< Pusher_type::Vay
                , Indexed_particle
                , double
                , double
                , double
                , std::function< Vec3<double>(double, Vec3<double>) >
                , std::function< Vec3<double>(double, Vec3<double>) > >::step(p, 1, dt, t, sine_field, zero_f);
};

int main() {
    Indexed_particle p;
    p.px = 0;
    p.py = 0;
    p.pz = 0;
    for (double t = 0; t < t_end; t += dt) {
        cout << p.px << '\t' << p.py << '\t' << p.pz << '\n';
        pusherf(t, p);
    }
    return 0;
}
