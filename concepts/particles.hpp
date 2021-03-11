/* Для компиляции нужно включить поддержку концептов:
 * $ g++ -fconcepts particles.hpp
 * Кстати, я понятия не имею, насколько хорошо концепты поддерживаются в doxygen, подозреваю, что
 * плохо.
 */

#include <functional>

// Vector2D - эта такой класс, у которого есть double члены x и y. Другие члены могут быть любыми
// (спин, метки, время рождения...)
template <typename T>
concept Vector2D = requires(T a) {
    { a.x } -> double;
    { a.y } -> double;
};

// Трёхмерный вектор - это двумерный вектор + ещё одна координата
template <typename T>
concept Vector3D = requires(T a) {
    { a }   -> Vector2D;
    { a.z } -> double;
};

// В простейшем случае продвижение координат не требует даже знания импульсов (как в
// асимптотической теории). В данном случае пушер - это функция большого числа аргументов (q, dt,
// глобальное t, E_x или даже функция E_x(x, t) и т. п., но последний аргумент - та частица,
// координаты которой мы продвигаем. См. pic.cpp для примеров конструкции пушеров частичной
// подстановкой параметров. Передавать частицу по не-константной ссылке может быть неудобно при
// вычисления токов, наводимых частицей; надо это проверить.
template <typename... Args>
void Pusher3D(Args... a, Vector3D& b) {
}
/* Сначала я хотел объявить пушер как концепт
template <typename T, typename P>
concept Pusher = requires(T a, P& b) {
    { a.type }                -> Pusher_type;
    { a.push_coordinates(b) } -> void;
    { a.push_momentum(b) }    -> void;
оставив возможности задать dt, E_x и т. п. в качестве остальных членов класса. Однако при таком
подходе для движения частицы в аналитически заданных полях естественным выглядит член класса t -
глобальное время (для передачи в функции поля E_x(x, t,..). Функция push_coordinates могла бы
продвигать t. Теперь вспомним, что наш пушер будет использоваться параллельно, что при таком
подходе сразу приведёт к его неправильной работе. Если же мы пушер в дальнейшем задаём как частично
применённую функцию, нам ничего не останется, как оставить t свободным параметром, "внешним" по
отношению к пушеру (и это правильно). Поэтому здесь variadic templates.
*/

// Можно отделить релятивизм (Лоренц-фактор) в отдельный концепт, но, мне кажется, пока это не
// принципиально
template <typename T>
concept Vector2P = requires(T a) {
    { a.px } -> double;
    { a.py } -> double;
    { a.g }  -> double; // Lorentz-factor
};

template <typename T>
concept Vector3P = requires(T a) {
    { a }    -> Vector2P;
    { a.pz } -> double;
};

// Пушер для импульса
template <typename... Args>
void Pusher3P(Args... a, Vector3P& b) {
}

// Небольшое лирическое отсупление.
// Обычный вектор, для векторных операций удобно его определить отдельно
template <typename T>
class Vec3 {
    public:
    T x;
    T y;
    T z;
    Vec3(T, T, T);
};

template <typename T>
Vec3<T>::Vec3(T a, T b, T c) {
    x = a;
    y = b;
    z = c;
}

// Немного duck typing-а; в принципе при желании можно ввести концепт для тех типов, которые можно
// приравнивать нулю
template <typename T>
Vec3<T> zero(0, 0, 0);

template <typename T>
Vec3<T> unit_x(1, 0, 0);

template <typename T>
Vec3<T> unit_y(0, 1, 0);

template <typename T>
Vec3<T> unit_z(0, 0, 1);

// Скалярное произведение
template <typename T>
T dot_product(Vec3<T> a, Vec3<T> b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Векторное произведение
template <typename T>
T cross_product(Vec3<T> a, Vec3<T> b) {
    Vec3<T> c;
    c.x =  a.y * b.z - a.z * b.y;
    c.y = -a.x * b.z + a.z * b.x;
    c.z =  a.x * b.y - a.y * b.x;
    return c;
}

// Обычная частица - содержит координаты (поэтому она Vector3D), импульсы и Лоренц-фактор (поэтому
// она одновременно и Vector3P). Аналогично можно определить Particle2D2P, 2D3P и т. п.
template <typename T>
concept Particle3D3P = requires(T a) {
    { a } -> Vector3D;
    { a } -> Vector3P;
};

// Мне кажется, можно попробовать записать в виде типов совместимость разных методов для решения
// уравнений движения и уравнений Максвелла. Но надо подумать, как бы это можно было бы сделать. В
// любом случае, нужен тип для индексации (то есть для параметризации, обозначения) пушеров (их
// ведь будет много).
enum class Pusher_type {
    Midpoint,
    Vay,
    Euler,
};

class qwe {
    public:
    double x, y, z, px, py, pz, g;
};

// Пушеры для частиц; из-за перегрузки функций для функций запрещено частичное определение (partial
// specialization) шаблонных параметров; в то же время для классов оно разрешена, чем и пользуемся
template <Pusher_type t, Particle3D3P P, typename... Args>
class Pusher3D3P {
    public:
    static void step(P&, Args... );
};

// Обычный midpoint-пушер для координат
template <Particle3D3P P>
class Pusher3D3P<Pusher_type::Midpoint, P, double> {
    public:
    static void step(P& a, double dt) {
        a.x += dt * a.px / a.g;
        a.y += dt * a.py / a.g;
        a.z += dt * a.pz / a.g;
    }
};

// Это не алгоритм Вэя, конечно, а просто заглушка, но она требует тех же параметров, что и Вэй
template <Particle3D3P P>
class Pusher3D3P< Pusher_type::Vay, P
                , double, double
                , double, double, double
                , double, double, double > {
    public:
    static void step( P& a
                    , double q, double dt
                    , double ex, double ey, double ez
                    , double bx, double by, double bz
                    ) {
        a.px += dt * q * ex;
        a.g = sqrt(1 + a.x * a.x + a.y * a.y + a.z * a.z);
    }
};

// Вэй для полей, заданных аналитически
template <Particle3D3P P>
class Pusher3D3P< Pusher_type::Vay
                , P
                , double
                , double
                , double
                , std::function< Vec3<double>(double, Vec3<double>) >
                , std::function< Vec3<double>(double, Vec3<double>) > > {
    public:
    static void step( P& a
                    , double q
                    , double dt
                    , double t
                    , std::function< Vec3<double>(double, Vec3<double>) > e_function
                    , std::function< Vec3<double>(double, Vec3<double>) > b_function
                    ) {
        Vec3<double> e_field = e_function(t, Vec3<double>(a.x, a.y, a.z));
        Vec3<double> b_field = e_function(t, Vec3<double>(a.x, a.y, a.z));
        Pusher3D3P< Pusher_type::Vay, P
                , double, double
                , double, double, double
                , double, double, double >::step(a, q, dt, e_field.x, e_field.y, e_field.z, b_field.x, b_field.y, b_field.z);
    }
};
