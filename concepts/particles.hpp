/* Для компиляции нужно включить поддержку концептов:
 * $ g++ -fconcepts particles.hpp
 * Кстати, я понятия не имею, насколько хорошо концепты поддерживаются в doxygen, подозреваю, что
 * плохо.
 */

#include <functional>
#include "vec.hpp"

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

// Обычная частица - содержит координаты (поэтому она Vector3D), импульсы и Лоренц-фактор (поэтому
// она одновременно и Vector3P). Аналогично можно определить Particle2D2P, 2D3P и т. п.
template <typename T>
concept Particle3D3P = requires(T a) {
    { a } -> Vector3D;
    { a } -> Vector3P;
};

// === ПУШЕРЫ === //
// В простейшем случае продвижение координат не требует даже знания импульсов (как в
// асимптотической теории). Но пока пушеры, которые используют только координаты и вообще не
// используют импульсы, вроде бы не нужны:
//
//     template <typename... Args>
//     void Pusher3D(Args... a, Vector3D& b) {
//     }
//
// В данном случае пушер - это функция большого числа аргументов (q, dt, глобальное t, E_x или даже
// в качестве аргумента может быть функция E_x(x, t) и т. п., но последний аргумент - та частица,
// координаты которой мы продвигаем.  См. pic.cpp для примеров конструкции пушеров частичной
// подстановкой параметров. Передавать частицу по не-константной ссылке может быть неудобно при
// вычисления токов, наводимых частицей; надо это проверить.
//
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
отношению к пушеру (и это правильно). Поэтому далее используются функции и variadic templates.
*/

// Мне кажется, можно попробовать записать в виде типов совместимость разных методов для решения
// уравнений движения и уравнений Максвелла. Но надо подумать, как бы это можно было бы сделать. В
// любом случае, нужен тип для индексации (то есть для параметризации, обозначения) пушеров (их
// ведь будет много).
enum class Pusher_type {
    Midpoint,
    Vay,
    Euler,
};

// Пушеры для частиц; (из-за перегрузки функций) для шаблонных функций запрещено частичное
// определение (partial specialization) шаблонных параметров; в то же время для классов оно
// разрешено, чем и пользуемся
template <Pusher_type t, Particle3D3P P, typename... Args>
class Pusher3D3P {
    public:
    static void step(P&, Args... );
};

// Обычный midpoint-пушер для координат
template <Particle3D3P P>
class Pusher3D3P< Pusher_type::Midpoint, P, double > {
    public:
    static void step(P& a, double dt) {
        a.x += dt * a.px / a.g;
        a.y += dt * a.py / a.g;
        a.z += dt * a.pz / a.g;
    }
};

// type alias, чтобы не нужно было всё время писать шаблонные параметры
template <Particle3D3P P>
using Midpoint_coord_pusher = Pusher3D3P<Pusher_type::Midpoint, P, double>;

// Это не алгоритм Вэя, конечно, а просто заглушка, но она требует тех же параметров, что и Вэй
template <Particle3D3P P>
class Pusher3D3P< Pusher_type::Vay, P , double, double, Fields > {
    public:
    static void step(P& a , double q, double dt , Fields f) {
        a.px += dt * q * f.ex;
    }
};

template <Particle3D3P P>
using Vay_pusher_simple = Pusher3D3P< Pusher_type::Vay, P, double, double, Fields >;

// Вэй для полей, заданных аналитически
template <Particle3D3P P>
class Pusher3D3P< Pusher_type::Vay
                , P
                , double
                , double
                , double
                , std::function< Fields(double, Vec3d) > > {
    public:
    static void step( P& a
                    , double q
                    , double dt
                    , double t
                    , std::function< Fields(double, Vec3d) > f
                    ) {
        Fields f_value = f(t, Vec3d(a.x, a.y, a.z));
        Vay_pusher_simple<P>::step(a, q, dt, f_value);
    }
};

template <Particle3D3P P>
using Vay_pusher_functional = Pusher3D3P < Pusher_type::Vay, P
                                         , double , double , double
                                         , std::function< Fields(double, Vec3d) > >;
