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

// В принципе при желании можно ввести концепт для тех типов, которые можно приравнивать нулю
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

// Чтобы не писать всегда шаблонные параметры, определим type alias
using Vec3d = Vec3<double>;

class Fields {
    public:
    double ex, ey, ez, bx, by, bz;
    Fields(double, double, double, double, double, double);
    Fields(Vec3d e_field, Vec3d b_field);
};

Fields::Fields( double a, double b, double c
              , double d, double e, double f
              ) {
    ex = a;
    ey = b;
    ez = c;
    bx = d;
    by = e;
    bz = f;
}

Fields::Fields(Vec3d e_field, Vec3d b_field) {
    ex = e_field.x;
    ey = e_field.y;
    ez = e_field.z;
    bx = b_field.x;
    by = b_field.y;
    bz = b_field.z;
}
