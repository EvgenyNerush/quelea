extern crate libc;
extern crate rand;
use libc::{c_int, c_double, size_t};

use std::f64;
use std::f64::consts;
use std::rc::Rc;
use rand::distributions::{IndependentSample, Range};

// нормировка в данной библиотеке - 1 / omega

// 3D vector
type V = (f64, f64, f64);

// EM fields
struct EM {
    ex: f64,
    ey: f64,
    ez: f64,
    bx: f64,
    by: f64,
    bz: f64,
}

// particle coordinates and momentum
#[repr(C)]
#[derive(Clone)]
struct Particle {
    x: f64,
    y: f64,
    z: f64,
    ux: f64,
    uy: f64,
    uz: f64,
}

//Constant homogeneous magnetic field, parallel to z axis with a0 amplitude
fn const_B(a0: f64) -> EM {
    EM {ex: 0.0, ey: 0.0, ez: 0.0, bx: 0.0, by: 0.0, bz: a0}
}

/* Constant homogeneous electric and magnetic field. If *theta* == 0, electric
field is parallel to y-axis, magnetic field is parallel to z-axis. E magnitude is a0 / 2, B
magnitude is a0.*/
fn const_EB(a0: f64, theta: f64) -> EM {
    EM {
        ex: 0.0,
        ey: 0.5 * a0 * f64::cos(theta),
        ez: 0.5 * a0 * f64::sin(theta),
        bx: 0.0,
        by: -a0 * f64::sin(theta),
        bz: a0 * f64::cos(theta),
    }
}

// Гармоническая плоская электромагнитная волна. Частота = 1, амплитуда электрического поля такова,
// что средняя интенсивность волны равна средней интенсивности линейно поляризованной волны
// амплитудой a0.
// Направление распространения волны задаётся волновым вектором *k*, эллипс поляризации
// определяется направлением вектора *E0* и "углом" *theta* так, что большая полуось эллипса
// поляризации лежит в одной плоскости с векторами *E0* и *k*, и, например, при k = (1, 0, 0) и E0
// = (0, 1, 0):
// E = a0 * (0, cos(theta) * cos(xi), sin(theta) * sin(xi)), где xi = x - t
fn plane_wave(r: (f64, f64, f64), t: f64, a0: f64, kx: f64, ky: f64, kz: f64, E0x: f64, E0y: f64,
              E0z: f64, theta: f64, omega: f64) -> EM {
    let k_times_E0 = (ky * E0z - kz * E0y,
                      kz * E0x - kx * E0z,
                      kx * E0y - ky * E0x);
    let k_times_E0_mod = f64::sqrt(k_times_E0.0 * k_times_E0.0 +
                                   k_times_E0.1 * k_times_E0.1 +
                                   k_times_E0.2 * k_times_E0.2);
    if k_times_E0_mod != 0.0 {
        let k_mod = f64::sqrt(kx * kx + ky * ky + kz * kz);
        let k = (kx / k_mod, ky / k_mod, kz / k_mod);
        // [k, E0] with norm = 1
        let kE0 = (k_times_E0.0 / k_times_E0_mod,
                   k_times_E0.1 / k_times_E0_mod,
                   k_times_E0.2 / k_times_E0_mod);
        // [[k, E0], k] with norm = 1
        let kE0k = (kE0.1 * kz - kE0.2 * ky,
                    kE0.2 * kx - kE0.0 * kz,
                    kE0.0 * ky - kE0.1 * kx);
        let sigma = 50.0;
        let xi = omega * (k.0 * r.0 + k.1 * r.1 + k.2 * r.2 - t) + sigma;
        let mut envelope = 0.0;
        if xi > sigma || xi < -sigma {
            envelope = 1.0;
            //EM{ex:0.0, ey:0.0, ez:0.0, bx:0.0, by:0.0, bz:0.0}
        } else {
            envelope = 1.0;  //f64::cos(0.5 * f64::consts::PI * (xi / sigma).powi(10)).powi(10);
        }
        let c = f64::cos(xi - sigma) * f64::cos(theta) * envelope;
        let s = f64::sin(xi - sigma) * f64::sin(theta) * envelope;
        EM { ex: a0 * (kE0k.0 * c + kE0.0 * s),
            ey: a0 * (kE0k.1 * c + kE0.1 * s),
            ez: a0 * (kE0k.2 * c + kE0.2 * s),
            bx: a0 * (kE0.0 * c - kE0k.0 * s),
            by: a0 * (kE0.1 * c - kE0k.1 * s),
            bz: a0 * (kE0.2 * c - kE0k.2 * s),
        }
    } else {
        panic!("*E0* and *k* are parallel to each other!");
    }
}


//Гауссов импульс с огибающей-косинусом
//E = a_0 \cos(\omega_0 t) \exp(-r^2/r_0^2)
pub const E: f64 = 2.71828182845904523536028747135266250f64;
fn gauss(r: (f64,f64, f64), t: f64, width: f64, a0: f64) -> EM {
    let x = r.0;
    //println!("{}", );
    let fields =
    EM {
        //let E0 = a0 * f64::cos(x -t)
        //let exponent = f64::pow(E, f64::powf(r/r0,2.0));
        ex: a0 * f64::cos(x - t),//* f64::powf(E, -1.0*f64::powf(x/width,2.0)),
        ey: 0.0,
        ez: 0.0,
        bx: 0.0,
        by: a0 * f64::cos(x - t), // f64::powf(E, -1.0*f64::powf(x/width,2.0)),
        bz: 0.0,
    };
    //println!("{}", fields.ex);
    return fields;
}
// Вращающееся электрическое поле для сравнения с аналитическим решением из Зельдовича
// Вращение происходит вокруг оси, задаваемой вектором k
// a0 - безразмерная амплитуда поля
fn rot_E(t: f64, a0: f64, kx: f64, ky: f64, kz: f64, E0x: f64, E0y: f64, E0z: f64) -> EM {
    let k_times_E0 = (ky * E0z - kz * E0y,
                      kz * E0x - kx * E0z,
                      kx * E0y - ky * E0x);
    let k_times_E0_mod = f64::sqrt(k_times_E0.0 * k_times_E0.0 +
                                   k_times_E0.1 * k_times_E0.1 +
                                   k_times_E0.2 * k_times_E0.2);
    if k_times_E0_mod != 0.0 {
        let k_mod = f64::sqrt(kx * kx + ky * ky + kz * kz);
        let k = (kx / k_mod, ky / k_mod, kz / k_mod);
        // [k, E0] with norm = 1
        let kE0 = (k_times_E0.0 / k_times_E0_mod,
                   k_times_E0.1 / k_times_E0_mod,
                   k_times_E0.2 / k_times_E0_mod);
        // [[k, E0], k] with norm = 1
        let kE0k = (kE0.1 * kz - kE0.2 * ky,
                    kE0.2 * kx - kE0.0 * kz,
                    kE0.0 * ky - kE0.1 * kx);
        let c = f64::cos(t);
        let s = f64::sin(t);
        EM { ex: a0 * (kE0k.0 * c + kE0.0 * s),
             ey: a0 * (kE0k.1 * c + kE0.1 * s),
             ez: a0 * (kE0k.2 * c + kE0.2 * s),
             bx: 0.0, by: 0.0,bz: 0.0}
    } else {
        panic!("*E0* and *k* are parallel to each other!");
    }
}

// E_y = a0, B_z = a_0 * x
fn const_E_lin_B(r: V, a0: f64) -> EM {
    EM {
        ex: 0.0,
        ey: a0,
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: a0 * r.0,
    }
}

// линейно поляризованная стоячая волна, a0 - максимальная амплитуда поля
fn lin_pol_stw(r: V, t: f64, a0: f64) -> EM {
    EM {
        ex: 0.0,
        ey: 0.5 * a0 * (f64::cos(r.0 - t) + f64::cos(r.0 + t)),
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: 0.5 * a0 * (f64::cos(r.0 - t) - f64::cos(r.0 + t)),
    }
}

// Две сходящихся линейно поляризованных волны с поглощающим слоем между ними. R - коэффициент
// "отражения" (он же коэффициент прохождения) волны от слоя. Волны в такой системе могут быть
// представлены как суперпозиция стоячей и бегущей волны. a0 - максимальное значение поля в стоячей
// волне, без "бегущей" части.
fn lin_pol_absorber(r: V, t: f64, a0: f64, R: f64) -> EM {
    if r.0 < 0.0 {
        return EM {
            ex: 0.0,
            ey: a0 * (f64::cos(r.0) * f64::cos(t) - 0.5 * (1.0 - R) * f64::cos(r.0 + t)),
            ez: 0.0,
            bx: 0.0,
            by: 0.0,
            bz: a0 * (f64::sin(r.0) * f64::sin(t) + 0.5 * (1.0 - R) * f64::cos(r.0 + t)),
        };
    } else if r.0 > 0.0 {
        return EM {
            ex: 0.0,
            ey: a0 * (f64::cos(r.0) * f64::cos(t) - 0.5 * (1.0 - R) * f64::cos(r.0 - t)),
            ez: 0.0,
            bx: 0.0,
            by: 0.0,
            bz: a0 * (f64::sin(r.0) * f64::sin(t) - 0.5 * (1.0 - R) * f64::cos(r.0 - t)),
        };
    } else {
        return EM {
            ex: 0.0,
            ey: a0 * 0.5 * (1.0 + R) * f64::cos(t),
            ez: 0.0,
            bx: 0.0,
            by: 0.0,
            bz: 0.0,
        };
    }
}

// TE11 мода волновода, периодически продлённая по y и z. _ly_, _lz_ --- размеры волновода.
// Максимум поля E_y равен a0. Стоит отметить, что если lz > ly, то max(E_z) > max(E_y) = a0.
// Предполагается, что k = 1.
fn te11(r: V, t: f64, a0: f64, ly: f64, lz: f64) -> EM {
    let ky = f64::consts::PI / ly;
    let kz = f64::consts::PI / lz;
    let kx2 = 1.0 - ky * ky - kz * kz;
    if kx2 < 0.0 {
        panic!("te11: закритический волновод!");
    }
    let kx = f64::sqrt(kx2);
    let ay = a0;
    let az = ay * ky / kz;
    // ay / kz = az / ky
    let alpha = ay / kz * (kz * kz + ky * ky);
    let (x, y, z) = r;
    let ex = 0.0;
    let ey =    ay * f64::cos(ky * y) * f64::sin(kz * z) * f64::cos(t - kx * x);
    let ez =   -az * f64::sin(ky * y) * f64::cos(kz * z) * f64::cos(t - kx * x);
    let bx = alpha * f64::cos(ky * y) * f64::cos(kz * z) * f64::sin(t - kx * x);
    let by =   -kx * ez;
    let bz =    kx * ey;
    EM { ex: ex, ey: ey, ez: ez, bx: bx, by: by, bz: bz }
}

/* Сумма плоской гармонической волны, падающей из области (x < 0, y < 0), с её зеркальным
 * отражением, так что: 1) E_y(x = 0) = 0, 2) угол между осью x и k есть theta, 3) k = 1. */
fn dms_field(r: V, t: f64, a0: f64, theta: f64) -> EM {
    let kx = f64::cos(theta);
    let ky = f64::sin(theta);
    let (x, y, z) = r;
    let El = a0 * f64::cos(t - kx * x - ky * y);
    let Er = -a0 * f64::cos(t + kx * x - ky * y);
    EM {
        ex: (-El + Er) * f64::sin(theta), // -2 a0 sin(theta) * cos(t - ky * y) * cos(kx * x)
        ey: (El + Er) * f64::cos(theta), // 2 a0 cos(theta) * sin(t - ky * y) * sin(kx * x)
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: El - Er, // -ex / sin(theta)
    }
}

/* Периодическое (во времени) электромагнитное поле, задаваемое потенциалом A_phi(r, z) (phi - угол
 * цилиндрической системы координат, относительно вращения по phi потенциал симметричен; потенциал
 * также считается симметричным (чётным) относительно преобразовния z -> -z).
 * Потенциал задаётся на прямоугольной сетке (nr, nz), соответствующей области [0, r_m] по r и
 * [0, z_m] по z (крайние точки лежат на границе области), через матрицу мнимой и вещественной
 * своей части, без exp(-i \omega t). Поля вычисляются через численное дифференцирование
 * потенциала. Матрицы для представления в виде одномерных массивов разбиваются на строки (длиной
 * nr), а не на столбцы. omega = 1 */
fn fields_from_A_phi(r: (f64, f64, f64), t: f64, n: (usize, usize), limits: (f64, f64), A_phi:
                     &(Vec<f64>, Vec<f64>)) -> EM {
    let (nr, nz) = n;
    let (r_m, z_m) = limits;
    let (ref re_A, ref im_A) = *A_phi;
    let (x, y, z0) = r;
    let z = f64::abs(z0);
    let r = f64::sqrt(x * x + y * y);
    let dr = r_m / (nr - 1) as f64;
    let dz = z_m / (nz - 1) as f64;
    let i: usize = f64::floor(r  / dr) as usize;
    let j: usize = f64::floor(f64::abs(z) / dz) as usize;
    // сетка для B_z смещена на -0.5 dr вдоль r, для B_r - на -0.5 dz вдоль z.
    let k: usize = f64::floor(r / dr + 0.5) as usize;
    let l: usize = f64::floor(z / dz + 0.5) as usize;
    if k < nr - 1usize && l < nz - 1usize { // на верхних пределах нет точек для линейной интерполяции
        let re_exp = f64::cos(t);
        let im_exp = -f64::sin(t);
        // zb 01 -- 11
        // :  |      |
        // Za 00 -- 10
        //     ra..rb
        let ra = r - (i as f64) * dr;
        let rb = 1.0 - ra;
        let za = z - (j as f64) * dz;
        let zb = 1.0 - za;
        let re_A_00 = re_A[i + j * nr];
        let re_A_10 = re_A[i + 1usize + j * nr];
        let re_A_01 = re_A[i + (j + 1usize) * nr];
        let re_A_11 = re_A[i + 1usize + (j + 1usize) * nr];
        let im_A_00 = im_A[i + j * nr];
        let im_A_10 = im_A[i + 1usize + j * nr];
        let im_A_01 = im_A[i + (j + 1usize) * nr];
        let im_A_11 = im_A[i + 1usize + (j + 1usize) * nr];
        // E_phi = i k A_phi_0 e^{-i \omega_t}
        let e00 = -(im_A_00 * re_exp + re_A_00 * im_exp);
        let e01 = -(im_A_01 * re_exp + re_A_01 * im_exp);
        let e10 = -(im_A_10 * re_exp + re_A_10 * im_exp);
        let e11 = -(im_A_11 * re_exp + re_A_11 * im_exp);
        let E_phi = e00 * rb * zb + e10 * ra * zb + e01 * rb * za + e11 * ra * za;
        // B_r = -\partial A_\phi / \partial z,
        // B_z = \partial A_\phi / \partial r.
        let (re_dAdz_00, re_dAdz_10, re_dAdz_01, re_dAdz_11,
             im_dAdz_00, im_dAdz_10, im_dAdz_01, im_dAdz_11) = 
            if l != 0 {
                ((re_A[i + l * nr] - re_A[i + (l - 1usize) * nr]) / dz,
                 (re_A[i + 1usize + l * nr] - re_A[i + 1usize + (l - 1usize) * nr]) / dz,
                 (re_A[i + (l + 1usize) * nr] - re_A[i + l * nr]) / dz,
                 (re_A[i + 1usize + (l + 1usize) * nr] - re_A[i + 1usize + l * nr]) / dz,
                 (im_A[i + l * nr] - im_A[i + (l - 1usize) * nr]) / dz,
                 (im_A[i + 1usize + l * nr] - im_A[i + 1usize + (l - 1usize) * nr]) / dz,
                 (im_A[i + (l + 1usize) * nr] - im_A[i + l * nr]) / dz,
                 (im_A[i + 1usize + (l + 1usize) * nr] - im_A[i + 1usize + l * nr]) / dz,
                )
            } else {
                ((re_A[i] - re_A[i + nr]) / dz,
                 (re_A[i + 1usize] - re_A[i + 1usize + nr]) / dz,
                 (re_A[i + nr] - re_A[i]) / dz,
                 (re_A[i + 1usize + nr] - re_A[i + 1usize]) / dz,
                 (im_A[i] - im_A[i + nr]) / dz,
                 (im_A[i + 1usize] - im_A[i + 1usize + nr]) / dz,
                 (im_A[i + nr] - im_A[i]) / dz,
                 (im_A[i + 1usize + nr] - im_A[i + 1usize]) / dz,
                )
            };
        let (re_dAdr_00, re_dAdr_10, re_dAdr_01, re_dAdr_11,
             im_dAdr_00, im_dAdr_10, im_dAdr_01, im_dAdr_11) = 
            if k != 0 {
                ((re_A[k + nr * j] - re_A[k - 1usize + nr * j]) / dr,
                 (re_A[k + 1usize + nr * j] - re_A[k + nr * j]) / dr,
                 (re_A[k + (j + 1usize) * nr] - re_A[k - 1usize + (j + 1usize) * nr]) / dr,
                 (re_A[k + 1usize + (j + 1usize) * nr] - re_A[k + (j + 1usize) * nr]) / dr,
                 (im_A[k + nr * j] - im_A[k - 1usize + nr * j]) / dr,
                 (im_A[k + 1usize + nr * j] - im_A[k + nr * j]) / dr,
                 (im_A[k + (j + 1usize) * nr] - im_A[k - 1usize + (j + 1usize) * nr]) / dr,
                 (im_A[k + 1usize + (j + 1usize) * nr] - im_A[k + (j + 1usize) * nr]) / dr,
                )
            } else {
                ((re_A[nr * j] + re_A[1usize + nr * j]) / dr,
                 (re_A[1usize + nr * j] - re_A[nr * j]) / dr,
                 (re_A[(j + 1usize) * nr] + re_A[1usize + (j + 1usize) * nr]) / dr,
                 (re_A[1usize + (j + 1usize) * nr] - re_A[(j + 1usize) * nr]) / dr,
                 (im_A[nr * j] + im_A[1usize + nr * j]) / dr,
                 (im_A[1usize + nr * j] - im_A[nr * j]) / dr,
                 (im_A[(j + 1usize) * nr] + im_A[1usize + (j + 1usize) * nr]) / dr,
                 (im_A[1usize + (j + 1usize) * nr] - im_A[(j + 1usize) * nr]) / dr,
                )
            };
        // B_r
        let za = z + 0.5 - (l as f64) * dz;
        let zb = 1.0 - za;
        let B_r = -re_exp * (re_dAdz_00 * rb * zb + re_dAdz_01 * ra * zb + re_dAdz_10 * rb * za + re_dAdz_11 * ra * za) +
                   im_exp * (im_dAdz_00 * rb * zb + im_dAdz_01 * ra * zb + im_dAdz_10 * rb * za + im_dAdz_11 * ra * za);
        // B_z
        let ra = r + 0.5 - (k as f64) * dr;
        let rb = 1.0 - ra;
        let za = z - (j as f64) * dz;
        let zb = 1.0 - za;
        let B_z = re_exp * (re_dAdr_00 * rb * zb + re_dAdr_01 * ra * zb + re_dAdr_10 * rb * za + re_dAdr_11 * ra * za) -
                  im_exp * (im_dAdr_00 * rb * zb + im_dAdr_01 * ra * zb + im_dAdr_10 * rb * za + im_dAdr_11 * ra * za);
        let em = if r != 0.0 {
        // При замене z -> -z меняется только B_r -> -B_r; B_z не меняется
            EM {
                ex: -E_phi * y / r,
                ey: E_phi * x / r,
                ez: 0.0,
                bx: f64::signum(z0) * B_r * x / r,
                by: f64::signum(z0) * B_r * y / r,
                bz: B_z,
            }
        } else {
            EM {
                ex: 0.0,
                ey: 0.0,
                ez: 0.0,
                bx: 0.0,
                by: 0.0,
                bz: B_z,
            }
        };
        em
    } else {
        EM {
            ex: 0.0,
            ey: 0.0,
            ez: 0.0,
            bx: 0.0,
            by: 0.0,
            bz: 0.0,
        }
    }
}

/* Тестовый потенциал для использования в fields_from_A_phi, A_phi = (1 - i) * r * cos z,
 * приводящий к E_phi = r cos(z) (sin t + cos t), B_z = cos(z) (cos t - sin t),
 * B_r = r sin(z) (cos t - sin t) */
fn test_A_phi(n: (usize, usize), limits: (f64, f64)) -> (Vec<f64>, Vec<f64>) {
    let (nr, nz) = n;
    let (r_m, z_m) = limits;
    let dr = r_m / (nr - 1) as f64;
    let dz = z_m / (nz - 1) as f64;
    let num_of_elements = nr * nz;
    let mut re_A = Vec::with_capacity(num_of_elements);
    let mut im_A = Vec::with_capacity(num_of_elements);
    for k in 0usize..num_of_elements {
        let i = k % nr;
        let j = k / nr;
        let r = i as f64 * dr;
        let z = j as f64 * dz;
        re_A.push(r * f64::cos(z));
        im_A.push(-r * f64::cos(z));
    }
    (re_A, im_A)
}

/* Потенциал (без exp(-i t)) витка с с переменным током. Амплитуда тока, текущего по витку = J0,
 * радиус витка = barlambda. Возвращает матрицу для дальнейшего использования в fields_from_A_phi.
 * n_int --- число отрезков для взятия интеграла.
 * */
fn current_loop_A_phi(n: (usize, usize), limits: (f64, f64), J0: f64, barlambda: f64, n_int: i32) -> (Vec<f64>, Vec<f64>) {
    let (nr, nz) = n;
    let (r_m, z_m) = limits;
    let dr = r_m / (nr - 1) as f64;
    let dz = z_m / (nz - 1) as f64;
    let num_of_elements = nr * nz;
    let mut re_A = Vec::with_capacity(num_of_elements);
    let mut im_A = Vec::with_capacity(num_of_elements);
    for k in 0usize..num_of_elements {
        let i = k % nr;
        if i != 0 {
            let j = k / nr;
            let r = i as f64 * dr;
            let z = j as f64 * dz;
            //
            let barr = r + barlambda;
            let barz2 = z * z + barr * barr;
            let kappa = f64::sqrt(4.0 * r * barlambda / barz2);
            // в нормировке данной библиотеки в уравнениях Максвелла $4 pi / c$ перед j исчезает
            let ampl = 2.0 * J0 * barlambda / (consts::PI * kappa * kappa * f64::sqrt(barz2));
            // Интегрирование методом трапеций (прямоугольников в средней точке); при этом мы избегаем
            // проблемы с theta = 0, хотя подинтегральная функция имеет особенность и точность такого
            // метода может быть не очень большой.
            let mut re_int = 0f64;
            let mut im_int = 0f64;
            let dtheta = 0.5 * consts::PI / n_int as f64;
            for l in 0..n_int {
                let theta = 0.5 * consts::PI * (l as f64 + 0.5) / n_int as f64;
                let sq = f64::sqrt(1.0 - kappa * kappa * f64::sin(theta) * f64::sin(theta));
                let kR = f64::sqrt(barz2) * sq;
                let in_braces = (1.0 - 0.5 * kappa * kappa) / sq - sq;
                re_int += f64::cos(kR) * in_braces;
                im_int += f64::sin(kR) * in_braces;
            }
            re_int *= ampl * dtheta;
            im_int *= ampl * dtheta;
            //
            re_A.push(re_int);
            im_A.push(im_int);
        } else {
            // точно на оси z потенциал равен нулю из-за симметрии тока
            re_A.push(0.0);
            im_A.push(0.0);
        }
    }
    (re_A, im_A)
}

// Плоская стоячая волна. Задаётся аналогично plane_wave
fn standing_wave(r: (f64, f64, f64), t: f64, a0: f64, kx: f64, ky: f64, kz: f64, E0x: f64, E0y: f64,
              E0z: f64, theta: f64) -> EM {
    let k_times_E0 = (ky * E0z - kz * E0y,
                      kz * E0x - kx * E0z,
                      kx * E0y - ky * E0x);
    let k_times_E0_mod = f64::sqrt(k_times_E0.0 * k_times_E0.0 +
                                   k_times_E0.1 * k_times_E0.1 +
                                   k_times_E0.2 * k_times_E0.2);
    if k_times_E0_mod != 0.0 {
        let k_mod = f64::sqrt(kx * kx + ky * ky + kz * kz);
        let k = (kx / k_mod, ky / k_mod, kz / k_mod);
        // [k, E0] with norm = 1
        let kE0 = (k_times_E0.0 / k_times_E0_mod,
                   k_times_E0.1 / k_times_E0_mod,
                   k_times_E0.2 / k_times_E0_mod);
        // [[k, E0], k] with norm = 1
        let kE0k = (kE0.1 * kz - kE0.2 * ky,
                    kE0.2 * kx - kE0.0 * kz,
                    kE0.0 * ky - kE0.1 * kx);
        let x = k.0 * r.0 + k.1 * r.1 + k.2 * r.2;
        let cE = 2.0 * f64::cos(x) * f64::cos(t) * f64::cos(theta);
        let sE = 2.0 * f64::cos(x) * f64::sin(t) * f64::sin(theta);
        let cB = 2.0 * f64::sin(x) * f64::cos(t) * f64::cos(theta);
        let sB = 2.0 * f64::sin(x) * f64::sin(t) * f64::sin(theta);
        EM { ex: a0 * (kE0k.0 * cE + kE0.0 * sE),
             ey: a0 * (kE0k.1 * cE + kE0.1 * sE),
             ez: a0 * (kE0k.2 * cE + kE0.2 * sE),
             bx: a0 * (kE0.0 * cB - kE0k.0 * sB),
             by: a0 * (kE0.1 * cB - kE0k.1 * sB),
             bz: a0 * (kE0.2 * cB - kE0k.2 * sB),
        }
    } else {
        panic!("*E0* and *k* are parallel to each other!");
    }
}

// Плоская стоячая циркулярно поляризованная волна. E и B - в плоскости YZ
fn cir_pol_stw(r: (f64, f64, f64), t: f64, a0: f64) -> EM {
    let x = r.0;
    let cosx = f64::cos(x);
    let sinx = f64::sin(x);
    let cost = f64::cos(t);
    let sint = f64::sin(t);
    EM {ex: 0.0,
        ey: a0 * cosx * sint,
        ez: a0 * cosx * cost,
        bx: 0.0,
        by: -a0 * sinx * sint,
        bz: -a0 * sinx * cost}
}

// Сумма двух плоских линейно поляризованных волн
// первая: на основной частоте, сонаправленые Е. противоположно направленные В, амплитуда а0
// Вторая: на частоте harmonic, сонаправленые В. противоположно направленные Е, амплитуда а1
fn double_standing_wave(r: (f64, f64, f64), t:f64, a0: f64, a1: f64, phase0: f64, phase1: f64, harmonic: f64) -> EM {
    let x = r.0;
    /*EM {ex: 0.0,
        ey: a0 * (f64::cos(x - t + phase0) + f64::cos(x + t + phase0)) + a1 * (f64::cos(harmonic * (x - t) + phase1) - f64::cos(harmonic * (x + t) + phase1)),
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: a0 * (f64::cos(x - t + phase0) - f64::cos(x + t + phase0)) + a1 * (f64::cos(harmonic * (x - t) + phase1) + f64::cos(harmonic * (x + t) + phase1)),
        }
        */
    EM {ex: 0.0,
        ey: a0 * f64::cos(x - t + phase0) + a1 * f64::cos(harmonic * (x + t) + phase1),
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: a0 * f64::cos(x - t + phase0) - a1 * f64::cos(harmonic * (x + t) + phase1),
        }
}

// аналог f_init_cos из QUILL, но без поперечной огибающей; при t = 0 центр - в нуле,
// распространяется вдоль x.
fn cosine_pulse(r: (f64, f64, f64), t:f64, a0: f64, xs: f64) -> EM {
    let xi = r.0 - t;
    let cosx = f64::cos(0.5 * f64::consts::PI * xi / xs);
    let sinx = f64::sin(0.5 * f64::consts::PI * xi / xs);
    let ey = a0 * (f64::cos(xi) * cosx * cosx - f64::consts::PI / xs * f64::sin(xi) * cosx * sinx);
    EM {ex: 0.0,
        ey: ey,
        ez: 0.0,
        bx: 0.0,
        by: 0.0,
        bz: ey,
        }
}

fn calc_envelope(r: [f64; 3], orts: [[f64; 3]; 3], etype: i32, center: [f64; 3], size: [f64; 3]) -> f64 {
    if etype != 1 && etype != 2 {
        1.0
    }
    else {
        let mut x = r[0] - center[0];
        let mut y = r[1] - center[1];
        let mut z = r[2] - center[2];
        x = x * orts[0][0] + y * orts[0][1] + z * orts[0][2];
        y = x * orts[1][0] + y * orts[1][1] + z * orts[1][2]; //projections to orts
        z = x * orts[2][0] + y * orts[2][1] + z * orts[2][2];
        if etype == 1 {
            if x < 0.5 * size[0] && x > - 0.5 * size[0] &&
               y < 0.5 * size[1] && y > - 0.5 * size[1] &&
               z < 0.5 * size[2] && z > - 0.5 * size[2] {
            return 1.0;
            }
            else {
                return 0.0;
            }
        }

        if etype == 2 {
            if x < 0.5 * size[0] && x > - 0.5 * size[0] &&
               y < 0.5 * size[1] && y > - 0.5 * size[2] &&
               z < 0.5 * size[2] && y > - 0.5 * size[2] {
                return f64::cos(0.5 * consts::PI * (2.0 * x / size[0]).powi(4)).powi(2) * 
                       f64::cos(0.5 * consts::PI * (2.0 * y / size[1]).powi(4)).powi(2) *  
                       f64::cos(0.5 * consts::PI * (2.0 * z / size[2]).powi(4)).powi(2);
            }
            else {
                return 0.0;
            }
        }
        else {
            return 0.0;
        }
    }
}

// *s* is the initial coordinates and momentum of a particle (at t = 0), *f* is a function of r and
// t that returns value of EM fields, *dt* is the time step, *n* is the number of steps. Note that
// len of returned Vec is n + 1.
// Vay's algorithm is used (Phys. Plasmas 2008)
fn drive_particle_new(s0: Particle, f: Box<Fn(V, f64) -> EM>, dt: f64, n: usize) -> Vec<Particle> {
    // issue: should Rust delete[] s0?
    let mut vec: Vec<Particle> = Vec::new();
    // для начала сделаем пол шага назад по координатам, поскольку в алгоритме Вэя координаты и
    // импульсы разнесены по времени
    let g0 = f64::sqrt(1.0 + s0.ux * s0.ux + s0.uy * s0.uy + s0.uz * s0.uz);
    let r0 = (s0.x - 0.5 * dt * s0.ux / g0,
              s0.y - 0.5 * dt * s0.uy / g0,
              s0.z - 0.5 * dt * s0.uz / g0);
    // теперь делаем n + 1 шагов
    let mut s = s0;
    s.x = r0.0;
    s.y = r0.1;
    s.z = r0.2;
    vec.push(s.clone());
    let mut _1_g = 1.0 / g0;
    for i in 0..(n + 1) {
        // r step
        s.x += dt * s.ux * _1_g;
        s.y += dt * s.uy * _1_g;
        s.z += dt * s.uz * _1_g;
        // u step
        let em = f((s.x, s.y, s.z), (i as f64 + 0.5) * dt);
        let e = (-em.ex * dt, -em.ey * dt, -em.ez * dt);
        let dt_2 = 0.5 * dt;
        let b = (-em.bx * dt_2, -em.by * dt_2, -em.bz * dt_2);
        let w = (s.ux + e.0 + _1_g * (s.uy * b.2 - s.uz * b.1),
                 s.uy + e.1 + _1_g * (s.uz * b.0 - s.ux * b.2),
                 s.uz + e.2 + _1_g * (s.ux * b.1 - s.uy * b.0));
        let mut bb = b.0 * b.0 + b.1 * b.1 + b.2 * b.2;
        let mut bw = b.0 * w.0 + b.1 * w.1 + b.2 * w.2;
        let h = 1.0 + w.0 * w.0 + w.1 * w.1 + w.2 * w.2 - bb;
        _1_g = 1.0 / f64::sqrt(0.5 * (h + f64::sqrt(h * h + 4.0 * (bb + bw * bw)))); // at {i+1}
        bb = 1.0 / (1.0 + bb * _1_g * _1_g);
        bw *= _1_g * _1_g;
        s.ux = bb * (w.0 + (w.1 * b.2 - w.2 * b.1) * _1_g + b.0 * bw);
        s.uy = bb * (w.1 + (w.2 * b.0 - w.0 * b.2) * _1_g + b.1 * bw);
        s.uz = bb * (w.2 + (w.0 * b.1 - w.1 * b.0) * _1_g + b.2 * bw);
        vec.push(s.clone());
    }
    // теперь усредняем (x, y, z), чтобы получить значения в те же моменты времени, что и (ux, uy,
    // uz)
    for i in 0..(vec.len() - 1) {
        vec[i].x = 0.5 * (vec[i].x + vec[i + 1].x);
        vec[i].y = 0.5 * (vec[i].y + vec[i + 1].y);
        vec[i].z = 0.5 * (vec[i].z + vec[i + 1].z);
    }
    vec.pop();
    vec
}

fn W(g: f64, r: f64, chi: f64, coef: f64) -> f64 {
    let alpha = 1.0 / 137.0;
    if (1.0 - r) * g < 1.0 || chi == 0.0 || r == 0.0 {
        return 0.0;
    } 
    else if chi > 0.13333 {
        let d = r / ((1.0 - r) * chi);
        let a = d.powf(1.0 / 3.0);
        let kp = a * a;
        let pi = f64::consts::PI;
        return -(alpha / (2.0 * f64::sqrt(pi))) * (coef / (g * g)) * f64::exp(- 2.0 / 3.0 * d) * ((kp + 0.80049).powf(-0.75) - chi * a * (2.0 * (1.0 - r) + r * r) / r * (kp + 0.70861 * (1.0 - 0.65 * kp / (1.0 + kp * kp))).powf(0.25));
    }
    else {
        let rm = 1.0 / (1.0 + 2.0 / (15.0 * chi));
        let r = r * rm;
        let d = r / ((1.0 - r) * chi);
        let a = d.powf(1.0 / 3.0);
        let kp = a * a;
        return -rm * (alpha / (2.0 * f64::sqrt(f64::consts::PI))) * (coef / (g * g)) * f64::exp(- 2.0 / 3.0 * d) * ((kp + 0.80049).powf(-0.75) - chi * a * (2.0 * (1.0 - r) + r * r) / r * (kp + 0.70861 * (1.0 - 0.65 * kp / (1.0 + kp * kp))).powf(0.25));
    }
}

fn w_test(g: f64, r: f64, chi: f64, coef: f64) -> f64 {
    // тестовое распределение вероятности d W / d gamma,
    // такое, что W = M / tau.
    let r0 = 0.3;
    let dr = 0.05;
    let tau = 22.0;
    let M = 0.5;
    let beta = if r > r0 - 0.5 * dr && r < r0 + 0.5 * dr { M / (g * dr * tau) } else { 0.0 };
    beta
}

fn W_new(g: f64, r: f64, chi: f64, omega: f64) -> f64 {
    if (1.0 - r) * g < 1.0 || chi == 0.0 || r == 0.0 {
        return 0.0;
    }
    else if chi > 0.13333 { // 0.13(3) corresponds to a cut off at 5 omega_c
        let alpha = 1.0 / 137.036; // e^2 / hbar c
        let coef = 9.1e-28 * 3e10 * 3e10 / (1.05e-27 * omega); // mc^2 / hbar omega
        let y = r / ((1.0 - r) * chi); // hbar * omega / (chi * (gamma * mc^2 - hbar * omega))
        let x = f64::powf(y, 2.0 / 3.0);
        let exp = f64::exp(-2.0 * y / 3.0);
        let iAiryapp = exp / (2.0 * f64::sqrt(f64::consts::PI)) * f64::powf(x + 0.80049, -0.75);
        let dAiryapp = - exp / (2.0 * f64::sqrt(f64::consts::PI)) * f64::powf((x + 0.70861 * (1.0 - 0.65 * x / (1.0 + x * x))), 0.25);
        return - alpha * coef / (g * g) * (iAiryapp + (2.0 / x + chi * r * f64::sqrt(x)) * dAiryapp);
    }
    else {
        //let rm = 1.0 / (1.0 + 2.0 / (15.0 * chi));
        let rm = chi / 0.13333;
        //let rm = 1.0;
        let r = r * rm;
        let alpha = 1.0 / 137.036; // e^2 / hbar c
        let coef = 9.1e-28 * 3e10 * 3e10 / (1.05e-27 * omega); // mc^2 / hbar omega
        let y = r / ((1.0 - r) * chi); // hbar * omega / (chi * (gamma * mc^2 - hbar * omega))
        let x = f64::powf(y, 2.0 / 3.0);
        let exp = f64::exp(-2.0 * y / 3.0);
        let iAiryapp = exp / (2.0 * f64::sqrt(f64::consts::PI)) * f64::powf(x + 0.80049, -0.75);
        let dAiryapp = - exp / (2.0 * f64::sqrt(f64::consts::PI)) * f64::powf((x + 0.70861 * (1.0 - 0.65 * x / (1.0 + x * x))), 0.25);
        return - rm * alpha * coef / (g * g) * (iAiryapp + (2.0 / x + chi * r * f64::sqrt(x)) * dAiryapp);
    }
}

// Drives particle as in drive_particle_new but includes radiation friction; lambda задаётся в
// микронах.
fn drive_particle(p0: Particle, f: Box<Fn(V, f64) -> EM>, t0: f64, dt: f64, n: usize, lambda: f64, rf_scheme: i32) -> Vec<Particle> {
    let range = Range::new(0.0, 1.0);
    let mut rng = rand::thread_rng();
    let mut Num: Vec<Particle> = Vec::new();
    let mut p = p0;
    let mut g: f64 = f64::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
    let mut chi = 0.0;
    //Пол шага назад
    let mut inv_g = 1.0 / g;
    p.x -= 0.5 * dt * p.ux * inv_g;
    p.y -= 0.5 * dt * p.uy * inv_g;
    p.z -= 0.5 * dt * p.uz * inv_g;
    Num.push(p.clone());

    // coef = \frac{2}{3} \frac{e^2}{\hbar c} \frac{\hbar \omega}{m c^2}
    let coef = (2.0 / 3.0) / 137.036 * 1.0546e-27 * (2.0 * consts::PI * 3e10 / (1e-4 * lambda)) /
        (9.109e-28 * 3e10 * 3e10);
    let omega = 3e10 * f64::consts::PI * 2.0 / (lambda * 1e-4);
    let q_coef = (9.109e-28 * 3e10 * 3e10) / (omega * 1.054e-27); // mc^2 / hbar * omega
    for i in 0..(n + 1) {
        //Vay relativistic particle pusher from PhysPlasmas 2008
        p.x += dt * p.ux * inv_g;
        p.y += dt * p.uy * inv_g;
        p.z += dt * p.uz * inv_g;
        let em = f((p.x, p.y, p.z), t0 + (i as f64 + 0.5) * dt);
        // найдём силу Лоренца, нужную для вычисления силы трения, здесь. Раз уж для F_rr
        // используется метод Эйлера, то нет смысла использовать "новый" импульс и т. п.
        let lor_f = [-em.ex - (p.uy * em.bz - p.uz * em.by) * inv_g, 
                     -em.ey - (p.uz * em.bx - p.ux * em.bz) * inv_g,
                     -em.ez - (p.ux * em.by - p.uy * em.bx) * inv_g];
        g = 1.0 / inv_g;
        let frr = coef * dt * (
            g * g * (lor_f[0] * lor_f[0] + lor_f[1] * lor_f[1] + lor_f[2] * lor_f[2]) -
            (p.ux * lor_f[0] + p.uy * lor_f[1] + p.uz * lor_f[2]).powi(2)
            );
        let rr = (frr * p.ux * inv_g, frr * p.uy * inv_g, frr * p.uz * inv_g);
        // Vay's pusher, PoP, 2008
        let e = [-em.ex * dt, -em.ey * dt, -em.ez * dt];
        let b = [-em.bx * 0.5 * dt, -em.by * 0.5 * dt, -em.bz * 0.5 * dt];
        let wx = p.ux + e[0] + (p.uy * b[2] - p.uz * b[1]) * inv_g;
        let wy = p.uy + e[1] + (p.uz * b[0] - p.ux * b[2]) * inv_g;
        let wz = p.uz + e[2] + (p.ux * b[1] - p.uy * b[0]) * inv_g;
        let mut tmp = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        let mut bw = b[0] * wx + b[1] * wy + b[2] * wz;
        g = 1.0 + wx * wx + wy * wy + wz * wz - tmp; //sigma из статьи
        inv_g = 1.0 / f64::sqrt(0.5 * (g + f64::sqrt(g * g + 4.0 * (tmp + bw * bw)))); //gamma i+1
        g = 1.0 / inv_g;
        tmp = 1.0 / (1.0 + tmp * inv_g * inv_g); //s из статьи
        bw = bw * inv_g * inv_g;
        let ux = (wx + (wy * b[2] - wz * b[1]) * inv_g + b[0] * bw) * tmp;
        let uy = (wy + (wz * b[0] - wx * b[2]) * inv_g + b[1] * bw) * tmp;
        let uz = (wz + (wx * b[1] - wy * b[0]) * inv_g + b[2] * bw) * tmp;

        if rf_scheme == 0 {
            //Classic radiation reaction
            p.ux = ux - rr.0;
            p.uy = uy - rr.1;
            p.uz = uz - rr.2;
            inv_g += frr * inv_g * inv_g; // d(1 / gamma) / dt = F_RR  /gamma^2
            Num.push(p.clone());
        } else {
            //Quantum radiation friction
            let g_tmp = 0.5 * (f64::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz) + g); //g_i+1/2
            let u_tmp = [0.5 * (ux + p.ux),
                         0.5 * (uy + p.uy), //u_i+1/2
                         0.5 * (uz + p.uz)];
            let lor_f2 = [-em.ex - (u_tmp[1] * em.bz - u_tmp[2] * em.by) / g_tmp, 
                          -em.ey - (u_tmp[2] * em.bx - u_tmp[0] * em.bz) / g_tmp,
                          -em.ez - (u_tmp[0] * em.by - u_tmp[1] * em.bx) / g_tmp];
            let frr2 = g_tmp * g_tmp * (lor_f2[0] * lor_f2[0] + lor_f2[1] * lor_f2[1] + lor_f2[2] * lor_f2[2]) - (u_tmp[0] * lor_f2[0] + u_tmp[1] * lor_f2[1] + u_tmp[2] * lor_f2[2]).powi(2);
            chi = (1.0 / q_coef) * f64::sqrt(frr2);
            let mut r = range.ind_sample(&mut rng);
            let w = W_new(g_tmp, r, chi, omega);
            /*if chi < 0.13333 {
                r = r * chi / 0.13333;
            }*/
            if range.ind_sample(&mut rng) < w * dt * g_tmp {
                if chi < 0.13333 {
                    let rm = chi / 0.13333;
                    r = r * rm;
                }
                p.ux = ux - r * u_tmp[0];
                p.uy = uy - r * u_tmp[1];
                p.uz = uz - r * u_tmp[2];
                inv_g = 1.0 / f64::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
            } else {
                p.ux = ux;
                p.uy = uy;
                p.uz = uz;
            }
            Num.push(p.clone());
        }
    }
    //усреднение координат i и i+1, чтобы получить i+1/2 
    for i in 0..n {
        Num[i].x = 0.5 * (Num[i].x + Num[i + 1].x);
        Num[i].y = 0.5 * (Num[i].y + Num[i + 1].y);
        Num[i].z = 0.5 * (Num[i].z + Num[i + 1].z);
    }
    Num.pop();
    Num
}

// calls a function *f* with list of parameters x1, x2, ..., z[0], z[1], ...
// Examples:
// call!(f()() _) yields f()
// call!(f(x,y)() _) yields f(x,y)
// call!(f()(_) vec![137]) yields f(137)
// call!(f()(_,_) vec![137]) yields f(137, 0)
// call!(f(x1, x2)(_,_) vec![1,2,3,4,5]) yields f(x1, x2, 1, 2)
// Note that trailing *,* are ignored in function calls: f(a,b,) = f(a, b)
macro_rules! call {
    ($f: ident ($($x: expr), *) ($($y: pat), *) $z: expr) => {
        {
            let mut iter = $z.iter();
            $f($($x,)* $(
                            // nth(0) yields head of a list and moves iterator to the next value
                            match iter.nth(0usize) {
                                Some(t) => *t, // dereference
                                $y      => 0f64,}, // parameter = 0 if it is not set in z
                       )*
            )
        }
    };
}

// Конструктор электромагнитных полей. Возвращает функцию пространственных координат и времени,
// возвращающую напряжённости полей. ftype - тип функции, описывающей электромагнитное поле fparam
// - массив её параметров, помимо r и t.
fn fields(ftype: i32, fparam: Vec<f64>) -> Box<Fn(V, f64) -> EM> {
    // howto https://doc.rust-lang.org/beta/book/closures.html
    // see also https://github.com/rust-lang/rust/issues/24036
    let f: Box<Fn(V, f64) -> EM> = match ftype {
        0 => Box::new(move |r, t| call!(const_EB()(_,_) fparam)),
        1 => Box::new(move |r, t| call!(plane_wave(r,t)(_,_,_,_,_,_,_,_,_) fparam)),
        2 => Box::new(move |r, t| call!(const_B()(_) fparam)),
        3 => Box::new(move |r, t| call!(rot_E(t)(_,_,_,_,_,_,_) fparam)),
        4 => Box::new(move |r, t| call!(const_E_lin_B(r)(_) fparam)),
        5 => Box::new(move |r, t| call!(cir_pol_stw(r,t)(_) fparam)),
        6 => Box::new(move |r, t| call!(lin_pol_stw(r, t)(_) fparam)),
        7 => Box::new(move |r, t| call!(te11(r, t)(_, _, _) fparam)),
        8 => Box::new(move |r, t| call!(lin_pol_absorber(r, t)(_,_) fparam)),
        9 => Box::new(move |r, t| call!(dms_field(r, t)(_,_) fparam)),
        10 => {
            let nr = fparam[0] as usize;
            let nz = fparam[1] as usize;
            let r_m = fparam[2];
            let z_m = fparam[3];
            let a = test_A_phi((nr, nz), (r_m, z_m));
            Box::new(move |r, t| {
                fields_from_A_phi(r, t, (nr, nz), (r_m, z_m), &a)})
        },
        11 => {
            let nr = fparam[0] as usize;
            let nz = fparam[1] as usize;
            let r_m = fparam[2];
            let z_m = fparam[3];
            // J0 = fparam[4], barlambda = fparam[5], n_int = fparam[6]
            let a = current_loop_A_phi((nr, nz), (r_m, z_m), fparam[4], fparam[5], fparam[6] as i32);
            // если fparam[7] задан (любым числом), то производится обращение t -> -t, B -> -B.
            if fparam.len() < 7 {
                return Box::new(move |r, t| {
                    fields_from_A_phi(r, t, (nr, nz), (r_m, z_m), &a)});
            } else {
                return Box::new(move |r, t| {
                    let mut f = fields_from_A_phi(r, -t, (nr, nz), (r_m, z_m), &a);
                    f.bx = -f.bx;
                    f.by = -f.by;
                    f.bz = -f.bz;
                    f})
            }
        },
        12 => Box::new(move |r, t| call!(double_standing_wave(r, t)(_,_,_,_,_) fparam)), 
        13 => Box::new(move |r, t| call!(cosine_pulse(r, t)(_,_) fparam)),
        14 => Box::new(move |r, t| call!(gauss(r,t)(_,_) fparam)),
        _ => Box::new(     |r, t| EM { ex: 0.0, ey: 0.0, ez: 0.0, bx: 0.0, by: 0.0, bz: 0.0 }),
    };
    f
}

// Обёртка для drive_particle_new, возвращает указатель на массив координат и импульсов длиной (n +
// 1) * 6. s0 указывает на начальные условия.
#[no_mangle]
pub extern "C"
fn drive_particle_new_binding(s0: *const c_double, dt: c_double, n: size_t, ftype: c_int, nfparam:
                              size_t, fparam: *const c_double) -> *const c_double {
    let slice_s0 = unsafe { std::slice::from_raw_parts(s0 as *const f64, 6usize) };
    let s = Particle { x: slice_s0[0], y: slice_s0[1], z: slice_s0[2],
                       ux: slice_s0[3], uy: slice_s0[4], uz: slice_s0[5] };
    let slice_fparam = unsafe { std::slice::from_raw_parts(fparam as *const f64, nfparam as usize) };
    let mut fp: Vec<f64> = Vec::new();
    for x in slice_fparam {
        fp.push(*x)
    }
    let vec = drive_particle_new(s, fields(ftype as i32, fp), dt as f64, n as usize);
    let mut plain_vec: Vec<f64> = Vec::new();
    for p in vec {
        plain_vec.push(p.x);
        plain_vec.push(p.y);
        plain_vec.push(p.z);
        plain_vec.push(p.ux);
        plain_vec.push(p.uy);
        plain_vec.push(p.uz);
    }
    let p = plain_vec.as_ptr();
    std::mem::forget(plain_vec); // деструктор v не будет вызван; небольшая утечка памяти на length,
                                 // capacity etc.
    p as *const c_double
}

#[no_mangle]
pub extern "C" fn drive_particle_binding(s0: *const c_double, t0: c_double, dt: c_double, n: size_t, ftype: c_int, nfparam:
                              size_t, fparam: *const c_double, lambda: c_double, rf_scheme: c_int) -> *const c_double {
    let slice_s0 = unsafe { std::slice::from_raw_parts(s0 as *const f64, 6usize) };
    let s = Particle { x: slice_s0[0], y: slice_s0[1], z: slice_s0[2],
                       ux: slice_s0[3], uy: slice_s0[4], uz: slice_s0[5] };
    let slice_fparam = unsafe { std::slice::from_raw_parts(fparam as *const f64, nfparam as usize) };
    let mut fp: Vec<f64> = Vec::new();
    for x in slice_fparam {
        fp.push(*x)
    }
    let vec = drive_particle(s, fields(ftype as i32, fp), t0 as f64, dt as f64, n as usize, lambda as f64, rf_scheme as i32);
    let mut plain_vec: Vec<f64> = Vec::new();
    for p in vec {
        plain_vec.push(p.x);
        plain_vec.push(p.y);
        plain_vec.push(p.z);
        plain_vec.push(p.ux);
        plain_vec.push(p.uy);
        plain_vec.push(p.uz);
    }
    let p = plain_vec.as_ptr();
    std::mem::forget(plain_vec); // деструктор v не будет вызван; небольшая утечка памяти на length,
                                 // capacity etc.
    p as *const c_double
}

// variable-length vector
#[repr(C)]
pub struct Vlv {
    n: usize,
    ptr: *const c_double
}

// binding for EM struct
#[repr(C)]
pub struct c_EM {
    ex: c_double,
    ey: c_double,
    ez: c_double,
    bx: c_double,
    by: c_double,
    bz: c_double,
}

#[no_mangle]
pub extern "C"
// эта функция очень неэффективна в том смысле, что при каждом её вызове closure f строится заново;
// для массивов точек r и t поэтому "правильнее" использовать другую функцию, c_map_fields.
fn fields_binding(r: *const c_double, t: c_double, ftype: c_int, nfparam: size_t, fparam: *const
                  c_double) -> *const c_EM {
    let slice_r = unsafe { std::slice::from_raw_parts(r as *const f64, 3usize) };
    let r1 = (slice_r[0], slice_r[1], slice_r[2]);
    let slice_fparam = unsafe { std::slice::from_raw_parts(fparam as *const f64, nfparam as usize) };
    let mut fp: Vec<f64> = Vec::new();
    for x in slice_fparam {
        fp.push(*x);
    }
    let f = fields(ftype as i32, fp);
    let em = f(r1, t as f64);
    let f_boxed = Box::new(c_EM { ex: em.ex as c_double,
                                  ey: em.ey as c_double,
                                  ez: em.ez as c_double,
                                  bx: em.bx as c_double,
                                  by: em.by as c_double,
                                  bz: em.bz as c_double, });
    let p = f_boxed.as_ref() as *const c_EM;
    std::mem::forget(f_boxed);
    p
}

#[no_mangle]
pub extern "C"
// см. примечание к fields_binding; rt - массив вида [x0, y0, z0, t0, x1, y1, z1, t1, ...];
// возвращает массив [ex0, ey0, ez0, bx0, by0, bz0, ex1, ey1, ...]
fn c_map_fields(npoints: size_t, rt: *const c_double, ftype: c_int, nfparam: size_t,
                fparam: *const c_double) -> *const c_double {
    let slice_rt = unsafe { std::slice::from_raw_parts(rt as *const f64, npoints as usize * 4usize) };
    let mut rts: Vec<f64> = Vec::new();
    for x in slice_rt {
        rts.push(*x as f64);
    }
    let slice_fparam = unsafe { std::slice::from_raw_parts(fparam as *const f64, nfparam as usize) };
    let mut fp: Vec<f64> = Vec::new();
    for x in slice_fparam {
        fp.push(*x);
    }
    let f = fields(ftype as i32, fp);
    let mut result: Vec<f64> = Vec::new();
    for i in 0usize..(npoints as usize) {
        let j = 4usize * i;
        let r1 = (rts[j], rts[j + 1], rts[j + 2]);
        let t1 = rts[j + 3];
        let em = f(r1, t1);
        result.push(em.ex);
        result.push(em.ey);
        result.push(em.ez);
        result.push(em.bx);
        result.push(em.by);
        result.push(em.bz);
    }
    let p = result.as_ptr();
    std::mem::forget(result); // деструктор v не будет вызван; небольшая утечка памяти на length,
                                 // capacity etc.
    p as *const c_double
}

