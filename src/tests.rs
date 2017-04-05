use vec2d::*;
use consts::*;
use num;

trait ApproxEq {
    fn approx_eq(self, other: Self) -> bool;
}

// Only works well for numbers around 0.1-10
// Tiny numbers always pass, huge never pass
const APPROX_ERROR: f32 = 0.0001;
impl ApproxEq for f32 {
    fn approx_eq(self, other: Self) -> bool {
        (self - other).abs() < APPROX_ERROR
    }
}

impl ApproxEq for Vector {
    fn approx_eq(self, other: Self) -> bool {
        self.x.approx_eq(other.x) && self.y.approx_eq(other.y)
    }
}

impl ApproxEq for Point {
    fn approx_eq(self, other: Self) -> bool {
        self.x.approx_eq(other.x) && self.y.approx_eq(other.y)
    }
}

impl ApproxEq for Rotation {
    fn approx_eq(self, other: Self) -> bool {
        self.cos.approx_eq(other.cos) && self.sin.approx_eq(other.sin)
    }
}

impl ApproxEq for Transform {
    fn approx_eq(self, other: Self) -> bool {
        self.mat[0][0].approx_eq(other.mat[0][0]) && 
        self.mat[0][1].approx_eq(other.mat[0][1]) &&
        self.mat[1][0].approx_eq(other.mat[1][0]) && 
        self.mat[1][1].approx_eq(other.mat[1][1]) &&
        self.vec.approx_eq(other.vec)
    }
}

fn assert_approx<T>(a: T, b: T)
    where T: ApproxEq + ::std::fmt::Debug + Copy
{
    assert!(a.approx_eq(b),
            "Assert failed:\r\n{:?} does not approximate {:?}.",
            a, b);
}

#[test]
fn vec_basic_arithmetic() {
    assert_approx(
        Vector::new(2.0, 3.0) + Vector::new(4.0, 5.0),
        Vector::new(6.0, 8.0));

    assert_approx(
        Vector::new(2.0, 3.0) - Vector::new(4.0, 5.0),
        Vector::new(-2.0, -2.0));

    assert_approx(
        Vector::new(2.0, 3.0) * 4.0,
        Vector::new(8.0, 12.0));

    assert_approx(
        Vector::new(8.0, 12.0) / 4.0,
        Vector::new(2.0, 3.0));
}

// Tests the types of the partly macro-generated impls of the
// operator traits. If this compiles, the test is successful.
#[test]
fn impl_op_types() {
    let mut v = VEC_ZERO;
    let mut p = POINT_ORIGIN;
    let s = 1.0;

    v = v + v;
    p = p + v;
    p = v + p;

    v += v;
    p += v;

    v = v - v;
    p = p - v;
    v = p - p;

    v -= v;
    p -= v;

    v = v * s;
    v *= s;
    
    v = v / s;
    v /= s;

    v = -v;

    assert_approx(p.as_vec(), v);
}

#[test]
fn vec_angle_deg() {
    assert_approx(Vector::new(2.0, 2.0).angle_deg(), 45.0);
}

#[test]
fn vec_angle_rad() {
    assert_approx(VEC_UP.angle_rad(), TAU / 4.0);
}

#[test]
fn vec_angle_between() {
    let v1 = Vector::new(-1.0, 2.0);
    let v2 = Vector::new(3.0, -2.0);
    let angle = v1.angle_between_deg(v2);
    assert_approx(v1.rotate_deg(angle).angle_deg(), v2.angle_deg());
}

#[test]
fn vec_normalize() {
    let v = Vector::new(-4.0, 3.0);
    assert_approx(v.normalize().magnitude(), 1.0);
}

#[test]
fn vec_normalize_mut() {
    let mut v1 = Vector::new(-5.0, 2.0);
    let copy = v1.clone();
    v1.normalize_mut();
    assert_approx(v1, copy.normalize());
}

#[test]
fn vec_unit_from_deg() {
    assert_approx(Vector::unit_from_degrees(45.0), VEC_UP_RIGHT);
}

#[test]
fn vec_to_radians_and_back() {
    let vec = Vector::new(2.0, 3.0);
    assert_approx(vec, Vector::new_polar_rad(vec.magnitude(), vec.angle_rad()));
}

#[test]
fn vec_to_degrees_and_back() {
    let vec = Vector::new(2.0, 3.0);
    assert_approx(vec, Vector::new_polar_deg(vec.magnitude(), vec.angle_deg()));
}

#[test]
fn vec_scale() {
    let v1 = Vector::new(2.0, 3.0);
    let v2 = Vector::new(-3.0, 2.0);
    assert_approx(v1.scale(v2),
    Vector::new(-6.0, 6.0));
}

#[test]
fn vec_scale_mut() {
    let mut v1 = Vector::new(2.0, 3.0);
    let v2 = v1.clone();
    let factor = Vector::new(-1.0, 2.0);
    v1.scale_mut(factor);
    assert_approx(v1, v2.scale(factor));
}

#[test]
fn vec_magnitude() {
    assert_approx(Vector::new(3.0, 4.0).magnitude(), 5.0);
}

#[test]
fn vec_rotate_direction() {
    assert_approx(VEC_UP, VEC_RIGHT.rotate_deg(90.0));
}

#[test]
fn vec_dot() {
    assert_approx(Vector::new(-6.0, 8.0).dot(Vector::new(5.0, 12.0)), 66.0);
}

#[test]
fn vec_dot_perpendicular() {
    assert_approx(Vector::new(2.0, 3.0).dot(Vector::new(3.0, -2.0)), 0.0);
}

#[test]
fn vec_dist() {
    let v1 = Vector::new(-2.0, -2.0);
    let v2 = Vector::new(1.0, 2.0);
    assert_approx(v1.dist(v2), 5.0);
}

#[test]
fn vec_grid_dist() {
    let v1 = Vector::new(-1.0, -2.0);
    let v2 = Vector::new(4.0, 3.0);
    assert_approx(v1.grid_dist(v2), 10.0);
}

#[test]
fn vec_map() {
    assert_approx(
        Vector::new(2.0, 3.0).map(&|x| -x),
        Vector::new(-2.0, -3.0));
}

#[test]
fn vec_map_mut() {
    let mut v1 = Vector::new(6.0, 3.0);
    let v2 = v1.clone();
    let func = |x: f32| -x;
    v1.map_mut(&func);
    assert_approx(v1, v2.map(&func));
}

#[test]
fn vec_lerp_average() {
    let v1 = Vector::new(2.0, 5.0);
    let v2 = Vector::new(-3.0, 2.0);
    assert_approx(v1.lerp(v2, 0.5), (v1 + v2) / 2.0);
}

#[test]
fn vec_lerp_one() {
    let v1 = Vector::new(2.0, 5.0);
    let v2 = Vector::new(-3.0, 2.0);
    assert_approx(v1.lerp(v2, 1.0), v2);
}

#[test]
fn point_rotate_pivot() {
    let pivot = Point::new(5.0, 2.0);
    assert_approx((VEC_UP_LEFT + pivot).rotate_deg(pivot, 45.0),
                  (VEC_LEFT + pivot));
}

#[test]
fn rot_add_sub() {
    let r1 = Rotation::new_deg(60.0);
    let r2 = Rotation::new_deg(45.0);
    assert_approx(r1.add(r2), Rotation::new_deg(105.0));
    assert_approx(r1.sub(r2), Rotation::new_deg(15.0));
}

#[test]
fn rot_to_from_angle() {
    let angle = 34.0;
    assert_approx(Rotation::new_deg(angle).angle_deg(), angle);
}

#[test]
fn rot_consts() {
    assert_approx(ROTATION_0_DEG, Rotation::new_deg(0.0));
    assert_approx(ROTATION_45_DEG, Rotation::new_deg(45.0));
    assert_approx(ROTATION_90_DEG, Rotation::new_deg(90.0));
    assert_approx(ROTATION_135_DEG, Rotation::new_deg(135.0));
    assert_approx(ROTATION_180_DEG, Rotation::new_deg(180.0));
    assert_approx(ROTATION_225_DEG, Rotation::new_deg(225.0));
    assert_approx(ROTATION_270_DEG, Rotation::new_deg(270.0));
    assert_approx(ROTATION_315_DEG, Rotation::new_deg(315.0));
}

#[test]
fn rot_90_deg() {
    assert_approx(VEC_UP.rotate(ROTATION_90_DEG), VEC_LEFT);
}

#[test]
fn rot_between() {
    let v1 = Vector::new(4.0, 5.0);
    let v2 = Vector::new(-2.0, 3.0);
    assert_approx(v2.angle_deg(), v1.rotate(Rotation::between(v1, v2)).angle_deg());
}

#[test]
fn rot_negate() {
    let angle = 50.0;
    let rot = Rotation::new_deg(angle);
    assert_approx(rot.negate().angle_deg(), -angle);
}

#[test]
fn transform_point_identity() {
    let p = Point::new(-3.0, 2.0);
    assert_approx(TRANSFORM_IDENTITY.transform(p), p);
}

#[test]
fn transform_new_rad_identity() {
    let tf = Transform::new_rad(VEC_ZERO, 0.0, VEC_ONE);
    assert_approx(tf, TRANSFORM_IDENTITY);
}

#[test]
fn transform_point() {
    let tf = Transform::new_deg(
        Vector::new(2.0, 2.0),
        90.0,
        Vector::new(3.0, 2.0));
    let p = Point::new(4.0, 3.0);
    // Result (-4, 14) derived by hand
    assert_approx(tf.transform(p), Point::new(-4.0, 14.0));
}

#[test]
fn transform_inverse_linear_vector_point() {
    let v = Vector::new(-3.0, -2.0);
    let p = v.as_point();
    let tf = Transform::new_matrix(2.0, 3.0, 4.0, 5.0);
    let new_v = tf.inverse_transform_checked(v).unwrap();
    let new_p = tf.inverse_transform_checked(p).unwrap();
    assert_approx(new_v, new_p.as_vec());
}

#[test]
fn transform_global_to_local_vec() {
    let v = Vector::new(5.0, 2.0);
    let tf = Transform::new_deg(
        Vector::new(3.0, 2.0),
        32.0,
        Vector::new(3.0, 2.0));
    let transformed = tf.transform(v);
    let back = tf.inverse_transform_checked(transformed).unwrap();
    assert_approx(v, back);
}

#[test]
fn transform_global_to_local_point() {
    let p = Point::new(-2.0, 3.0);
    let tf = Transform::new_deg(
        Vector::new(6.0, 5.0),
        43.0,
        Vector::new(2.0, 1.0));
    let transformed = tf.transform(p);
    let back = tf.inverse_transform_checked(transformed).unwrap();
    assert_approx(p, back);
}

#[test]
fn transform_translate() {
    let tf1 = Transform::new_deg(
        Vector::new(2.0, 3.0),
        45.0,
        Vector::new(6.0, 7.0));
    let tf2 = Transform::new_deg(
        Vector::new(4.0, 5.0),
        45.0,
        Vector::new(6.0, 7.0));
    assert_approx(tf1.translate(Vector::new(2.0, 2.0)), tf2);
}

#[test]
fn transform_translate_mut() {
    let mut tf1 = Transform::new_deg(
        Vector::new(-3.0, -2.0),
        -54.0,
        Vector::new(-7.0, -6.0));
    let tf2 = tf1.clone();
    let v = Vector::new(4.0, -3.0);
    tf1.translate_mut(v);
    assert_approx(tf1, tf2.translate(v));
}

#[test]
fn transform_translate_local_mut() {
    let mut tf1 = Transform::new_deg(
        Vector::new(1.0, 2.0),
        34.0,
        Vector::new(5.0, 6.0));
    let tf2 = tf1.clone();
    let v = Vector::new(3.0, 4.0);
    tf1.translate_local_mut(v);
    assert_approx(tf1, tf2.translate_local(v));
}

#[test]
fn transform_combine_identity() {
    let tf = Transform::new_deg(
        Vector::new(2.0, 3.0),
        45.0,
        Vector::new(-5.0, -6.0));
    let id = TRANSFORM_IDENTITY;
    assert_approx(id.transform(tf.transform(id)), tf);
}

#[test]
fn transform_combine() {
    // Tests the associative property, that for two matrices M, N
    // and a column vector v it holds that (M*N)*v = M*(N*v)
    let p = Point::new(3.0, 4.0);
    let coordsys1 = Transform::new_deg(
        Vector::new(2.0, 3.0),
        76.0,
        Vector::new(3.0, 2.0));
    let coordsys2 = Transform::new_deg(
        Vector::new(-1.0, 2.0),
        -64.0,
        Vector::new(6.0, 2.0));
    let comb = coordsys1.transform(coordsys2);
    assert_approx(coordsys1.transform(coordsys2.transform(p)),
                  comb.transform(p));
}

#[test]
fn transform_combine_mut() {
    let coordsys1 = Transform::new_deg(
        Vector::new(2.0, 3.0),
        76.0,
        Vector::new(3.0, 2.0));
    let mut coordsys2 = Transform::new_deg(
        Vector::new(-1.0, 2.0),
        -64.0,
        Vector::new(6.0, 2.0));
    let copy = coordsys2.clone();
    coordsys1.transform_mut(&mut coordsys2);
    assert_approx(coordsys2, coordsys1.transform(copy));
}

#[test]
fn num_lerp_average() {
    let a = 1.231234;
    let b = 5.43123;
    assert_approx(num::lerp(a, b, 0.5), (a + b) / 2.0);
}

#[test]
fn num_clamp_commut() {
    let a = 2.0;
    let b = -3.0;
    assert_approx(num::clamp(100.0, a, b), num::clamp(100.0, b, a));
}

#[test]
fn num_range_map() {
    let x = 7.0;
    let y = num::range_map(x, 0.0, 5.0, 0.0, 1.0);
    assert_approx(y, 1.4);
}

#[test]
fn num_modulo_negative() {
    assert_approx(num::modulo(-167.0, 100.0), 33.0);
}
