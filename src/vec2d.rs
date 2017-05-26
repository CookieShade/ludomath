//! Functions and other utilities for vector math.
//!
//! The conventions of this library are:
//! The positive X-axis points rightward, and the positive Y-axis
//! points upward. Angles are 0 at (1, 0), increasing counter-clockwise.
//!
//! Many methods have a _mut variant, these are the assignment operators.
//! That is, if the standard variant is +, the _mut variant is +=.

use std::fmt;
use std::ops::{Add, AddAssign, Sub, SubAssign, Neg, Mul, MulAssign, Div, DivAssign};
use std::mem;
use consts::*;
use num;

/// Represents a cartesian vector in the 2D euclidean plane.
#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct Vector {
    /// The X component of the vector.
    pub x: f32,
    /// The Y component of the vector.
    pub y: f32,
}

/// Represents a cartesian point in the 2D euclidean plane.
#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct Point {
    /// The X coordinate of the point.
    pub x: f32,
    /// The Y coordinate of the point.
    pub y: f32,
}

impl Vector {
    /// Constructs a new vector.
    #[inline]
    pub fn new(x: f32, y: f32) -> Vector {
        Vector { x: x, y: y }
    }

    /// Constructs a unit-length vector with a certain angle in radians.
    #[inline]
    pub fn unit_from_radians(angle: f32) -> Vector {
        let (sin, cos) = angle.sin_cos();
        Vector { x: cos, y: sin }
    }

    /// Constructs a unit-length vector with a certain angle in degrees.
    #[inline]
    pub fn unit_from_degrees(angle: f32) -> Vector {
        let (sin, cos) = angle.to_radians().sin_cos();
        Vector { x: cos, y: sin }
    }

    /// Constructs a vector of a certain magnitude and angle (in radians).
    #[inline]
    pub fn new_polar_rad(radius: f32, angle: f32) -> Vector {
        let (sin, cos) = angle.sin_cos();
        Vector {
            x: cos * radius,
            y: sin * radius,
        }
    }

    /// Constructs a vector of a certain magnitude and angle (in degrees).
    #[inline]
    pub fn new_polar_deg(radius: f32, angle: f32) -> Vector {
        Vector::new_polar_rad(radius, angle.to_radians())
    }

    /// Computes the magnitude/length of a vector.
    #[inline]
    pub fn magnitude(self) -> f32 {
        self.x.hypot(self.y)
    }

    /// Computes the dot product of two vectors.
    #[inline]
    pub fn dot(self, other: Vector) -> f32 {
        self.x.mul_add(other.x, self.y * other.y)
    }

    /// Computes the component-wise product of two vectors.
    #[inline]
    pub fn scale(self, other: Vector) -> Vector {
        Vector {
            x: self.x * other.x,
            y: self.y * other.y,
        }
    }

    /// Computes the component-wise product of two vectors,
    /// assigning the result to the first one.
    #[inline]
    pub fn scale_mut(&mut self, other: Vector) {
        self.x *= other.x;
        self.y *= other.y;
    }

    /// Rotates a vector by a `Rotation` struct.
    #[inline]
    pub fn rotate(self, rotation: Rotation) -> Vector {
        Vector {
            x: self.x * rotation.cos - self.y * rotation.sin,
            y: self.x * rotation.sin + self.y * rotation.cos,
        }
    }
     
    /// Rotates a vector, mutating it.
    #[inline]
    pub fn rotate_mut(&mut self, rotation: Rotation) {
        *self = self.rotate(rotation);
    }

    /// Rotates a vector by an angle in radians.
    #[inline]
    pub fn rotate_rad(self, angle: f32) -> Vector {
        self.rotate(Rotation::new_rad(angle))
    }

    /// Rotates a vector by an angle in radians, mutating it.
    #[inline]
    pub fn rotate_rad_mut(&mut self, angle: f32) {
        *self = self.rotate_rad(angle);
    }

    /// Rotates a vector by an angle in degrees.
    #[inline]
    pub fn rotate_deg(self, angle: f32) -> Vector {
        self.rotate(Rotation::new_deg(angle))
    }

    /// Rotates a vector by an angle in degrees, mutating it.
    #[inline]
    pub fn rotate_deg_mut(&mut self, angle: f32) {
        *self = self.rotate_deg(angle);
    }

    /// Normalizes a vector to unit length.
    ///
    /// If the vector is zero, returns (1, 0).
    #[inline]
    pub fn normalize(self) -> Vector {
        if self == VEC_ZERO {
            VEC_RIGHT
        } else {
            self * self.magnitude().recip()
        }
    }

    /// Sets a vector to unit length.
    ///
    /// If the vector is zero, sets to (1, 0).
    #[inline]
    pub fn normalize_mut(&mut self) {
        if *self == VEC_ZERO {
            *self = VEC_RIGHT
        } else {
            *self *= self.magnitude().recip()
        }
    }

    /// Computes the euclidean distance between two vectors.
    #[inline]
    pub fn dist(self, other: Vector) -> f32 {
        (self - other).magnitude()
    }

    /// Computes the grid-based distance between two vectors.
    ///
    /// This is the length of the shortest path between the vectors,
    /// if the only movement allowed is parallel to the x- or y-axis.
    #[inline]
    pub fn grid_dist(self, other: Vector) -> f32 {
        (self.x - other.x).abs() + (self.y - other.y).abs()
    }

    /// Computes the angle between a vector and the positive x-axis, in radians.
    ///
    /// Return values are contained in (-τ/2, τ/2], aka. (-π, π].
    #[inline]
    pub fn angle_rad(self) -> f32 {
        self.y.atan2(self.x)
    }

    /// Computes the angle between a vector and the positive x-axis, in degrees.
    ///
    /// Return values are approximately contained in (-180, 180].
    #[inline]
    pub fn angle_deg(self) -> f32 {
        self.angle_rad().to_degrees()
    }

    /// Computes the angle between two vectors, in radians.
    ///
    /// Rotating the first vector by the resulting angle will give it
    /// the same angle as the second.
    ///
    /// Return values are contained in (-τ/2, τ/2], aka (-π, π].
    #[inline]
    pub fn angle_between_rad(self, other: Vector) -> f32 {
        let x = self.x * other.x + self.y * other.y;
        let y = self.x * other.y - self.y * other.x;
        y.atan2(x)
    }

    /// Computes the angle between two vectors, in degrees.
    ///
    /// Return values are approximately contained in (-180, 180].
    #[inline]
    pub fn angle_between_deg(self, other: Vector) -> f32 {
        self.angle_between_rad(other).to_degrees()
    }

    /// Calls a function on each component of a vector,
    /// returning the result as a new vector.
    #[inline]
    pub fn map<F: Fn(f32) -> f32>(self, func: &F) -> Vector {
        Vector {
            x: func(self.x),
            y: func(self.y),
        }
    }

    /// Sets each component to the result of calling a function
    /// on the component's value.
    #[inline]
    pub fn map_mut<F: Fn(f32) -> f32>(&mut self, func: &F) {
        self.x = func(self.x);
        self.y = func(self.y);
    }

    /// Linearly interpolates from one vector to another by `t`.
    ///
    /// `t` < 0 or `t` > 1 gives results outside the
    /// rectangular bounds of the vectors.
    #[inline]
    pub fn lerp(self, other: Vector, t: f32) -> Vector {
        self + (other - self) * t
    }

    /// Returns `true` if `self` is within the (inclusive) bounds
    /// of the rectangle defined by two vertices.
    #[inline]
    pub fn is_in_rect(self, vert1: Vector, vert2: Vector) -> bool {
        let min_x = vert1.x.min(vert2.x);
        let max_x = vert1.x.max(vert2.x);
        let min_y = vert1.y.min(vert2.y);
        let max_y = vert1.y.max(vert2.y);

        self.x >= min_x && self.x <= max_x &&
        self.y >= min_y && self.y <= max_y
    }

    /// Clamps a vector to within the (inclusive) bounds of the rectangle
    /// defined by two vertices.
    #[inline]
    pub fn clamp_to_rect(self, vert1: Vector, vert2: Vector) -> Vector {
        Vector {
            x: num::clamp(self.x, vert1.x, vert2.x),
            y: num::clamp(self.y, vert1.y, vert2.y),
        }
    }

    /// Reinterprets the x- and y-components of a vector as a point.
    #[inline]
    pub fn as_point(self) -> Point {
        Point::new(self.x, self.y)
    }

    /// Converts a OpenGL-style `vec2` to a `Vector`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn from_vec2(vec2: [f32; 2]) -> Vector {
        Vector::new(vec2[0], vec2[1])
    }

    /// Converts a `Vector` to a OpenGL-style `vec2`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn to_vec2(self) -> [f32; 2] {
        [self.x, self.y]
    }

    /// Converts a `Vector` to a OpenGL-style `vec3`.
    ///
    /// The result is padded to `[x, y, 0]`.
    #[inline]
    pub fn to_vec3(self) -> [f32; 3] {
        [self.x, self.y, 0.0]
    }

    /// Converts a `Vector` to a OpenGL-style `vec4`.
    ///
    /// The result is padded to `[x, y, 0, 0]`.
    #[inline]
    pub fn to_vec4(self) -> [f32; 4] {
        [self.x, self.y, 0.0, 0.0]
    }

    /// Reinterprets a slice of `vec2`s as `Vector`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_from_vec2s(arr: &[[f32; 2]]) -> &[Vector] {
        unsafe { mem::transmute(arr) }
    }

    /// Reinterprets a slice of `Vector`s as `vec2`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_to_vec2s(arr: &[Vector]) -> &[[f32; 2]] {
        unsafe { mem::transmute(arr) }
    }
}

impl Point {
    /// Constructs a new point from x- and y-components.
    #[inline]
    pub fn new(x: f32, y: f32) -> Point {
        Point { x: x, y: y }
    }

    /// Constructs a new point from polar coordinates.
    ///
    /// The angle is given in radians.
    #[inline]
    pub fn new_polar_rad(radius: f32, angle: f32) -> Point {
        let (sin, cos) = angle.sin_cos();
        Point::new(radius * cos, radius * sin)
    }

    /// Constructs a new point from polar coordinates.
    ///
    /// The angle is given in degrees.
    #[inline]
    pub fn new_polar_deg(radius: f32, angle: f32) -> Point {
        Point::new_polar_rad(radius, angle.to_radians())
    }

    /// Computes the euclidean distance between two points.
    #[inline]
    pub fn dist(self, other: Point) -> f32 {
        (self - other).magnitude()
    }

    /// Rotates a point around a pivot point.
    #[inline]
    pub fn rotate(self, pivot: Point, rotation: Rotation) -> Point {
        (self - pivot).rotate(rotation) + pivot
    }

    /// Rotates a point `angle` radians counter-clockwise around a pivot point.
    #[inline]
    pub fn rotate_rad(self, pivot: Point, angle: f32) -> Point {
        self.rotate(pivot, Rotation::new_rad(angle))
    }

    /// Rotates a point `angle` degrees counter-clockwise around a pivot point.
    #[inline]
    pub fn rotate_deg(self, pivot: Point, angle: f32) -> Point {
        self.rotate_rad(pivot, angle.to_radians())
    }

    /// Linearly interpolates from one point to another by `t`.
    ///
    /// `t` < 0 or `t` > 1 gives results outside the
    /// rectangular bounds of the points.
    #[inline]
    pub fn lerp(self, other: Point, t: f32) -> Point {
        self + (other - self) * t
    }

    /// Returns `true` if `self` is within the (inclusive) bounds
    /// of the rectangle defined by two points.
    #[inline]
    pub fn is_in_rect(self, vert1: Point, vert2: Point) -> bool {
        let min_x = vert1.x.min(vert2.x);
        let max_x = vert1.x.max(vert2.x);
        let min_y = vert1.y.min(vert2.y);
        let max_y = vert1.y.max(vert2.y);

        self.x >= min_x && self.x <= max_x &&
        self.y >= min_y && self.y <= max_y
    }

    /// Clamps a vector to within the (inclusive) bounds of the rectangle
    /// defined by two points.
    #[inline]
    pub fn clamp_to_rect(self, vert1: Point, vert2: Point) -> Point {
        Point {
            x: num::clamp(self.x, vert1.x, vert2.x),
            y: num::clamp(self.y, vert1.y, vert2.y),
        }
    }

    /// Reinterprets the x- and y-components of a point as a vector.
    #[inline]
    pub fn as_vec(self) -> Vector {
        Vector::new(self.x, self.y)
    }

    /// Changes the type of a OpenGL-style `vec2` array to a `Point`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    /// Useful mostly for interop with other libraries, like OpenGL wrappers.
    #[inline]
    pub fn from_gl_vec2(vec2: [f32; 2]) -> Point {
        Point::new(vec2[0], vec2[1])
    }

    /// Changes the type of a `Point` to a OpenGL-style `vec2` array.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    /// Useful mostly for interop with other libraries, like OpenGL wrappers.
    #[inline]
    pub fn to_gl_vec2(self) -> [f32; 2] {
        [self.x, self.y]
    }

    /// Converts a OpenGL-style `vec2` to a `Point`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn from_vec2(vec2: [f32; 2]) -> Point {
        Point::new(vec2[0], vec2[1])
    }

    /// Converts a `Point` to a OpenGL-style `vec2`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn to_vec2(self) -> [f32; 2] {
        [self.x, self.y]
    }

    /// Converts a `Point` to a OpenGL-style `vec3`.
    ///
    /// The result is padded to `[x, y, 0]`.
    #[inline]
    pub fn to_vec3(self) -> [f32; 3] {
        [self.x, self.y, 0.0]
    }

    /// Converts a `Point` to a OpenGL-style `vec4`.
    ///
    /// The result is padded to `[x, y, 0, 1]`.
    #[inline]
    pub fn to_vec4(self) -> [f32; 4] {
        [self.x, self.y, 0.0, 1.0]
    }

    /// Reinterprets a slice of `vec2`s as `Point`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_from_vec2s(arr: &[[f32; 2]]) -> &[Point] {
        unsafe { mem::transmute(arr) }
    }

    /// Reinterprets a slice of `Point`s as `vec2`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_to_vec2s(arr: &[Point]) -> &[[f32; 2]] {
        unsafe { mem::transmute(arr) }
    }
}

/// Represents a 2D rotation or direction.
///
/// Stores the sine and cosine of the angle of rotation, which makes
/// the rotation faster if the same angle is used multiple times.
/// The drawbacks are memory usage, and that 
/// combining multiple rotations is more expensive.
#[repr(C)]
#[derive(Clone, Copy, PartialEq)]
pub struct Rotation {
    /// The cosine of the angle.
    pub cos: f32,
    /// The sine of the angle.
    pub sin: f32,
}

impl Rotation {
    /// Constructs a rotation from the direction of a vector.
    ///
    /// The angle of the rotation is the angle from the positive
    /// x-axis to the vector.
    #[inline]
    pub fn new_in_direction(dir: Vector) -> Rotation {
        let normalized = dir.normalize();
        Rotation {
            cos: normalized.x,
            sin: normalized.y,
        }
    }

    /// Constructs a rotation of a specified angle, in radians.
    #[inline]
    pub fn new_rad(angle: f32) -> Rotation {
        let (sin, cos) = angle.sin_cos();
        Rotation {
            cos: cos,
            sin: sin,
        }
    }

    /// Constructs a rotation of a specified angle, in degrees.
    #[inline]
    pub fn new_deg(angle: f32) -> Rotation {
        Rotation::new_rad(angle.to_radians())
    }

    /// Gets the rotation from one vector to another.
    ///
    /// Rotating the first vector by the returned rotation yields
    /// a vector with the same direction as the second vector.
    #[inline]
    pub fn between(from: Vector, to: Vector) -> Rotation {
        Rotation::new_in_direction(Vector {
            x: from.x * to.x + from.y * to.y,
            y: from.x * to.y - from.y * to.x,
        })
    }

    /// Negates the angle of a rotation.
    #[inline]
    pub fn negate(self) -> Rotation {
        // Simply uses that: cos -x = cos x, sin -x = -sin x
        Rotation {
            cos: self.cos,
            sin: -self.sin,
        }
    }

    /// Adds the angles of two rotations together.
    #[inline]
    pub fn add(self, other: Rotation) -> Rotation {
        // Equivalent to complex multiplication where `cos`
        // is the real component and `sin` is the imaginary.
        Rotation {
            cos: self.cos * other.cos - self.sin * other.sin,
            sin: self.sin * other.cos + self.cos * other.sin,
        }
    }

    /// Subtracts the angle of one rotation from another.
    ///
    /// Equivalent to `self.add(other.negate())`.
    #[inline]
    pub fn sub(self, other: Rotation) -> Rotation {
        // Equivalent to `add()` but with the sign of `other.sin` reversed.
        // Uses that cos -x = cos x, sin -x = -sin x.
        Rotation {
            cos: self.cos * other.cos + self.sin * other.sin,
            sin: self.sin * other.cos - self.cos * other.sin,
        }
    }

    /// Returns the angle of rotation, in radians.
    #[inline]
    pub fn angle_rad(self) -> f32 {
        self.sin.atan2(self.cos)
    }

    /// Returns the angle of rotation, in degrees.
    #[inline]
    pub fn angle_deg(self) -> f32 {
        self.angle_rad().to_degrees()
    }
}

impl fmt::Debug for Rotation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
               "Rotation(angle: {} degrees, cos: {}, sin: {})",
               self.angle_deg(),
               self.cos,
               self.sin)
    }
}

/// Represents an affine transformation of the plane.
///
/// Can represent any linear transformation, plus translation.
/// Equivalent in layout to a GLSL (column-major) 3x2 matrix,
/// where the rightmost column contains the translation vector.
#[derive(Clone, Copy, PartialEq)]
pub struct Transform {
    /// The 2x2 linear transformation matrix.
    ///
    /// Stored in column-major order as a 2D array.
    pub mat: [[f32; 2]; 2],
    /// The translation vector.
    pub vec: Vector,
}

impl Transform {
    /// Constructs a new `Transform` from translation, rotation and scale.
    #[inline]
    pub fn new(translation: Vector, rotation: Rotation, scale: Vector) -> Transform {
        Transform {
            vec: translation,
            mat: [
                [scale.x * rotation.cos,
                 scale.x * rotation.sin],
                [scale.y * -rotation.sin,
                 scale.y * rotation.cos],
            ],
        }
    }

    /// Constructs a new `Transform` from position, rotation and scale.
    ///
    /// Rotation is given as an angle in radians.
    #[inline]
    pub fn new_rad(translation: Vector, angle: f32, scale: Vector) -> Transform {
        Transform::new(translation, Rotation::new_rad(angle), scale)
    }

    /// Constructs a new `Transform` from position, rotation and scale.
    ///
    /// Rotation is given as an angle in degrees.
    #[inline]
    pub fn new_deg(translation: Vector, angle: f32, scale: Vector) -> Transform {
        Transform::new(translation, Rotation::new_deg(angle), scale)
    }

    /// Constructs a `Transform` from a 2x2 linear tranformation matrix.
    ///
    /// The order of the arguments is their order in memory - column-major.
    #[inline]
    pub fn new_matrix(mat1_1: f32, mat2_1: f32, mat1_2: f32, mat2_2: f32) -> Transform {
        Transform {
            vec: VEC_ZERO,
            mat: [[mat1_1, mat2_1],
                  [mat1_2, mat2_2]],
        }
    }

    /// Constructs a scaling `Transform`.
    #[inline]
    pub fn new_scaling(scale: Vector) -> Transform {
        Transform::new_matrix(scale.x, 0.0, 0.0, scale.y)
    }

    /// Constructs a rotating `Transform` from a `Rotation`.
    #[inline]
    pub fn new_rotation(rotation: Rotation) -> Transform {
        Transform::new_matrix(rotation.cos, -rotation.sin, rotation.sin, rotation.cos)
    }

    /// Constructs a rotating `Transform` from an angle in radians.
    #[inline]
    pub fn new_rotation_rad(angle: f32) -> Transform {
        Transform::new_rotation(Rotation::new_rad(angle))
    }

    /// Constructs a rotating `Transform` from an angle in degrees.
    #[inline]
    pub fn new_rotation_deg(angle: f32) -> Transform {
        Transform::new_rotation(Rotation::new_deg(angle))
    }

    /// Constructs a translating `Transform` from an offset.
    #[inline]
    pub fn new_translation(offset: Vector) -> Transform {
        Transform {
            vec: offset,
            mat: TRANSFORM_IDENTITY.mat,
        }
    }

    /// Creates a new `Transform` by translating an existing one.
    ///
    /// The offset is relative to the global coordinate system.
    #[inline]
    pub fn translate(self, offset: Vector) -> Transform {
        Transform { vec: self.vec + offset, ..self }
    }

    /// Translates a `Transform`, mutating it.
    ///
    /// The offset is relative to the global coordinate system.
    #[inline]
    pub fn translate_mut(&mut self, offset: Vector) {
        self.vec += offset;
    }

    /// Creates a new `Transform` by translating an existing one.
    ///
    /// The offset is relative to the local coordinate system defined
    /// by the `Transform`.
    #[inline]
    pub fn translate_local(self, offset: Vector) -> Transform {
        Transform { vec: self.vec + self.transform(offset), ..self }
    }

    /// Translates a `Transform`, mutating it.
    ///
    /// The offset is relative to the local coordinate system defined
    /// by the `Transform`.
    #[inline]
    pub fn translate_local_mut(&mut self, offset: Vector) {
        self.vec += self.transform(offset);
    }

    /// Creates a new `Transform` by rotating an existing one.
    #[inline]
    pub fn rotate(self, rotation: Rotation) -> Transform {
        Transform {
            mat: [
                [self.mat[0][0] * rotation.cos + self.mat[1][0] * rotation.sin,
                 self.mat[0][1] * rotation.cos + self.mat[1][1] * rotation.sin],
                [self.mat[1][0] * rotation.cos - self.mat[0][0] * rotation.sin,
                 self.mat[1][1] * rotation.cos - self.mat[0][1] * rotation.sin],
            ],
            ..self
        }
    }

    /// Creates a new `Transform` by rotating an existing one
    /// by an angle in radians.
    #[inline]
    pub fn rotate_rad(self, angle: f32) -> Transform {
        self.rotate(Rotation::new_rad(angle))
    }

    /// Creates a new `Transform` by rotating an existing one
    /// by an angle in degrees.
    #[inline]
    pub fn rotate_deg(self, angle: f32) -> Transform {
        self.rotate(Rotation::new_deg(angle))
    }

    /// Rotates a `Transform`, mutating it.
    #[inline]
    pub fn rotate_mut(&mut self, rotation: Rotation) {
        *self = self.rotate(rotation);
    }

    /// Rotates a `Transform` by an angle in radians, mutating it.
    #[inline]
    pub fn rotate_rad_mut(&mut self, angle: f32) {
        self.rotate_mut(Rotation::new_rad(angle));
    }

    /// Rotates a `Transform` by an angle in degrees, mutating it.
    #[inline]
    pub fn rotate_deg_mut(&mut self, angle: f32) {
        self.rotate_mut(Rotation::new_deg(angle));
    }

    /// Creates a new `Transform` by scaling an existing one.
    #[inline]
    pub fn scale(self, factor: Vector) -> Transform {
        Transform {
            mat: [
                [self.mat[0][0] * factor.x, self.mat[0][1] * factor.x],
                [self.mat[1][0] * factor.y, self.mat[1][1] * factor.y],
            ],
            ..self
        }
    }

    /// Creates a new `Transform` by scaling an existing one.
    ///
    /// Scales by the same amount in both axes.
    #[inline]
    pub fn scale_by_num(self, factor: f32) -> Transform {
        Transform {
            mat: [
                [self.mat[0][0] * factor, self.mat[0][1] * factor],
                [self.mat[1][0] * factor, self.mat[1][1] * factor],
            ],
            ..self
        }
    }

    /// Scales a `Transform`, mutating it.
    #[inline]
    pub fn scale_mut(&mut self, factor: Vector) {
        self.mat[0][0] *= factor.x;
        self.mat[0][1] *= factor.x;
        self.mat[1][0] *= factor.y;
        self.mat[1][1] *= factor.y;
    }

    /// Scales a `Transform` equally in both axes, mutating it.
    #[inline]
    pub fn scale_by_num_mut(&mut self, factor: f32) {
        self.mat[0][0] *= factor;
        self.mat[0][1] *= factor;
        self.mat[1][0] *= factor;
        self.mat[1][1] *= factor;
    }

    /// Transforms a list of transformable items into a new `Vec`.
    #[inline]
    pub fn transform_arr<T>(self, arr: &[T]) -> Vec<T> where
        Self: Transformation<T>, T: Copy 
    {
        arr.iter().map(|el| self.transform(*el)).collect()
    }

    /// Transforms a list of transformable items, mutating them.
    #[inline]
    pub fn transform_arr_mut<T>(self, arr: &mut [T]) where
        Self: Transformation<T>, T: Copy 
    {
        for el in arr {
            self.transform_mut(el);
        }
    }

    /// Converts a OpenGL-style `mat3x2` to a `Transform`.
    ///
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn from_mat3x2(mat3x2: [[f32; 2]; 3]) -> Transform {
        Transform {
            mat: [mat3x2[0], mat3x2[1]],
            vec: Vector::new(mat3x2[2][0], mat3x2[2][1]),
        }
    }

    /// Converts a `Transform` to a OpenGL-style `mat3x2`.
    ///
    /// 
    /// Optimizes into a simple copy, since the types have the same layout.
    #[inline]
    pub fn to_mat3x2(self) -> [[f32; 2]; 3] {
        [
            self.mat[0],
            self.mat[1],
            [self.vec.x, self.vec.y],
        ]
    }

    /// Converts a `Transform` to a OpenGL-style `mat4x4`.
    ///
    /// The matrix is inserted in the top-left corner, and
    /// the translation vector in the top of the rightmost column.
    /// The rest of the elements are set to the identity matrix.
    #[inline]
    pub fn to_mat4x4(self) -> [[f32; 4]; 4] {
        [
            [self.mat[0][0], self.mat[0][1], 0.0, 0.0],
            [self.mat[1][0], self.mat[1][1], 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [self.vec.x, self.vec.y, 0.0, 1.0],
        ]
    }

    /// Reinterprets a slice of OpenGL-style `mat4x4`s as `Transform`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_from_mat3x2s(arr: &[[[f32; 4]; 4]]) -> &[Transform] {
        unsafe { mem::transmute(arr) }
    }

    /// Reinterprets a slice of `Transform`s as OpenGL-style `mat4x4`s.
    ///
    /// Does not allocate a new array, only changes the type of the pointer.
    #[inline]
    pub fn arr_to_mat3x2s(arr: &[Transform]) -> &[[[f32; 4]; 4]] {
        unsafe { mem::transmute(arr) }
    }
}

impl fmt::Debug for Transform {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
            "Transform:\n[{}, {}]\n[{}, {}],\n{:?}",
            self.mat[0][0], self.mat[1][0],
            self.mat[0][1], self.mat[1][1],
            self.vec)
    }
}

/// A trait for possibly-invertible transformation.
///
/// The inversion is analogous to 1/x, in that it might
/// sometimes fail numerically. For example, scaling the
/// plane by zero can't be reversed, since information is lost.
/// 
/// An implementation of this is required for implementing `Transformation`.
pub trait MaybeInvertible: Sized + Copy {
    /// Returns `true` if the object is invertible, `false` otherwise.
    fn is_invertible(self) -> bool;

    /// Returns the inverse of an object.
    ///
    /// This method must not panic, even if the object can't be
    /// reasonably inverted.
    fn invert(self) -> Self;

    /// Inverts an object, mutating it.
    #[inline]
    fn invert_mut(&mut self) {
        *self = self.invert();
    }

    /// Attempts to compute the inverse of a transformation.
    ///
    /// Returns `Some(inverse)` if the transformation is invertible,
    /// returns `None` otherwise.
    #[inline]
    fn invert_checked(self) -> Option<Self> {
        if self.is_invertible() {
            Some(self.invert()) 
        } else { 
            None 
        }
    }
}

impl MaybeInvertible for Transform {
    #[inline]
    fn is_invertible(self) -> bool {
        // Returns whether the determinant of the matrix is zero
        self.mat[0][0] * self.mat[1][1] != self.mat[0][1] * self.mat[1][0]
    }

    #[inline]
    fn invert(self) -> Self {
        let inv_det = (
            self.mat[0][0] * self.mat[1][1] - self.mat[0][1] * self.mat[1][0]
            ).recip();
        
        Transform {
            mat: [
                [self.mat[1][1] * inv_det, -self.mat[0][1] * inv_det],
                [-self.mat[1][0] * inv_det, self.mat[0][0] * inv_det],
            ],
            vec: inv_det * Vector {
                x: self.mat[0][1] * self.vec.y - self.mat[1][1] * self.vec.x,
                y: self.mat[1][0] * self.vec.x - self.mat[0][0] * self.vec.y,
            },
        }
    }
}

/// A trait describing the transformation of some type `T`.
pub trait Transformation<T>: MaybeInvertible {
    /// Returns a transformed copy of an object.
    ///
    /// Implementors only need to implement this function,
    /// as the other have default implementations based on this.
    fn transform(self, operand: T) -> T;

    /// Transforms an object, mutating it.
    #[inline]
    fn transform_mut(self, operand: &mut T) where T: Copy {
        *operand = self.transform(*operand);
    }

    /// Transforms an object by the inverse of a transformation.
    ///
    /// Should not panic even if the transformation isn't invertible.
    #[inline]
    fn inverse_transform(self, operand: T) -> T {
        self.invert().transform(operand)
    }

    /// Transforms an object by the inverse of a transformation.
    ///
    /// Returns `Some(result)` if the transformation is invertible,
    /// returns `None` otherwise.
    #[inline]
    fn inverse_transform_checked(self, operand: T) -> Option<T> {
        if self.is_invertible() {
            Some(self.inverse_transform(operand))
        }
        else { None }
    }
}

impl Transformation<Vector> for Transform {
    #[inline]
    fn transform(self, vec: Vector) -> Vector {
        Vector {
            x: self.mat[0][0] * vec.x + self.mat[1][0] * vec.y,
            y: self.mat[0][1] * vec.x + self.mat[1][1] * vec.y,
        }
    }

    #[inline]
    fn inverse_transform(self, vec: Vector) -> Vector {
        let inv_det = (
            self.mat[0][0] * self.mat[1][1] - self.mat[0][1] * self.mat[1][0]
        ).recip();

        Vector {
            x: inv_det * (self.mat[1][1] * vec.x - self.mat[1][0] * vec.y),
            y: inv_det * (self.mat[0][0] * vec.y - self.mat[0][1] * vec.x),
        }
    }
}

impl Transformation<Point> for Transform {
    #[inline]
    fn transform(self, pt: Point) -> Point {
        (self.transform(pt.as_vec()) + self.vec).as_point()
    }

    #[inline]
    fn inverse_transform(self, pt: Point) -> Point {
        self.inverse_transform(pt.as_vec() - self.vec).as_point()
    }
}

impl Transformation<Transform> for Transform {
    #[inline]
    fn transform(self, rhs: Transform) -> Transform {
        Transform {
            mat: [[self.mat[0][0] * rhs.mat[0][0] + self.mat[1][0] * rhs.mat[0][1],
                   self.mat[0][1] * rhs.mat[0][0] + self.mat[1][1] * rhs.mat[0][1]],
                  [self.mat[0][0] * rhs.mat[1][0] + self.mat[1][0] * rhs.mat[1][1],
                   self.mat[0][1] * rhs.mat[1][0] + self.mat[1][1] * rhs.mat[1][1]]],
            vec: self.vec + self.transform(rhs.vec),
        }
    }
}

// --OPERATOR IMPLS-- \\
// in short, implements
// vector + vector = vector
// point  + vector = point
// vector + point  = point
// vector - vector = vector
// point  - vector = point
// point  - point  = vector
// vector += vector
// point  += vector
// vector -= vector
// point  -= vector
// vector * f32 = vector
// vector / f32 = vector
// vector *= f32
// vector /= f32
// -vector = vector

macro_rules! impl_elementwise {
    ($lhs:ty, $rhs:ty, $res:ident, $trai:ident, $fun:ident) => (
        impl $trai<$rhs> for $lhs {
            type Output = $res;
            #[inline]
            fn $fun(self, other: $rhs) -> $res {
                $res {
                    x: self.x.$fun(other.x),
                    y: self.y.$fun(other.y),
                }
            }
        }
    )
}

impl_elementwise!(Vector, Vector, Vector, Add, add);
impl_elementwise!(Point,  Vector, Point,  Add, add);
impl_elementwise!(Vector, Point,  Point,  Add, add);

impl_elementwise!(Vector, Vector, Vector, Sub, sub);
impl_elementwise!(Point,  Vector, Point,  Sub, sub);
impl_elementwise!(Point,  Point,  Vector, Sub, sub);

macro_rules! impl_elwise_assign {
    ($lhs:ty, $rhs:ty, $trai:ident, $fun:ident) => (
        impl $trai<$rhs> for $lhs {
            #[inline]
            fn $fun(&mut self, other: $rhs) {
                self.x.$fun(other.x);
                self.y.$fun(other.y);
            }
        }
    )
}

impl_elwise_assign!(Vector, Vector, AddAssign, add_assign);
impl_elwise_assign!(Point,  Vector, AddAssign, add_assign);
impl_elwise_assign!(Vector, Vector, SubAssign, sub_assign);
impl_elwise_assign!(Point,  Vector, SubAssign, sub_assign);

// v1 * s = v2
impl Mul<f32> for Vector {
    type Output = Vector;
    #[inline]
    fn mul(self, scalar: f32) -> Vector {
        Vector {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

// s * v1 = v2
impl Mul<Vector> for f32 {
    type Output = Vector;
    #[inline]
    fn mul(self, vec: Vector) -> Vector {
        Vector {
            x: vec.x * self,
            y: vec.y * self,
        }
    }
}

// v1 *= s
impl MulAssign<f32> for Vector {
    #[inline]
    fn mul_assign(&mut self, scalar: f32) {
        self.x *= scalar;
        self.y *= scalar;
    }
}

// v1 / s = v2
impl Div<f32> for Vector {
    type Output = Vector;
    #[inline]
    fn div(self, scalar: f32) -> Vector {
        Vector {
            x: self.x / scalar,
            y: self.y / scalar,
        }
    }
}

// v1 /= s
impl DivAssign<f32> for Vector {
    #[inline]
    fn div_assign(&mut self, scalar: f32) {
        self.x /= scalar;
        self.y /= scalar;
    }
}

// -v1 = v2
impl Neg for Vector {
    type Output = Vector;
    #[inline]
    fn neg(self) -> Vector {
        Vector {
            x: -self.x,
            y: -self.y,
        }
    }
}
