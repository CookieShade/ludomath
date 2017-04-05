//! Various relevant mathematical constants.
//!
//! The constants here in the module would be on their respective types,
//! but that [isn't supported yet in stable Rust.](https://github.com/rust-lang/rust/issues/29646)

use vec2d::{Rotation, Vector, Point, Transform};

/// The circle constant, τ. Historically known as 2π.
///
/// The ratio between the circumpherence of a circle and the radius,
/// or the number of radians in a full turn.
///
/// Seen in Euler's identity: e^iτ = 1.
pub const TAU: f32 = 6.2831853;

/// The origin of the coordinate system, (0, 0).
pub const POINT_ORIGIN: Point = Point { x: 0.0, y: 0.0 };

/// The vector (1, 1).
pub const VEC_ONE: Vector = Vector { x: 1.0, y: 1.0 };
/// The zero vector, (0, 0).
pub const VEC_ZERO: Vector = Vector { x: 0.0, y: 0.0 };

/// A unit vector pointing upwards, (0, 1).
pub const VEC_UP: Vector = Vector { x: 0.0, y: 1.0 };
/// A unit vector pointing downwards, (0, -1).
pub const VEC_DOWN: Vector = Vector { x: 0.0, y: -1.0 };
/// A unit vector pointing left, (-1, 0).
pub const VEC_LEFT: Vector = Vector { x: -1.0, y: 0.0 };
/// A unit vector pointing right, (1, 0).
pub const VEC_RIGHT: Vector = Vector { x: 1.0, y: 0.0 };

const SQRT1_2: f32 = 0.70710678;

/// A unit vector pointing up and left, at a 45 degree angle.
pub const VEC_UP_LEFT   : Vector = Vector { x: -SQRT1_2, y:  SQRT1_2 };
/// A unit vector pointing up and right, at a 45 degree angle.
pub const VEC_UP_RIGHT  : Vector = Vector { x:  SQRT1_2, y:  SQRT1_2 };
/// A unit vector pointing down and left, at a 45 degree angle.
pub const VEC_DOWN_LEFT : Vector = Vector { x:  SQRT1_2, y: -SQRT1_2 };
/// A unit vector pointing down and right, at a 45 degree angle.
pub const VEC_DOWN_RIGHT: Vector = Vector { x: -SQRT1_2, y: -SQRT1_2 };

/// A rotation of zero degrees.
pub const ROTATION_0_DEG  : Rotation = Rotation { cos: 1.0, sin: 0.0 };
/// A rotation of 45 degrees.
pub const ROTATION_45_DEG : Rotation = Rotation { cos: SQRT1_2, sin: SQRT1_2 };
/// A rotation of 90 degrees.
pub const ROTATION_90_DEG : Rotation = Rotation { cos: 0.0, sin: 1.0 };
/// A rotation of 135 degrees.
pub const ROTATION_135_DEG: Rotation = Rotation { cos: -SQRT1_2, sin: SQRT1_2 };
/// A rotation of 180 degrees.
pub const ROTATION_180_DEG: Rotation = Rotation { cos: -1.0, sin: 0.0 };
/// A rotation of 225 degrees.
pub const ROTATION_225_DEG: Rotation = Rotation { cos: -SQRT1_2, sin: -SQRT1_2 };
/// A rotation of 270 degrees.
pub const ROTATION_270_DEG: Rotation = Rotation { cos: 0.0, sin: -1.0 };
/// A rotation of 315 degrees.
pub const ROTATION_315_DEG: Rotation = Rotation { cos: SQRT1_2, sin: -SQRT1_2 };

/// A transformation which has no effect,
/// representing the global coordinate system.
pub const TRANSFORM_IDENTITY: Transform = Transform {
    vec: VEC_ZERO,
    mat: [[1.0, 0.0],
          [0.0, 1.0]],
};
