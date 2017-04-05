//! A small collection of numeric functions.

/// Clamps `num` to the range [`a`, `b`].
///
/// If `num` is in the range [`a`, `b`], returns `num`.
/// If `num` is outside the range, returns the nearest end of the range.
///
/// The ends of the range are interchangeable, so that
/// `clamp(num, a, b) == clamp(num, b, a)`
pub fn clamp(num: f32, a: f32, b: f32) -> f32 {
    let min = a.min(b);
    let max = a.max(b);

    num.max(min).min(max)
}

/// Linearly interpolates from `start` to `end` by `t`.
///
/// `t` < 0 or `t` > 1 gives results outside [`start`, `end`].
pub fn lerp(start: f32, end: f32, t: f32) -> f32 {
    start + (end - start) * t
}

/// Linearly inpterpolates from `start` to `end` by `t`.
///
/// `t` is clamped to the range [0, 1].
pub fn lerp_clamped(start: f32, end: f32, t: f32) -> f32 {
    lerp(t.max(0.0).min(1.0), start, end)
}

/// Maps a value relative to some range to
/// the corresponding value relative to another range.
///
/// Computes how far `x` is (percentage-wise) from the start to the end 
/// of the range [`from_start`, `from_end`], and returns the value with
/// the same percentage in the range [`to_start`, `to_end`].
///
/// Still works if `x` is outside the range.
pub fn range_map(x: f32, from_start: f32, from_end: f32, to_start: f32, to_end: f32) -> f32 {
    let offset = x - from_start;
    let from_span = from_end - from_start;
    let to_span = to_end - to_start;

    (offset / from_span) * to_span + to_start
}

/// Calculates `num` modulo `denom`.
///
/// Equivalent to num % denom, except that the result
/// is corrected to never be negative. No guarantees made
/// about the result if `denom <= 0`, however.
///
/// This gives the property that (discounting rounding errors)
/// `modulo(n, d) == modulo(n + d, d)` for all `n` and positive `d`.
///
/// Returns NaN if `denom` is zero.
#[inline]
pub fn modulo(num: f32, denom: f32) -> f32 {
    let rem = num % denom;
    if rem < 0.0 {
        rem + denom
    }
    else { 
        rem
    }
}
