//! Functions for pseudo-random number generation.
//!
//! All the magic happens in `impl Rng`, I'd suggest
//! looking there for specific documentation.
//!
//! **DISCLAIMER:** The functions contained within are
//! not even close to cryptographically secure.
//!
//! As a result, they are much faster than they otherwise would be
//! but sacrifice perfect distribution and entropy. Should be good
//! enough for most use cases only needing to *look* random,
//! but not for anything requiring actual, proper randomness.

use std::mem;
use rand;

use vec2d::{Vector, Point};

/// Converts a u32 to a float in [min, max).
///
/// Helper function, as described by Saito & Matsumoto (2009).
/// dx.doi.org/10.1007/978-3-642-04107-5_38
#[inline]
fn u32_to_float(uint: u32, min: f32, max: f32) -> f32 {
    // 23 bits
    let mantissa = uint >> 9;
    // sign=0, exponent=127 (bit pattern of exponent is offset by +127)
    let bit_pattern = mantissa | 0b0_01111111__0000000_00000000_00000000;
    let float_1_to_2: f32 = unsafe { mem::transmute(bit_pattern) };
    (float_1_to_2 - 1.0) * (max - min) + min
}

/// A pseudorandom number generator with a 128-bit seed.
///
/// Manually setting all bits to zero generates a sequence of zeroes.
///
/// An xorshift128+ generator, as described
/// in (Vigna, 2014, arXiv:1404.0390v3 [cs.DS]).
/// **Not cryptographically secure**, but very fast,
/// and good enough to look random.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rng(pub u64, pub u64);

impl Rng {
    /// Constructs a new `Rng`, getting the seed from the OS.
    #[inline]
    pub fn new() -> Rng {
        use rand::Rng as Rand; // for access to thread_rng
        let (s0, s1) = rand::thread_rng().gen();
        Rng::new_seeded(s0, s1)
    }

    /// Constructs a new `Rng` from predetermined seeds.
    ///
    /// If the seed is all zeroes, sets it to another fixed seed.
    #[inline]
    pub fn new_seeded(seed_a: u64, seed_b: u64) -> Rng {
        match (seed_a, seed_b) {
            (0, 0) => Rng(0x988bafaf1dc4899c, 0x330ef5df7a448957),
            _ => Rng(seed_a, seed_b),
        }
    }

    /// Generates 64 pseudo-random bits and sets the next seed.
    ///
    /// All other RNG functions on `Rng` use this function for randomness.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.0.wrapping_add(self.1);
        self.skip();
        result
    }

    /// Sets the next seed without generating a result.
    #[inline]
    pub fn skip(&mut self) {
        let mut s1 = self.0;
        let s0 = self.1;
        self.0 = s0;
        s1 ^= s1 << 23;
        self.1 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    }

    /// Moves the sequence ahead 2^64 steps.
    ///
    /// Equivalent to calling `skip()` or some other RNG function
    /// 2^64 times, but obviously faster.
    /// Useful for deterministically generating a new seed
    /// from an existing one, with practically no risk of
    /// the sequences overlapping.
    ///
    /// This algorithm is also described in the xorshift128+
    /// paper referenced in `next_u64()`.
    pub fn jump(&mut self) {
        const BIT_MASK: [u64; 2] = [0x8a5cd789635d2dff, 0x121fd2155c472f96];

        let mut s0 = 0;
        let mut s1 = 0;
        for bits in BIT_MASK.iter() {
            for bit in 0..64 {
                if bits & (1 << bit) != 0 {
                    s0 ^= self.0;
                    s1 ^= self.1;
                }
                self.skip();
            }
        }

        self.0 = s0;
        self.1 = s1;
    }

    /// Generates an integer in the range [`min`, `max`).
    ///
    /// Panics if min <= max.
    #[inline]
    pub fn rand_int(&mut self, min: i64, max: i64) -> i64 {
        assert!(max > min, "The maximum must be smaller than the minimum!");

        let delta = (max.wrapping_sub(min)) as u64;
        let random = self.next_u64();
        min + (random % delta) as i64
    }

    /// Generates an unsigned integer in the range [`min`, `max`).
    #[inline]
    pub fn rand_uint(&mut self, min: u64, max: u64) -> u64 {
        assert!(max > min, "The maximum must be smaller than the minimum!");

        let delta = max - min;
        let random = self.next_u64();
        min + random % delta
    }

    /// Generates a float in the range [`min`, `max`).
    ///
    /// Uses a technique described by Saito & Matsumoto (2009).
    /// dx.doi.org/10.1007/978-3-642-04107-5_38
    #[inline]
    pub fn rand_float(&mut self, min: f32, max: f32) -> f32 {
        u32_to_float(self.next_u64() as u32, min, max)
    }

    /// Has a `probability` chance of returning `true`.
    ///
    /// `probability` is effectively clamped to the range [0, 1],
    /// so p >= 1 is guaranteed `true` and p <= 0 is guaranteed `false`.
    #[inline]
    pub fn rand_bool(&mut self, probability: f32) -> bool {
        self.rand_float(0.0, 1.0) < probability
    }

    /// Generates a Vector within the (inclusive) bounds of a rectangle.
    ///
    /// The rectangle's bounds are defined by two opposite corners.
    ///
    /// Since a Vector is made up of two components, the RNG sequence
    /// is moved forward two steps.
    #[inline]
    pub fn rand_vec(&mut self, bound1: Vector, bound2: Vector) -> Vector {
        let min = [bound1.x.min(bound2.x), bound1.y.min(bound2.y)];
        let max = [bound1.x.max(bound2.x), bound1.y.max(bound2.y)];

        let random = self.next_u64();
        let split: [u32; 2] = unsafe { mem::transmute(random) };

        Vector::new(u32_to_float(split[0], min[0], max[0]),
                    u32_to_float(split[1], min[1], max[1]))
    }

    /// Generates a Point within the (inclusive) bounds of a rectangle.
    ///
    /// The rectangle's bounds are defined by two opposite corners.
    ///
    /// Since a Point is made up of two components, the RNG sequence
    /// is moved forward two steps.
    #[inline]
    pub fn rand_point(&mut self, bound1: Point, bound2: Point) -> Point {
        self.rand_vec(bound1.as_vec(), bound2.as_vec()).as_point()
    }

    /// Chooses an element in an array slice.
    ///
    /// Panics if the slice's length is zero.
    #[inline]
    pub fn rand_el<'a, T>(&mut self, arr: &'a [T]) -> &'a T {
        assert!(arr.len() > 0, "A random element? There are no elements!");
        let random = self.next_u64() as usize;
        let index = random % arr.len();
        &arr[index]
    }

    /// Gets a `mut` pointer to a randomly chosen element in an array slice.
    ///
    /// Panics if the slice's length is zero.
    #[inline]
    pub fn rand_el_mut<'a, T>(&mut self, arr: &'a mut [T]) -> &'a mut T {
        assert!(arr.len() > 0, "A random element? There are no elements!");
        let random = self.next_u64() as usize;
        let index = random % arr.len();
        &mut arr[index]
    }

    /// Shuffles the elements of an array into a random order.
    ///
    /// Assumes (fairly safely) that the length of the array fits in a u64.
    #[inline]
    pub fn shuffle_mut<T>(&mut self, arr: &mut [T]) {
        for i in (1..arr.len()).rev() {
            let j = self.rand_uint(0, i as u64 + 1) as usize;
            arr.swap(i, j);
        }
    }

    /// Creates a new `Vec<T>` with the elements of another array in random order.
    ///
    /// Clones each element in the given array to create a new one.
    /// To avoid allocation and potentially-expensive `.clone()` calls,
    /// use `shuffle_mut()`.
    ///
    /// Assumes (fairly safely) that the length of the array fits in a u64.
    #[inline]
    pub fn shuffle<T: Clone>(&mut self, arr: &[T]) -> Vec<T> {
        let mut out = Vec::from(arr);
        self.shuffle_mut(&mut out);
        out
    }
}
