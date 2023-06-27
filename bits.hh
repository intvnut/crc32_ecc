#ifndef BITS_HH_
#define BITS_HH_
// Copyright 2019, J. Zbiciak <joe.zbiciak@leftturnonly.info>
// Author:  Joe Zbiciak <joe.zbiciak@leftturnonly.info>
// SPDX-License-Identifier:  CC-BY-SA-4.0

// C++17 version of my generalized CRC.
static_assert(__cplusplus >= 201703L, "Requires C++17 or later.");

#include <cstdint>
#include <limits>
#include <type_traits>

namespace bits {

// Returns the fastest type that can hold at least 'Bits' bits.
template <int Bits> struct UintFast {
  using type = typename UintFast<Bits + 1>::type;
};

template <> struct UintFast<64> { using type = uint_fast64_t; };
template <> struct UintFast<32> { using type = uint_fast32_t; };
template <> struct UintFast<16> { using type = uint_fast16_t; };
template <> struct UintFast<8>  { using type = uint_fast8_t;  };

// Returns the most closely-sized type that can hold at least 'Bits' bits.
template <int Bits> struct UintLeast {
  using type = typename UintLeast<Bits + 1>::type;
};

template <> struct UintLeast<64> { using type = uint64_t; };
template <> struct UintLeast<32> { using type = uint32_t; };
template <> struct UintLeast<16> { using type = uint16_t; };
template <> struct UintLeast<8>  { using type = uint8_t;  };

// Shorthand aliases for uint_fast<>::type and uint_size<>::type
template <int Bits> using UintFastT  = typename UintFast<Bits>::type;
template <int Bits> using UintLeastT = typename UintLeast<Bits>::type;

// Returns the maximum unsigned type we currently support.
using UintMaxT = UintLeastT<64>;

// Returns the number of bits in a type in its unsigned form.
template <typename T>
  constexpr auto TypeBits =
    std::numeric_limits<std::make_unsigned_t<T>>::digits;

// Computes mu bitmasks constants from Knuth, TAOCP Volume 4.
//
//  Mu<0,64> = 0x5555555555555555ull  Mu<0,32> = 0x55555555ull
//  Mu<1,64> = 0x3333333333333333ull  Mu<1,32> = 0x33333333ull
//                     ...                           ...
//  Mu<4,64> = 0x0000FFFF0000FFFFull  Mu<4,32> = 0x0000FFFFull
//  Mu<5,64> = 0x00000000FFFFFFFFull  Mu<5,32> = 0xFFFFFFFFull
//
// Note that Mu<K,Bits> isn't useful for K >= Log2(Bits), as it becomes
// all 1s.  But, due to how Mu<> is defined recursively, those definitions
// are available if you ask for them.
template <int K, int Bits, typename T = UintFastT<Bits>>
constexpr T kMu = (T((((Bits - 1) >> K) & 1) ^ 1) << (Bits - 1)) 
                | kMu<K, Bits - 1, T>;

template <int K, typename T> constexpr T kMu<K,0,T> = T(1);

template <int Bits, typename T = UintFastT<Bits>>
constexpr T Mu(const int k) {
  switch (k) {
    case 0: return kMu<0, Bits>;
    case 1: return kMu<1, Bits>;
    case 2: return kMu<2, Bits>;
    case 3: return kMu<3, Bits>;
    case 4: return kMu<4, Bits>;
    case 5: return kMu<5, Bits>;
    case 6: return kMu<6, Bits>;
    case 7: return kMu<7, Bits>;
    default: return k > 0 ? kMu<8, Bits> : 0;
  }
}

// Smears the bits of N rightward.
template <UintMaxT N>
  constexpr auto kSmearRight = N | kSmearRight<N / 2>;
template <> constexpr auto kSmearRight<0> = 0;

template <typename T, unsigned Bits = TypeBits<T>>
constexpr T SmearRight(const T n) {
  using U = std::make_unsigned_t<T>;
  auto result = U(n);
  auto shift = Bits / 2;
  while (shift) {
    n |= n >> shift;
    shift /= 2;
  }
  return T(result);
}

// Computes the nearest power of 2 >= N.
template <UintMaxT N>
  constexpr auto kNearestGeqPow2 = kSmearRight<N - 1> + 1;

template <typename T>
constexpr T NearestGeqPow2(const T n) {
  return SmearRight(n - 1) + 1;
}

// Computes ceil(Log2(N)); map Log2(0) to 0.
// Computed based on the unsigned image of the number.
//
// Example:
//  kLog2<2> == Log2(2) == 1
//  kLog2<3> == Log2(3) == 2
//  kLog2<4> == Log2(4) == 2
template <UintMaxT N>
  constexpr auto kLog2 = 1 + kLog2<kNearestGeqPow2<N> / 2>;
template <> constexpr auto kLog2<1> = 0;
template <> constexpr auto kLog2<0> = 0;

template <typename T, unsigned Bits = TypeBits<T>>
constexpr int Log2(const T n) {
  using U = std::make_unsigned_t<T>;
  auto un = U(n);  

  // Currently only supports up to 64 bits.
  static_assert(Bits <= 64);

  // Compute the leftmost 1 bit of (n - 1) in log2 steps, and then add 1.
  auto result = 0;

  if (un <= 1) return 0;  // Log2(0) == Log2(1) == 0.
  un -= 1;
  if (un & ~kMu<7, Bits>) { un &= ~kMu<7, Bits>; result |= U(1) << 7; }
  if (un & ~kMu<6, Bits>) { un &= ~kMu<6, Bits>; result |= U(1) << 6; }
  if (un & ~kMu<5, Bits>) { un &= ~kMu<5, Bits>; result |= U(1) << 5; }
  if (un & ~kMu<4, Bits>) { un &= ~kMu<4, Bits>; result |= U(1) << 4; }
  if (un & ~kMu<3, Bits>) { un &= ~kMu<3, Bits>; result |= U(1) << 3; }
  if (un & ~kMu<2, Bits>) { un &= ~kMu<2, Bits>; result |= U(1) << 2; }
  if (un & ~kMu<1, Bits>) { un &= ~kMu<1, Bits>; result |= U(1) << 1; }
  if (un & ~kMu<0, Bits>) { un &= ~kMu<0, Bits>; result |= U(1) << 0; }

  return result + 1;
}

// Computes a right-aligned LEN-bit bitmask. (bit 0 to LEN-1, aka. "right mask")
//
// Constant form tries to pick smallest appropriate unsigned type, unless
// overridden.
//
// Function form must be provided a type as a template arg.
//
// Example:
//   kRMask<13>          == uint16_t{0b0001'1111'1111'1111}
//   RMask<uint16_t>(13) == uint16_t{0b0001'1111'1111'1111}
template <unsigned LEN, typename T = UintLeastT<LEN>> 
  constexpr T kRMask = T{2} * kRMask<LEN - 1> + 1;
template <typename T> constexpr auto kRMask<1, T> = UintLeastT<1>{1};
template <typename T> constexpr auto kRMask<0, T> = UintLeastT<1>{0};

template <typename T>
constexpr T RMask(const unsigned len) {
  using U = std::make_unsigned_t<T>;
  if (!len) return 0;
  return T((U{2} << (len - 1)) - 1);
}

// Computes a general bitmask of length LEN, with rightmost bit at LSB.
//
// Constant form tries to pick smallest appropriate unsigned type, unless
// overridden.
//
// Function form must be provided a type as a template arg.
//
// Example:
//   kBitMask<13,2>           == uint16_t{0b0111'1111'1111'1100}
//   BitMask<uint16_t>(13, 2) == uint16_t{0b0111'1111'1111'1100}
template <unsigned LEN, unsigned LSB, typename T = UintLeastT<LEN + LSB>>
  constexpr T kBitMask = UintLeastT<LEN + LSB>{kRMask<LEN>} << LSB;

template <typename T>
constexpr T BitMask(const unsigned len, const unsigned lsb) {
  using U = std::make_unsigned_t<T>;
  if (!len) return 0;
  return T(((U{2} << (len - 1)) - 1) << lsb);
}

// Computes a left-aligned bitmask of length LEN in word of field WIDTH.
//
// Constant form tries to pick smallest appropriate unsigned type, unless
// overridden.
//
// Function form must be provided a type as a template arg.  Function form
// field widths are limited to what's expressable in a C++ fundamental type.
//
// Example:
//   kLMask<13,16>       == uint16_t{0b1111'1111'1111'1000}
//   kLMask<10,13>       == uint16_t{0b0001'1111'1111'1000}
//   LMask<uint16_t>(16) == uint16_t{0b1111'1111'1111'1000}
template <unsigned LEN, unsigned WIDTH, typename T = UintLeastT<WIDTH>>
  constexpr T kLMask = UintLeastT<WIDTH>{kRMask<LEN>} << (WIDTH - LEN);

template <typename T>
constexpr T LMask(const unsigned len) {
  using U = std::make_unsigned_t<T>;
  constexpr auto bits = TypeBits<U>;
  if (!len) return 0;
  auto nlen = bits - len;
  return T(~((U{2} << (nlen - 1)) - 1));
}

// Bit-reverses an N-bit value quickly, via a Log2-structured method.
//
// Although Hacker's Delight and TAOCP offer fancier looking algorithms
// which presumably require fewer steps, benchmarking with a recent
// compiler on AMD and Intel x86 processors suggests the Log2-structured
// method actually _performs_ the best, by a decent margin.
template <int Bits, typename ValType = UintFastT<Bits>,
          int BitsPow2 = kNearestGeqPow2<Bits>, int Rank = kLog2<Bits> - 1>
constexpr ValType BitRev(ValType val) noexcept {
  if constexpr (Rank < 0) {
    return val >> (BitsPow2 - Bits);
  } else {
    const ValType rev = (val &  kMu<Rank, BitsPow2>) << (1 << Rank)
                      | (val & ~kMu<Rank, BitsPow2>) >> (1 << Rank);
    return BitRev<Bits, ValType, BitsPow2, Rank - 1>(rev);
  }
}

}  // namespace bits
#endif // BITS_HH_
