#ifndef CRC_HH_
#define CRC_HH_
// Copyright 2019, J. Zbiciak <joe.zbiciak@leftturnonly.info>
// Author:  Joe Zbiciak <joe.zbiciak@leftturnonly.info>
// SPDX-License-Identifier:  CC-BY-SA-4.0

// C++17 version of my generalized CRC.
static_assert(__cplusplus >= 201703L, "Requires C++17 or later.");

#include <cstdint>
#include <limits>
#include <type_traits>
#include "bits.hh"

namespace crc {
//  CRC Calculation Template Class
//
//  Notation follows RockSoft's definition.
//
//  RockSoft's model is a big endian model, where the MSB of both the
//  CRC and the data represent the earliest bit in the bitstream.  This
//  fits the view of CRC as polynomial division, with the MSB of the CRC
//  corresponding to the highest powered value in the field.
//
//  The RefIn flag indicates that the LSB of each datum is the 'earliest'
//  bit, rather than the MSB.
//
//  The RefOut flag indicates that the final CRC should be bit-reversed,
//  so that the 'earliest' bit of the CRC ends up in the LSB rather than
//  the MSB.
//
//  The truncated polynomial assumes a big-endian representation.
//
//  The "CRC-32" algorithm, for example, sets RefIn and RefOut both to 'true'.
//  It's a 'right-shifting' CRC, where the 'earliest' bit of both the data
//  and the checksum are in the LSB.  The Rocksoft model, however, represeents
//  this as a reflected input and reflected output.  The polynomial remains
//  big endian.  As a result, you need to use the polynomial 0x04C11DB7u, not
//  0xEDB88320 as many C implementations use.

// Default CRC paramters.
struct CrcDefaults {
  constexpr static int              kWidth  = 32;
  constexpr static ::bits::UintMaxT kPoly   = 0x04C11DB7u;
  constexpr static ::bits::UintMaxT kInit   = 0xFFFFFFFFu;
  constexpr static bool             kRefIn  = true;
  constexpr static bool             kRefOut = true;
  constexpr static ::bits::UintMaxT kXorOut = 0xFFFFFFFFu;
};

// Specifies polynomial width, in bits.
template <int W> struct Width : virtual CrcDefaults {
  constexpr static int kWidth = W;
};

// Specifies truncated polynomial.  Only lowest 'width' bits are significant.
template <::bits::UintMaxT P> struct Poly : virtual CrcDefaults {
  constexpr static ::bits::UintMaxT kPoly = P;
};

// Specifies initial CRC value, if none provided at construction time.
template <::bits::UintMaxT I> struct Init : virtual CrcDefaults {
  constexpr static ::bits::UintMaxT kInit = I;
};

// Specifies whether CRC should reflect the input (in terms of bit ordering).
// false:  MSB is first bit.  true:  LSB is first.
template <bool RI> struct RefIn : virtual CrcDefaults {
  constexpr static bool kRefIn = RI;
};

// Specifies whether CRC should reflect the output (in terms of bit ordering).
// false:  MSB is first bit;  true:  LSB is first.
template <bool RO> struct RefOut : virtual CrcDefaults {
  constexpr static bool kRefOut = RO;
};

// Specifies XOR value to apply on the output at the end, before any bit
// order reflection.
template <::bits::UintMaxT XO> struct XorOut : virtual CrcDefaults {
  constexpr static auto kXorOut = XO;
};

// RightShifting: Shorthand for RefIn<true>, RefOut<true>
// LeftShifting:  Shorthand for RefIn<false>, RefOut<false;
struct RightShifting : virtual CrcDefaults {
  constexpr static bool kRefIn  = true;
  constexpr static bool kRefOut = true;
};
struct LeftShifting : virtual CrcDefaults {
  constexpr static bool kRefIn  = false;
  constexpr static bool kRefOut = false;
};

// Invert:   Shorthand for Init<all 1s>, XorOut<all 1s>
// NoInvert: Shorthand for Init<0>, XorOut<0>
struct Invert : virtual CrcDefaults {
  constexpr static ::bits::UintMaxT kInit   = ~::bits::UintMaxT{0};
  constexpr static ::bits::UintMaxT kXorOut = ~::bits::UintMaxT{0};
};
struct NoInvert : virtual CrcDefaults {
  constexpr static ::bits::UintMaxT kInit   = 0;
  constexpr static ::bits::UintMaxT kXorOut = 0;
};

// Field:  Shorthand for specifying Width and Polynomial.
template <int W, ::bits::UintMaxT P>
struct Field : virtual CrcDefaults {
  constexpr static int              kWidth = W;
  constexpr static ::bits::UintMaxT kPoly  = P;
};

// Internal helpers.
namespace internal {

// The Discriminator and ParameterSelector classes allow us to build up
// parameters piecemeal in any order, overriding the default parameters
// as needed.
//
// See C++ Templates Book, 2d Ed, Section 21.4.
template <typename B, int I> class Discriminator : public B { };

template <typename P1, typename P2, typename P3,
          typename P4, typename P5, typename P6>
struct ParameterSelector : Discriminator<P1, 1>, Discriminator<P2, 2>,
                           Discriminator<P3, 3>, Discriminator<P4, 4>,
                           Discriminator<P5, 5>, Discriminator<P6, 6> { };

struct CrcDefaultArgs : virtual CrcDefaults { };

}  // namespace internal


template<typename P1 = internal::CrcDefaultArgs,
         typename P2 = internal::CrcDefaultArgs,
         typename P3 = internal::CrcDefaultArgs,
         typename P4 = internal::CrcDefaultArgs,
         typename P5 = internal::CrcDefaultArgs,
         typename P6 = internal::CrcDefaultArgs>
class Crc {
 public:
  using Params     = ::crc::internal::ParameterSelector<P1, P2, P3, P4, P5, P6>;
  using ResultType = ::bits::UintLeastT<Params::kWidth>;
  using FastType   = ::bits::UintFastT<Params::kWidth>;

  constexpr static int      kWidth  = Params::kWidth;
  constexpr static FastType kMask   = ::bits::kRMask<Params::kWidth>;
  constexpr static FastType kPoly   = FastType(Params::kPoly & kMask);
  constexpr static FastType kInit   = FastType(Params::kInit & kMask);
  constexpr static bool     kRefIn  = Params::kRefIn;
  constexpr static bool     kRefOut = Params::kRefOut;
  constexpr static FastType kXorOut = FastType(Params::kXorOut & kMask);

  constexpr Crc() noexcept : accum_{kRefInit} {}

  explicit constexpr Crc(const ResultType init, bool reflect = kRefIn) noexcept
    : accum_{reflect ? ::bits::BitRev<kWidth>(init) : init} {}

  [[nodiscard]] constexpr ResultType GetChecksum() const noexcept {
    if constexpr (kRefIn == kRefOut) {
      return static_cast<ResultType>(accum_ ^ kRefXorOut);
    } else {
      const auto xored = accum_ ^ kRefXorOut;
      return static_cast<ResultType>(::bits::BitRev<kWidth>(xored));
    }
  }

  template <typename T,
        unsigned Bits = std::numeric_limits<std::make_unsigned_t<T>>::digits>
  [[nodiscard]] constexpr auto Update(const T value) const noexcept {
    // Note:  u_value is used on all paths, but GCC 9.2 is not convinced.
    [[maybe_unused]] const auto u_value = FastType{value};
    Crc result = *this;
    if constexpr (Bits > kWidth) {
      // Assume for now kWidth % kLutBits == 0.  It will still work otherwise;
      // however, it may not be optimal.
      const auto head_cnt = kWidth;
      const auto tail_cnt = Bits - head_cnt;
      const auto [head, tail] = SplitBits<head_cnt, Bits>(u_value);
      result = result.Update<FastType, head_cnt>(head);
      result = result.Update<FastType, tail_cnt>(tail);
    } else if constexpr (Bits >= kLutBits) {
      const auto head_cnt = Bits - (Bits > kLutBits ? Bits % kLutBits : 0);
      const auto tail_cnt = Bits - head_cnt;
      const auto [head, tail] = SplitBits<head_cnt, Bits>(u_value);
      result.accum_ = LutUpdate<head_cnt>(accum_, head);
      result = result.Update<FastType, tail_cnt>(tail);
    } else if constexpr (Bits != 0) {
      result.accum_ = IterUpdate<Bits>(accum_, u_value);
    }
    return result;
  }

  template <typename ForwardIt>
  [[nodiscard]] constexpr Crc Update(ForwardIt b, ForwardIt e) const noexcept {
    using ET = std::remove_reference_t<decltype(*b)>;
    using UT = std::make_unsigned_t<ET>;
    constexpr auto bits_per_elem = std::numeric_limits<UT>::digits;
    // Note:  Using the size of FastType seems slower than using kWidth.
    //constexpr auto word_size = std::numeric_limits<FastType>::digits;
    constexpr auto word_size = kWidth;
    constexpr auto iters_per_word = std::max(1, word_size / bits_per_elem);
    constexpr auto bits_per_word = bits_per_elem * iters_per_word;
    auto crc = *this;
    auto word = FastType{0};
    auto iter = 0u;
    
    // Update in 'word' size chunks.  Allow us to completely fill up a
    // FastType word before checksumming it.
    while (b != e) {
      if constexpr (!kRefIn) {
        if constexpr (bits_per_elem < word_size) {
          word = (word << bits_per_elem) | *b;
        } else {
          word = *b;
        }
      } else {
        word |= FastType{UT(*b)} << (iter * bits_per_elem);
      }
      ++b;
      if (++iter == iters_per_word) {
        crc = crc.template Update<FastType, bits_per_word>(word);
        word = 0;
        iter = 0;
      }
    }

    // Handle the partial word at the end one element at a time.
    if (bits_per_elem < word_size) {
      if constexpr (!kRefIn) {
        while (iter-- != 0) {
          const auto elem = (word >> (bits_per_elem * iter)) 
                          & ::bits::kRMask<bits_per_elem>;
          crc = crc.template Update<FastType, bits_per_elem>(elem);
        }
      } else {
        for (auto i = 0u; i != iter; ++i) {
          const auto elem = word & ::bits::kRMask<bits_per_elem>;
          word >>= bits_per_elem;
          crc = crc.template Update<FastType, bits_per_elem>(elem);
        }
      }
    }
    return crc;
  }

  // Debug interface.
  static const auto& GetLut() noexcept { return lut; }

 private:
  FastType accum_ = kRefInit;

  // Lookup-table size for fast implementation.  For CRCs smaller than 8
  // bits, limit the table size to the CRC width.
  constexpr static int kLutBits = std::min(kWidth, 8);
  constexpr static int kLutSize = 1 << kLutBits;
  constexpr static int kLutCount = std::min(4, std::max(1, kWidth/kLutBits));

  // Possibly reflected Poly and Init, if RefIn set.
  constexpr static FastType kRefPoly =
    kRefIn ? ::bits::BitRev<kWidth>(kPoly) : kPoly;
  constexpr static FastType kRefInit =
    kRefIn ? ::bits::BitRev<kWidth>(kInit) : kInit;

  // Possibly reflected XorOut, if RefOut set.
  constexpr static FastType kRefXorOut =
    kRefOut ? ::bits::BitRev<kWidth>(kXorOut) : kXorOut;

  // Splits the remaining TotBits bits into Head and Tail.
  // Assumes lower TotBits bits are valid of the input, regardless of RefIn.
  template <unsigned HeadBits, unsigned TotBits>
  constexpr static auto SplitBits(FastType value) noexcept {
    if constexpr (!kRefIn) {
      const auto head = HeadBits ? value >> (TotBits - HeadBits) : FastType{0};
      const auto tail = value & ::bits::kRMask<TotBits - HeadBits>;
      return std::make_pair(head, tail);
    } else {
      const auto head = value & ::bits::kRMask<HeadBits>;
      const auto tail = HeadBits < TotBits ? value >> HeadBits : FastType{0};
      return std::make_pair(head, tail);
    }
  }

  // Computes CRC iteratively, one bit at a time.
  template <unsigned Bits>
  constexpr static FastType
      IterUpdate(FastType accum, FastType value) noexcept {
    if constexpr (!kRefIn) {
      // If we're not reflecting the input, implement in a big-endian manner.
      accum ^= value << (kWidth - Bits);
      for (unsigned i = 0; i != Bits; ++i) {
        const auto bit = bool(1 & (accum >> (kWidth - 1)));
        accum <<= 1;
        if (bit) accum ^= kRefPoly;
        accum &= kMask;
      }
      return accum;
    } else {
      // If we're reflecting the input, implement in a little-endian manner.
      accum ^= value;
      for (unsigned i = 0; i != Bits; ++i) {
        const auto bit = bool(1 & accum);
        accum >>= 1;
        if (bit) accum ^= kRefPoly;
      }
      return accum;
    }
  }

  // Computes CRC with a lookup table (LUT), up to kLutBits per iteration,
  // up to a total of kWidth bits.  Bits must be a multiple of kLutBits.
  template <unsigned UpdateBits>
  constexpr static FastType LutUpdate(FastType accum, FastType value) noexcept {
    constexpr auto lut_batch   = kLutBits * kLutCount;
    constexpr auto batch_iters = UpdateBits / lut_batch;
    constexpr auto end_iters   = (UpdateBits % lut_batch) / kLutBits;
    static_assert(UpdateBits <= kWidth);
    if constexpr (!kRefIn) {
      constexpr auto idx_mask = ::bits::kRMask<kLutBits>;
      // If we're not reflecting the input, implement in a big-endian manner.
      accum ^= value << (kWidth - UpdateBits);
      for (auto iters = 0u; iters < batch_iters; ++iters) {
        if constexpr (kLutCount == 4) {
          const auto lut_idx3 = (accum >> (kWidth - 1*kLutBits)) & idx_mask;
          const auto lut_idx2 = (accum >> (kWidth - 2*kLutBits)) & idx_mask;
          const auto lut_idx1 = (accum >> (kWidth - 3*kLutBits)) & idx_mask;
          const auto lut_idx0 = (accum >> (kWidth - 4*kLutBits)) & idx_mask;
          const auto lut_val3 = lut[3][lut_idx3];
          const auto lut_val2 = lut[2][lut_idx2];
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? (accum << lut_batch) & kMask : 0;
          accum ^= (lut_val3 ^ lut_val2) ^ (lut_val1 ^ lut_val0);
        } else if constexpr (kLutCount == 3) {
          const auto lut_idx2 = (accum >> (kWidth - 1*kLutBits)) & idx_mask;
          const auto lut_idx1 = (accum >> (kWidth - 2*kLutBits)) & idx_mask;
          const auto lut_idx0 = (accum >> (kWidth - 3*kLutBits)) & idx_mask;
          const auto lut_val2 = lut[2][lut_idx2];
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? (accum << lut_batch) & kMask : 0;
          accum ^= lut_val2 ^ (lut_val1 ^ lut_val0);
        } else if constexpr (kLutCount == 2) {
          const auto lut_idx1 = (accum >> (kWidth - 1*kLutBits)) & idx_mask;
          const auto lut_idx0 = (accum >> (kWidth - 2*kLutBits)) & idx_mask;
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? (accum << lut_batch) & kMask : 0;
          accum ^= lut_val1 ^ lut_val0;
        } else {
          const auto lut_idx0 = (accum >> (kWidth - 1*kLutBits)) & idx_mask;
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? (accum << lut_batch) & kMask : 0;
          accum ^= lut_val0;
        } 
      }
      for (auto iters = 0u; iters < end_iters; ++iters) {
        const auto lut_idx = (accum >> (kWidth - kLutBits)) ;
        accum = (accum << kLutBits) & kMask;
        accum ^= lut[0][lut_idx];
      }
      return accum;
    } else {
      // If we're reflecting the input, implement in a little-endian manner.
      constexpr auto idx_mask = ::bits::kRMask<kLutBits>;
      accum ^= value;
      for (auto iters = 0u; iters < batch_iters; ++iters) {
        if constexpr (kLutCount == 4) {
          const auto lut_idx3 = (accum >> 0*kLutBits) & idx_mask;
          const auto lut_idx2 = (accum >> 1*kLutBits) & idx_mask;
          const auto lut_idx1 = (accum >> 2*kLutBits) & idx_mask;
          const auto lut_idx0 = (accum >> 3*kLutBits) & idx_mask;
          const auto lut_val3 = lut[3][lut_idx3];
          const auto lut_val2 = lut[2][lut_idx2];
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? accum >> lut_batch : 0;
          accum ^= lut_val3 ^ lut_val2 ^ (lut_val1 ^ lut_val0);
        } else if constexpr (kLutCount == 3) {
          const auto lut_idx2 = (accum >> 0*kLutBits) & idx_mask;
          const auto lut_idx1 = (accum >> 1*kLutBits) & idx_mask;
          const auto lut_idx0 = (accum >> 2*kLutBits) & idx_mask;
          const auto lut_val2 = lut[2][lut_idx2];
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? accum >> lut_batch : 0;
          accum ^= lut_val2 ^ lut_val1 ^ lut_val0;
        } else if constexpr (kLutCount == 2) {
          const auto lut_idx1 = (accum >> 0*kLutBits) & idx_mask;
          const auto lut_idx0 = (accum >> 1*kLutBits) & idx_mask;
          const auto lut_val1 = lut[1][lut_idx1];
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? accum >> lut_batch : 0;
          accum ^= lut_val1 ^ lut_val0;
        } else {
          const auto lut_idx0 = (accum >> 1*kLutBits) & idx_mask;
          const auto lut_val0 = lut[0][lut_idx0];
          accum = lut_batch < kWidth ? accum >> lut_batch : 0;
          accum ^= lut_val0;
        }
      }
      for (auto iters = 0u; iters < end_iters; ++iters) {
        const auto lut_idx = accum & idx_mask;
        accum >>= kLutBits;
        accum ^= lut[0][lut_idx];
      }
      return accum;
    }
  }

  // Generates the lookup table.
  constexpr static auto  LutGen() noexcept {
    auto lut_init = std::array<std::array<FastType, kLutSize>, kLutCount>{};
    for (FastType lut_idx = 0; lut_idx < kLutSize; lut_idx++) {
      auto xor_val = IterUpdate<kLutBits>(0, lut_idx);
      for (auto lut_num = 0; lut_num < kLutCount; lut_num++) {
        lut_init[lut_num][lut_idx] = xor_val;
        xor_val = IterUpdate<kLutBits>(xor_val, 0);
      }
    }
    return lut_init;
  }

  // Lookup-table for fast computation.
  constexpr static auto lut = LutGen();
};

// Common CRCs.
// These are complete Crc types, not parameter sets.
namespace common {
// Notes:
//  -- Crc32Ieee802 is also used by Ethernet and ZIP.
//  -- Crc16Ccitt is also used by Intellicart.
using Crc32Ieee802 = Crc<CrcDefaults>;  // AKA "zip";
using Crc16Ccitt   = Crc<Field<16, 0x1021>,
                         Init<0xFFFF>, XorOut<0>, LeftShifting>;
using Crc16Jlp     = Crc<Field<16, 0x4AB5>, NoInvert, RightShifting>;
} // namespace common

// The following parameter sets come from the paper "Optimization of Cyclic
// Redundancy Check Codes with 24 and 32 Parity Bits", by Castagnoli, Brauer,
// and Hermann (CBH).
//
// Note that the paper does not specify a pre/post inversion, nor does it
// specify reflection on input or output.  Thus, the first set of parameter
// structures merely define the field width and truncated polynomial.
//
// Combine these with RefIn/RefOut and Init/XorOut as needed.
namespace cbh {
// 24-bit fields.
//
// dmin     24/6.1        24/6.2       24/5.1       24/5.2         24/4
//------------------------------------------------------------------------------
//  16                   25 - 26
//  15                                  25         25 - 26
//  14        25
//  12        26                                   27 - 28        25 - 30
//  11                                             29 - 31
//  10     27 - 36       27 - 41      25 - 33      32 - 33        31 - 36
//   9                                             34 - 35
//   8     37 - 83       42 - 95                   36 - 41        37 - 61
//   7                                34 - 37      42 - 77
//   6     84 - 2050     96 - 2048    38 - 252     78 - 217       62 - 846
//   5                               253 - 4097   218 - 4095
//   4   2051 - 4098   2049 - 4094                               847 - 8388607
//   2      >= 4099      >= 4095       >= 4098     >= 4096        >= 8388608
//------------------------------------------------------------------------------
using Crc24_6_1 = Field<24, 0x5D6DCB>;
using Crc24_6_2 = Field<24, 0x7B01BD>;
using Crc24_5_1 = Field<24, 0x31FF19>;
using Crc24_5_2 = Field<24, 0x5BC4F5>;
using Crc24_4   = Field<24, 0x328B63>;

// 32-bit fields.
//
// dm  IEEE-802       32/8        32/6       32/5.1     32/5.2      32/4
//------------------------------------------------------------------------------
// 20                               33
// 18                             34-35                                 33
// 17                                                     33-34       34-38
// 16                               36
// 15    33-43                                33-35
// 14                 33-44         37                     35         39-40
// 13                                                     36-38
// 12    43-44        45-48       38-43       36-49                   41-52
// 11    45-53                                50-53       39-52
// 10    54-66        49-98       44-56       54-59       53-68       53-79
//  9    67-89                                            69-80
//  8    90-123       99-1024     57-306      60-90       81-110      80-209
//  7   124-203                               91-113     111-266
//  6   204-300                  307-32768   114-1092    267-1029    210-5275
//  5   301-3006                            1093-65537  1030-65535
//  4  3007-x*      1025-2046  32769-65534                          5276-2^31-1
//  3   x+1-2^32-1
//  2    >= 2^32      >= 2047     65535      >= 65538    >= 65536    >= 2^31
//------------------------------------------------------------------------------
// *x >= 64,000;  dm == minimum Hamming distance.
using Crc32_8   = Field<32, 0xF1922815>;
using Crc32_6   = Field<32, 0xF6ACFB13>;
using Crc32_5_1 = Field<32, 0xA833982B>;
using Crc32_5_2 = Field<32, 0x572D7285>;
using Crc32_4   = Field<32, 0x1EDC6F41>;
} // namespace cbh

// The following parameters sets were proposed by Merkey and Posner in
// "Optimum Cyclic Redundancy Codes for Noisy Channels."  The polynomials
// themselves, however, were referenced from in CBH's paper noted above.
//
// Note: The 32-bit codes appear to actually be 31 bit codes, given that
// the LSB is not set.
namespace mp {
using Crc24_6_1 = Field<24, 0x323009>;
using Crc24_6_2 = Field<24, 0x401607>;
using Crc32_8_1 = Field<32, 0x404098E2>;
using Crc32_8_2 = Field<32, 0x0884C512>;
} // namespace mp

} // namespace crc

#endif // CRC_HH_
