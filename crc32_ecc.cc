#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

#include "crc.hh"

// This program demonstrates a simple ECC using a 32-bit CRC on blocks of
// 1024 bits:  992 payload bits, 32 CRC bits.  This uses CRC32/8 from the
// paper "Optimization of Cyclic Redundancy Check Codes with 24 and 32
// Parity Bits", by Castagnoli, Brauer, and Hermann (CBH).
//
// The CRC32/8 polynomial has Hamming distance 8 for blocks up to 1024 bits.
// We seek to detect up to 3 bit errors, and that requires a minimum distance
// of 7.  The minimum distance of 8 allows us to reliably detect that there's
// a 4th bit error, if I'm not mistaken, but not uniquely identify the 4
// errors.
//
// Mathematically, this makes sense. There are nearly 2^30 different
// combinations of single, dual, and triple bit errors, as each error has a
// 10-bit bit-number within the block.  A 32-bit CRC syndrome can't uniquely
// identify four bit error positions.
//
// Bit-numbering:  The code below treats the codeword as a 1024-bit integer,
// with the CRC-32 itself in bits 31-0.  Thus, each bit within the codeword
// represents a distinct power of the CRC-32 polynomial equal to its bit
// position.  Actual memory layout and transmission order put the CRC-32 last.

using std::uint16_t;
using std::uint32_t;

namespace {

using ecc_code_word = std::array<uint32_t, 32>;
using ecc_crc = crc::Crc<crc::cbh::Crc32_8, crc::LeftShifting, crc::NoInvert>;

// Gets a bit at a particular position in an ECW.
bool get_ecw_bit(const ecc_code_word& ecw, const int bit_num) {
  const int windex = 31 - (bit_num >> 5);
  const int bindex = bit_num & 31;
  return (ecw[windex] >> bindex) & 1;
}

// Sets a bit at a particular position in an ECW to a particular value.
void set_ecw_bit(ecc_code_word& ecw, const int bit_num, bool bit) {
  const int windex = 31 - (bit_num >> 5);
  const int bindex = bit_num & 31;
  const uint32_t mask = uint32_t(1) << bindex;
  const uint32_t sbit = uint32_t(bit) << bindex;
  ecw[windex] = (ecw[windex] & ~mask) | sbit;
}

// Toggles a bit at a particular position in an ECW to a particular value.
void toggle_ecw_bit(ecc_code_word& ecw, const int bit_num) {
  const int windex = 31 - (bit_num >> 5);
  const int bindex = bit_num & 31;
  const uint32_t mask = uint32_t(1) << bindex;
  ecw[windex] = ecw[windex] ^ mask;
}

// Returns the CRC word / syndrome.
uint32_t get_crc_word(const ecc_code_word& ecw) {
  return ecw[31];
}

// Sets the CRC word to a particular value.
void set_crc_word(ecc_code_word& ecw, const uint32_t crc) {
  ecw[31] = crc;
}

// XORs the CRC word with a particular value.
void xor_crc_word(ecc_code_word& ecw, const uint32_t crc) {
  ecw[31] ^= crc;
}

// Computes the CRC word for an ECW.
uint32_t compute_crc(const ecc_code_word& ecw) {
  return ecc_crc().Update(ecw.cbegin(), ecw.cend() - 1).GetChecksum();
}

// Associates a syndrome with up to two error bit positions.  The first
// error bit position is always the smaller one.  A zero in the second entry
// means this syndrome represents a single-bit error.
struct synd_to_bits {
  synd_to_bits() : syndrome{0}, bit_0{0}, bit_1{0} {}
  synd_to_bits(uint32_t s, uint16_t b0, uint16_t b1)
  : syndrome{s}, bit_0{b0}, bit_1{b1} {}

  uint32_t syndrome;
  uint16_t bit_0;
  uint16_t bit_1;
};

// Comparison predicates for looking up syndromes.
constexpr bool operator<(const synd_to_bits& lhs, const synd_to_bits& rhs) {
  return lhs.syndrome < rhs.syndrome;
}

constexpr bool operator<(const uint32_t lhs, const synd_to_bits& rhs) {
  return lhs < rhs.syndrome;
}

constexpr bool operator<(const synd_to_bits& lhs, const uint32_t rhs) {
  return lhs.syndrome < rhs;
}

constexpr bool operator==(const synd_to_bits& lhs, const synd_to_bits& rhs) {
  return lhs.syndrome == rhs.syndrome;
}

constexpr bool operator==(const synd_to_bits& lhs, const uint32_t rhs) {
  return lhs.syndrome == rhs;
}

constexpr bool operator==(const uint32_t lhs, const synd_to_bits& rhs) {
  return lhs == rhs.syndrome;
}


// synd1_bit:    Single-bit error syndromes, sorted in terms of syndrome.
// bit_to_synd:  Syndromes associated with each bit.
const auto [synd_1bit, bit_to_synd] = []{
  std::array<uint32_t, 1024> b2s{};
  std::array<synd_to_bits, 1024> s{};
  ecc_code_word ecw{};

  for (int i = 0; i < 1024; ++i) {
    toggle_ecw_bit(ecw, i);
    const uint32_t crc = compute_crc(ecw);
    xor_crc_word(ecw, crc);

    b2s[i] = get_crc_word(ecw);
    s[i] = synd_to_bits(get_crc_word(ecw), i, 0);

    xor_crc_word(ecw, crc);
    toggle_ecw_bit(ecw, i);
  }

  std::sort(s.begin(), s.end());

  return std::pair{s, b2s};
}();

// Dual-bit error syndromes, sorted.
const auto synd_2bit = []{
  std::array<synd_to_bits, 1024*1023 / 2> s{};
  ecc_code_word ecw{};
  int k = 0;

  for (int i = 0; i < 1023; ++i) {
    toggle_ecw_bit(ecw, i);
    for (int j = i + 1; j < 1024; ++j) {
      toggle_ecw_bit(ecw, j);
      const uint32_t crc = compute_crc(ecw);
      xor_crc_word(ecw, crc);
      
      s[k++] = synd_to_bits(get_crc_word(ecw), i, j);

      xor_crc_word(ecw, crc);
      toggle_ecw_bit(ecw, j);
    }
    toggle_ecw_bit(ecw, i);
  }

  std::sort(s.begin(), s.end());

  return s;
}();

// Holds error info determined by check_ecw().
struct error_info {
  error_info()
  : no_errors{true}, uncorrectable{false}, bit_0{0}, bit_1{0}, bit_2{0} {}

  explicit error_info(bool u, uint16_t b0 = 0, uint16_t b1 = 0, uint16_t b2 = 0)
  : no_errors{false}, uncorrectable{u}, bit_0{b0}, bit_1{b1}, bit_2{b2} {}

  bool     no_errors;
  bool     uncorrectable;
  uint16_t bit_0;
  uint16_t bit_1;
  uint16_t bit_2;
};

// Prints an error_info.
std::ostream& operator<<(std::ostream& os, const error_info& ei) {
  os << "ei("
     << (ei.no_errors     ? "noerr," : "err,")
     << (ei.uncorrectable ? "uncor," : "cor,")
     << std::hex << ei.bit_0 << ',' << ei.bit_1 << ',' << ei.bit_2 << ')'
     << std::dec;
  return os;
}

error_info check_ecw(const ecc_code_word& ecw) {
  const uint32_t crc = compute_crc(ecw);
  const uint32_t synd = crc ^ get_crc_word(ecw);

  if (!synd) {
    return error_info{};  // no errors.
  }

  // Check for single-bit errors.
  auto [f1, l1] = std::equal_range(synd_1bit.cbegin(), synd_1bit.cend(), synd);
  if (f1 != synd_1bit.cend() && f1->syndrome == synd) {
    return error_info{false, f1->bit_0};
  }

  // Check for double-bit errors.
  auto [f2, l2] = std::equal_range(synd_2bit.cbegin(), synd_2bit.cend(), synd);
  if (f2 != synd_2bit.cend() && f2->syndrome == synd) {
    return error_info{false, f2->bit_0, f2->bit_1};
  }
    
  // Check for triple-bit errors.
  for (int b2 = 1023; b2 >= 2 ; --b2) {
    const uint32_t b2synd = synd ^ bit_to_synd[b2];
    auto [f3, l3] = std::equal_range(synd_2bit.cbegin(), synd_2bit.cend(),
                                     b2synd);
    if (f3 != synd_2bit.cend() && f3->syndrome == b2synd) {
      return error_info{false, f3->bit_0, f3->bit_1, uint16_t(b2)};
    }
  }

  return error_info{true};  // Uncorrectable.
}

// Random number generation.
std::random_device rand_dev;
std::default_random_engine rng{rand_dev()};
std::uniform_int_distribution<unsigned> uniform_unsigned_dist{};

// Generates a random ECW and sets its CRC.
ecc_code_word create_random_ecw() {
  ecc_code_word ecw;

  for (auto& w : ecw) {
    w = uniform_unsigned_dist(rng);
  }

  set_crc_word(ecw, compute_crc(ecw));
  return ecw;
}

// Tests 0, 1, 2, 3, 4 bit errors with a given ECW.
bool test_ecc(ecc_code_word ecw) {
  
  // Zero errors.
  error_info zero_error = check_ecw(ecw);
  if (zero_error.no_errors == false) {
    std::cout << "no err fail: " << zero_error << " crc: " << std::hex
              << get_crc_word(ecw) << '\n';
    return false;
  }

  // One error, exhaustive.
  for (int b0 = 0; b0 < 1024; ++b0) {
    toggle_ecw_bit(ecw, b0);
    error_info one_error = check_ecw(ecw);
    if (one_error.no_errors || one_error.uncorrectable ||
        one_error.bit_0 != b0) {
      std::cout << "one err fail: " << one_error << " crc: " << std::hex
                << get_crc_word(ecw) << " b=" << b0 << '\n' << std::dec;
      return false;
    }
    toggle_ecw_bit(ecw, b0);
  }

  // Two errors, exhaustive.
  for (int b0 = 0; b0 < 1023; ++b0) {
    toggle_ecw_bit(ecw, b0);
    for (int b1 = b0 + 1; b1 < 1024; ++b1) {
      toggle_ecw_bit(ecw, b1);
      error_info two_error = check_ecw(ecw);
      if (two_error.no_errors || two_error.uncorrectable ||
          two_error.bit_0 != b0 || two_error.bit_1 != b1) {
        std::cout << "two err fail: " << two_error << " crc: " << std::hex
                  << get_crc_word(ecw) << " b=" << b0 << ',' << b1 << '\n'
                  << std::dec;
        return false;
      }
      toggle_ecw_bit(ecw, b1);
    }
    toggle_ecw_bit(ecw, b0);
  }

  // Three errors, non-exhaustive.  There are ~2^30 combinations and it takes
  // too long to check them all.
  std::uniform_int_distribution<int> bit_num_dist{0, 1021};
  for (int i = 0; i < 100000; ++i) {
    int b0 = bit_num_dist(rng);
    int b1 = bit_num_dist(rng);
    int b2 = bit_num_dist(rng);

    if (b0 > b1) std::swap(b0, b1);
    if (b0 > b2) std::swap(b0, b2);
    if (b1 > b2) std::swap(b1, b2);

    if (b1 >= b0) b1++;
    if (b2 >= b0) b2++;
    if (b2 >= b1) b2++;

    toggle_ecw_bit(ecw, b0);
    toggle_ecw_bit(ecw, b1);
    toggle_ecw_bit(ecw, b2);

    error_info three_error = check_ecw(ecw);
    if (three_error.no_errors || three_error.uncorrectable ||
        three_error.bit_0 != b0 || three_error.bit_1 != b1 ||
        three_error.bit_2 != b2) {
      std::cout << "three err fail: " << three_error << " crc: " << std::hex
                << get_crc_word(ecw) << " b=" << b0 << ',' << b1 << ','
                << b2 << '\n' << std::dec;
      return false;
    }

    toggle_ecw_bit(ecw, b2);
    toggle_ecw_bit(ecw, b1);
    toggle_ecw_bit(ecw, b0);
  }

  // Four errors, non-exhaustive.  There are ~2^40 combinations and it takes
  // far, far too long to check them all.
  for (int i = 0; i < 100; ++i) {
    int b0 = bit_num_dist(rng);
    int b1 = bit_num_dist(rng);
    int b2 = bit_num_dist(rng);

    if (b0 > b1) std::swap(b0, b1);
    if (b0 > b2) std::swap(b0, b2);
    if (b1 > b2) std::swap(b1, b2);

    if (b1 >= b0) b1++;
    if (b2 >= b0) b2++;
    if (b2 >= b1) b2++;

    toggle_ecw_bit(ecw, b0);
    toggle_ecw_bit(ecw, b1);
    toggle_ecw_bit(ecw, b2);

    for (int b3 = 0; b3 < 1024; ++b3) {
      if (b3 == b0 || b3 == b1 || b3 == b2) continue;

      toggle_ecw_bit(ecw, b3);
      error_info four_error = check_ecw(ecw);
      if (four_error.no_errors || !four_error.uncorrectable) {
        std::cout << "four err fail: " << four_error << " crc: " << std::hex
                  << get_crc_word(ecw) << " b=" << b0 << ',' << b1 << ','
                  << b2 << ',' << b3 << '\n' << std::dec;
        return false;
      }
      toggle_ecw_bit(ecw, b3);
    }

    toggle_ecw_bit(ecw, b2);
    toggle_ecw_bit(ecw, b1);
    toggle_ecw_bit(ecw, b0);
  }

  return true;
}

}  // namespace

int main() {
  // Test.
  for (int i = 0; i < 100; ++i) {
    std::cout << "Test #" << i << std::endl;
    if (!test_ecc(create_random_ecw())) {
      std::cout << "FAIL\n";
      return 1;
    }
  }

  std::cout << "PASS\n";

  return 0;
}
