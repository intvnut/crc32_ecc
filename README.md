# ECC using CRC-32

Simple demonstration of using a 32-bit CRC as an error correcting code (ECC).

_Note:_  This code is a quick and dirty implementation.  I've made no attempt
to optimize the ECC code itself.  This is just draft quality code.

## License

Everything in this repository other than `crc_v3.txt` is authored by me
(Joe Zbiciak, joe.zbiciak@leftturnonly.info), and is licensed under the
Creative Commons Attribution-ShareAlike 4.0 International license, aka.
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

The file [`crc_v3.txt`](./crc_v3.txt) is authored by Ross N. Williams of
Rocksoft(TM) Pty Ltd.  See the file for details.

## Background

The CRC-32 code in this directory ([`crc.hh`](./crc.hh) and
[`bits.hh`](./bits.hh) is my stab at making a highly flexible, highly
`constexpr` "functional-style" CRC-32 class.  I've reused it here.

The CRC-32 ECC code in [`crc32_ecc.cc`](./crc32_ecc.cc) was inspired by a
paragraph I saw go by in an answer on Quora,
[here](https://www.quora.com/Is-it-possible-if-the-target-BER-is-10-32-or-what-can-be-the-maximum-BER-possible-if-we-want-to-use-64-bits-for-CRC/answer/Rcgldr):

> It is also possible to correct error bits in smaller messages,
> such as up to 3 error bits corrected with a specific 32 bit CRC for
> a 992 data bit + 32 bit CRC (1024 bits total), but I ended up using
> a 1.4 GB table to do this. I’m not aware of a mathematical method
> to use CRC for error correction, only for actual error correcting
> codes like Hamming or Reed Solomon.

After thinking about it a bit, I realized you could probably get away with
much smaller tables if you were willing to make up to 1024 probes in the
case of a 3-bit or 4+ bit error.

## CRC-32 Polynomial Selection

I selected the polynomial CRC32/8 from the paper "Optimization of Cyclic
Redundancy Check Codes with 24 and 32 Parity Bits", by Castagnoli, Brauer,
and Hermann (CBH).  This polynomial has Hamming distance 8 for all checksums
on blocks of 1024 bits or shorter.

This fits the parameters of the comment I quoted above.

With Hamming distance 8, we should be able to correct 1-bit, 2-bit, and 3-bit
errors, and detect 4-bit errors reliably.

I configured the CRC-32 as left-shifting with no pre/post inversion for
simplicity.  It should be possible to rework the code to support pre/post
inversion, etc., but I haven't taken the effort.

The code treats the ECC code word as a 1024-bit integer, stored as 32 x 32-bit
integers.  Bit 1023 is the MSB of the first word.  The CRC itself occupies 
the last word of the array, and corresponds to bits 31..0 of the 1024-bit
integer.

Within this large integer, the syndrome for bit #n is x^n, where x is the
polynomial.

## ECC Check/Correct Algorithm Summary

The code leverages three tables:

1. `synd_1bit` holds the syndromes associated with each single-bit error,
   and the corresponding bit number.  It is sorted by syndrome.  This table
   is 8192 bytes.

2. `bit_to_synd` holds the same syndromes as `synd_1bit`, but is indexed
   by bit number.  This table is 4096 bytes.

3. `synd_2bit` holds syndromes associated with each two-bit error, and the
   corresponding bit numbers.  Like `synd_1bit`, it is sorted by syndrome.
   This table is 4MiB - 4096 bytes.

The ECC check algorithm computes the syndrome (XOR of received CRC
with locally computed CRC).

*  If the syndrome is 0, return "no errors".
*  Binary search for the syndrome in `synd_1bit`.  If found, return the bit
   position of the single-bit error.
*  Binary search for the syndrome in `synd_2bit`.  If found, return the bit
   positions of the two errored bits.
*  For each syndrome in `bit_to_synd`, for b2 = 2..1023:
   * XOR the computed syndrome with `bit_to_synd[b2]`.
   * Binary search for the syndrome in `synd_2bit`  If found, return the bit
     positions of the three errored bits.
*  Return "uncorrectable", otherwise.

The correction step is trivial:  Just flip the indicated bits.  I have not
actually implemented that in this code.  The test-code adequately demonstrates
that the check routine finds the correct bit positions to flip.

## Potential Improvements

* If the syndrome only has 1, 2, or 3 bits set, the bit errors are in the CRC
  itself.  You can just replace the received CRC with the locally computed CRC.

* You should be able to create a perfect hashing function that will let you
  look up syndromes in `synd_1bit` and `synd_2bit` with a single probe.
  If you take advantage of the previous bullet to ignore bit errors in the CRC
  itself, then for `synd_1bit` that may be as simple as taking a prefix of
  _n_ bits from the syndrome.  

* Encapsulate the ECC code word in a class.

## Beware of Bugs!

I haven't tested this code thoroughly.  In particular, the CRC-32 code itself
hasn't been thoroughly tested.  I found and fixed a couple bugs while writing
this program.

____

Copyright © 2023, Joe Zbiciak <joe.zbiciak@leftturnonly.info>
`SPDX-License-Identifier:  CC-BY-SA-4.0`
