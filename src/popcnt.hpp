/**
 * @brief Contains population count functions
 * 
 */

#pragma once

#include <cstdint>

#if defined(_MSC_VER) // windows intel
#include <nmmintrin.h>
#define popcntll(X) _mm_popcnt_u64(X)
#else // unix (linux, osx) intel / arm
#define popcntll(X) __builtin_popcountll(X)
#endif

#if defined(HAVE_AVX512VPOPCNTDQ)
#include <immintrin.h>

/**
 * The bitwise_and_popcount function calculates the bitwise AND of corresponding elements in two arrays,
 * and then counts the number of set bits (1s) in the result.
 *
 * The function uses SIMD instructions (AVX512VPOPCOUNTDQ) for improved performance on 
 * systems with AVX512VPOPCOUNTDQ support.
 * If the system does not support , it falls back to a scalar implementation using the popcntll intrinsic.
 *
 * @param a: A pointer to the first array of 64-bit unsigned integers.
 * @param b: A pointer to the second array of 64-bit unsigned integers.
 * @param n: The number of elements in the arrays.
 *
 * @return: The total count of set bits (1s) in the result of the bitwise AND operation.
 *
 * Note: The function assumes that the input arrays are aligned and have a size that is a multiple of 8.
 * 
 * @author Rajendra Kumar
 */
inline uint64_t bitwise_and_popcount(uint64_t* a, uint64_t* b, const size_t n) {
    __m512i cnt;
    __m512i accum = _mm512_setzero_si512();
    for (size_t i = 0; i < n ; i+=8) {
        // bitwise AND
        const __m512i va = _mm512_loadu_si512((const __m512i *)(a + i));
        const __m512i vb = _mm512_loadu_si512((const __m512i *)(b + i));
        const __m512i vc = _mm512_and_si512(va, vb);

        // pocount
        cnt = _mm512_popcnt_epi64(vc);
        accum = _mm512_add_epi64(accum, cnt);
    }
    return _mm512_reduce_add_epi64(accum);
}
#else

/**
 * The bitwise_and_popcount function calculates the bitwise AND of corresponding elements in two arrays,
 * and then counts the number of set bits (1s) in the result.
 *
 * If the system does not support AVX512VPOPCOUNTDQ, it falls back to this function.
 *
 * @param a: A pointer to the first array of 64-bit unsigned integers.
 * @param b: A pointer to the second array of 64-bit unsigned integers.
 * @param n: The number of elements in the arrays.
 *
 * @return: The total count of set bits (1s) in the result of the bitwise AND operation.
 *
 * Note: The function assumes that the input arrays are aligned and have a size that is a multiple of 8.
 */
inline uint64_t bitwise_and_popcount(uint64_t* a, uint64_t* b, const size_t n) {
    uint64_t common_popcnt = 0;
    for (size_t i=0; i < n ; i++) {
        common_popcnt += popcntll(a[i] & b[i]);
    }
    return common_popcnt;
}
#endif