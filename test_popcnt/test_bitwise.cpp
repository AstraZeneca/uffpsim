#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <immintrin.h>

#define UINT64CONST(x) (x##ULL)

inline uint64_t bitwise_and_popcount_avx2(uint64_t* a, uint64_t* b, const size_t n) {
    size_t i = 0;
    uint64_t result = 0;
    __m256i accum = _mm256_setzero_si256();
    while (i + 4 <= n) {
        const __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        const __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        const __m256i vc = _mm256_and_si256(va, vb);
        
        // pocount
        result += __builtin_popcountll(_mm256_extract_epi64(vc, 0)); 
        result += __builtin_popcountll(_mm256_extract_epi64(vc, 1)); 
        result += __builtin_popcountll(_mm256_extract_epi64(vc, 2)); 
        result += __builtin_popcountll(_mm256_extract_epi64(vc, 3)); 

        i += 4;
    }
  return result;
}

inline uint64_t bitwise_avx512vpopcntdq(uint64_t* a, uint64_t* b, const size_t n) {
    size_t i = 0;
    __m512i cnt;
    __m512i accum = _mm512_setzero_si512();
    while (i + 8 <= n) {
        // bitwise AND
        const __m512i va = _mm512_loadu_si512((const __m512i *)(a + i));
        const __m512i vb = _mm512_loadu_si512((const __m512i *)(b + i));
        const __m512i vc = _mm512_and_si512(va, vb);

        // pocount
        cnt = _mm512_popcnt_epi64(vc);
        accum = _mm512_add_epi64(accum, cnt);

        i += 8;
    }
    return _mm512_reduce_add_epi64(accum);
}

uint64_t get_posix_clock_time ()
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (uint64_t) (ts.tv_sec * 1000000 + ts.tv_nsec );
    else
        return 0;
}

int main() {
    uint64_t prev_time_value, time_value;
    uint64_t time_diff;

    int nbits = 2048;
    size_t size = nbits/64;
    uint64_t* A = new uint64_t[size];
    uint64_t* B = new uint64_t[size];
    uint64_t* C = new uint64_t[size];
    uint64_t bits=0;

    memset(C, 0xff, size*sizeof(uint64_t));

    srand((unsigned) time(0));
    for (size_t i = 0; i < size; i++){
        uint64_t num = (uint64_t) rand();
        A[i] = num;
        printf("%u ", num);
    }
    printf("\n");

    for (size_t i = 0; i < size; i++){
        uint64_t num = (uint64_t) rand();
        B[i] = num;
        printf("%u ", num);
    }
    printf("\n");


    
    bits=0;
    prev_time_value = get_posix_clock_time();
    for(int i=0; i< size; i++)
        bits += __builtin_popcountll(A[i] & B[i]);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);

      
    bits=0;
    prev_time_value = get_posix_clock_time();
    for(int i=0; i< size; i++)
        bits += __builtin_popcountll(A[i] & B[i]);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);

    /*
    prev_time_value = get_posix_clock_time();
    andOp(A, B, C, size);
    bits = popcnt(C, sizeof(uint64_t)*size);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);
    */

    size_t size_8bit = nbits/8;
    prev_time_value = get_posix_clock_time();
    bits = bitwise_and_popcount_avx2(A, B, size);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);

    prev_time_value = get_posix_clock_time();
    bits = bitwise_and_popcount_avx2(A, B, size);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);
    
    prev_time_value = get_posix_clock_time();
    bits = bitwise_avx512vpopcntdq(A, B, size);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);

    prev_time_value = get_posix_clock_time();
    bits = bitwise_avx512vpopcntdq(A, B, size);
    time_value = get_posix_clock_time ();
    printf("%u %u \n", bits, time_value-prev_time_value);
}