#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include "/opt/homebrew/include/gmp.h"

// int modinv(int a, int m)
// {
//     int m0 = m;
//     int y = 0, x = 1;

//     if (m == 1)
//         return 0;

//     while (a > 1)
//     {
//         int q = a / m;
//         int t = m;

//         m = a % m;
//         a = t;
//         t = y;

//         y = x - q * y;
//         x = t;
//     }

//     if (x < 0)
//         x += m0;

//     return x;
// }

// int main(){
//     printf("5 %% 7 : %d\n", modinv(5, 7));
//     printf("20 %% 31 : %d\n", modinv(20, 31));
//     printf("58 %% 137 : %d\n", modinv(58, 137));
//     return 0;
// }

void modinv(mpz_t result, const mpz_t a, const mpz_t m)
{
    mpz_t m0, mc, ac, y, x, q, t;
    mpz_inits(m0, mc, ac, y, x, q, t, NULL);
    mpz_set(m0, m);
    mpz_set(mc, m);
    mpz_set_ui(y, 0);
    mpz_set_ui(x, 1);

    if (mpz_cmp_ui(m, 1) == 0)
    {
        mpz_set_ui(result, 0);
        return;
    }

    while (mpz_cmp_ui(ac, 1) > 0)
    {
        mpz_fdiv_q(q, ac, mc);
        mpz_set(t, mc);

        mpz_mod(mc, ac, mc);
        mpz_set(ac, t);
        mpz_set(t, y);

        mpz_submul(y, q, y);
        mpz_set(x, t);
    }

    if (mpz_cmp_ui(x, 0) < 0)
        mpz_add(x, x, m0);

    mpz_set(result, x);
    mpz_clears(m0, mc, ac, y, x, q, t, NULL);
}


// int main(){
//     mpz_t r, five, seven, t20, t31, t58, t137;
//     mpz_inits(r, five, seven, t20, t31, t58, t137,NULL);
//     mpz_set_ui(five, 5);
//     mpz_set_ui(seven, 7);
//     mpz_set_ui(t20, 20);
//     mpz_set_ui(t31, 31);
//     mpz_set_ui(t58, 58);
//     mpz_set_ui(t137, 137);
//     printf("got here!");

//     modinv(r, five, seven);
//     printf("5 %% 7 : %s\n", mpz_get_str(NULL, 10, r));
//     modinv(r, t20, t31);
//     printf("20 %% 31 : %s\n", mpz_get_str(NULL, 10, r));
//     modinv(r, t58, t137);
//     printf("58 %% 137 : %s\n", mpz_get_str(NULL, 10, r));
//     mpz_clears(r, five, seven, t20, t31, t58, t137, NULL);
//     return 0;
// }

bool isprime_mr(const mpz_t n)
{
    int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                   59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int i = 0; i < 25; i++){
        if (mpz_fdiv_ui(n, primes[i])==0) return false;
    }
    printf("passed small primes\n");
    
    bool ret = true;
    mpz_t a, d, x, y, n1;
    mpz_inits(a, d, x, y, n1, NULL);

    size_t s = 0;
    mpz_sub_ui(d, n, 1);
    mpz_set(n1, d);
    while (mpz_even_p(d)) {
        s++;
        mpz_fdiv_q_ui(d, d, 2);
    }
    assert(mpz_odd_p(d) > 0);
    mpz_t test1, test2;
    mpz_inits(test1, test2);
    mpz_sub_ui(test1, n, 1);
    mpz_mul_2exp(test2, d, s);
    assert(mpz_cmp(test1, test2)==0);
    printf("d and s good\n");
    gmp_printf("%d %Zd\n", s, d);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (size_t)time(NULL));
    for (int i = 0; i < 25; i++){
        mpz_urandomm(a, state, n1);
        mpz_add(a,a,2);
        mpz_powm(x, a, d, n);
        printf("made random and %d\n", ret);

        if (mpz_cmp_ui(x, 1) == 0) continue;

        for (int j = 0; j < s; j++){
            if (mpz_cmp(x, n1) == 0) break;
            mpz_powm_ui(x, x, 2, n);
            gmp_printf("%Zd\n", x);
            // if (mpz_cmp_ui(y, 1) == 0 && mpz_cmp_ui(x, 1) != 0
            //     && mpz_cmp(x, n1) != 0) {
            //     ret = false;
            //     break;
            // }
            //mpz_set(x, y);
        }
        if (mpz_cmp(x, n1) == 0) continue;
        //if (mpz_cmp_ui(y, 1) != 0) ret = false;
        if (ret == false) break;
    }
    printf("about to clear\n");


    mpz_clears(a, d, x, y, n1, test1, test2, NULL);
    printf("cleared\n");
    return ret;
}
int main(){
    mpz_t r, five, seven, t20, t31, t58, t137;
    mpz_inits(r, five, seven, t20, t31, t58, t137,NULL);
    mpz_set_str(five, "997", 10);
    assert(isprime_mr(five));
    mpz_set_str(seven, "19b80661205fe32a935cb517e40c9d2b", 16);
    assert(isprime_mr(seven));
    mpz_set_str(t20, "867cf8f60dff3627afdbfcf862a84f11", 16);
    assert(isprime_mr(t20));
    mpz_set_str(t20, "867cf8f60dff3627afdbfcf862a84fff", 16);
    assert(!isprime_mr(t20));
    mpz_set_str(t31, "c755817b317fefe30f9dd2e7ed4190ff", 16);
    assert(isprime_mr(t31));
    mpz_set_str(t58, "b83a19ffd027e84e961e2758506db3e72d3cc9fadfc2a7feed4eee5f1f76f8693441180b1b95f17e8ab61dba6aed0e6da2f9989cd9740ccb37fa4c35f605c288", 16);
    assert(!isprime_mr(t58));
    mpz_set_str(t137, "fa08906ecee4aff268f1089fa88d56ca404c4d1438b1e0537539d2af93a3496b3b7839bf9b88810353bd1203a39e8a3d01e05c76e757e91227cf44241c4de269", 16);
    assert(isprime_mr(t137));
    mpz_clears(r, five, seven, t20, t31, t58, t137, NULL);
    return 0;
}
