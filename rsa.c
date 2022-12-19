// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>

// // Computes the modular inverse of a (modulo m) using the Extended Euclidean algorithm.
// // Returns -1 if no modular inverse exists.
// int modinv(int a, int m) {
// int m0 = m;
// int y = 0, x = 1;

// if (m == 1) return 0;

// while (a > 1) {
// int q = a / m;
// int t = m;

// m = a % m, a = t;
// t = y;

// y = x - q * y;
// x = t;


// }

// if (x < 0) x += m0;

// return x;
// }

// // Returns 1 if n is prime, and 0 otherwise.
// int isprime(int n) {
// if (n <= 1) return 0;
// if (n == 2) return 1;
// if (n % 2 == 0) return 0;

// for (int i = 3; i <= sqrt(n); i += 2) {
// if (n % i == 0) return 0;
// }

// return 1;
// }

// // Generates a random prime number with the given number of bits.
// int randprime(int bits) {
// // Generate a random number with the given number of bits.
// int p = rand() % (1 << (bits - 1)) + (1 << (bits - 1));

// // Keep incrementing p until it is prime.
// while (!isprime(p)) p++;

// return p;
// }

// // Encrypts the given message using the RSA algorithm with the given public key (n, e).
// int encrypt(int m, int n, int e) {
// return (int)pow(m, e) % n;
// }

// // Decrypts the given message using the RSA algorithm with the given private key (n, d).
// int decrypt(int c, int n, int d) {
// return (int)pow(c, d) % n;
// }

// int main() {
// // Generate the RSA key pair.
// int p = randprime(32);
// int q = randprime(32);
// int n = p * q;
// int phi = (p - 1) * (q - 1);
// int e = 3; // Choose a small odd integer for the public key.
// int d = modinv(e, phi); // Compute the modular inverse of e (modulo phi).

// printf("p: %d\n", p);
// printf("q: %d\n", q);
// printf("n: %d\n", n);
// printf("phi: %d\n", phi);
// printf("e: %d\n", e);
// printf("d: %d\n", d);

// // Encode the message using a suitable scheme (e.g. ASCII or Unicode).
// int *m = (int*)"Hello, world!";

// printf("Message: %d\n", *m);

// printf("Message enc and dec: %d\n", decrypt(encrypt(*m,n,e),n,d));

// }





// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>

// // Computes the greatest common divisor of two integers using the Euclidean algorithm
// int gcd(int a, int b)
// {
//     if (b == 0)
//         return a;
//     else
//         return gcd(b, a % b);
// }

// // Computes the modular inverse of a (modulo m) using the Extended Euclidean algorithm
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

// // Generates an RSA key pair with the specified bit length
// void generate_key_pair(int nbits, int *n, int *e, int *d)
// {
//     // Generate two prime numbers p and q
//     int p = rand() % nbits;
//     int q = rand() % nbits;

//     // Compute the modulus n = p * q
//     *n = p * q;

//     // Compute the totient of n, which is equal to (p-1) * (q-1)
//     int totient = (p - 1) * (q - 1);

//     // Choose a small odd integer e that is relatively prime to the totient
//     *e = rand() % totient;
//     while (gcd(*e, totient) != 1)
//         *e += 2;

//     // Use the Extended Euclidean algorithm to compute the modular inverse of e (modulo the totient)
//     *d = modinv(*e, totient);
// }

// // Encrypts the specified message using the RSA algorithm
// void encrypt(int n, int e, const char *msg, char *enc)
// {
//     int m = 0;
//     int len = strlen(msg);

//     // Convert the message to a numerical representation using ASCII
//     for (int i = 0; i < len; i++)
//         m = m * 256 + (int)msg[i];

//     // Compute the encrypted message C = M^e (mod n)
//     int c = (int)pow(m, e);
//     c = c % n;

//     // Convert the encrypted message back to a string representation
//     sprintf(enc, "%d", c);
// }

// Decrypts the specified message using the RSA algorithm
// void decrypt(int n, int d, const char *enc, char *dec)
// {
//     int c = 0;
//     int len = strlen(enc);

//     // Convert the encrypted message to a numerical representation
//     for (int i = 0; i < len; i++)
//         c = c * 10 + (enc[i] - '0');

//     // Compute the decrypted message M = C^d (mod n)
//     int m = (int)pow(c, d);
//     m = m % n;

//     // Convert the decrypted message back to a string representation using ASCII
//     for (int i = 0; i < len; i++)
//     {
//         dec[i] = m % 256;
//         m = m / 256;
//     }

//     // Add a null terminator to the decrypted message string
//     dec[len] = '\0';
// }


// int main()
// {
//     // Generate an RSA key pair with a bit length of 1024
//     int n, e, d;
//     generate_key_pair(1024, &n, &e, &d);

//     // Encrypt the message "Hello world" using the RSA algorithm
//     char enc[1024];
//     encrypt(n, e, "Hello world", enc);

//     // Decrypt the encrypted message using the RSA algorithm
//     char dec[1024];
//     decrypt(n, d, enc, dec);

//     // Print the original message, the encrypted message, and the decrypted message
//     printf("Original message: %s\n", "Hello world");
//     printf("Encrypted message: %s\n", enc);
//     printf("Decrypted message: %s\n", dec);

//     return 0;
// }





#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <assert.h>
//#include <gmp.h>
#include "/opt/homebrew/include/gmp.h"

typedef enum { VERBOSE_MODE, QUIET_MODE } verb_t;

// Computes the greatest common divisor of two integers using the Euclidean algorithm
// Alternatively use mpz_gcd()
void gcd(mpz_t result, const mpz_t a, const mpz_t b)
{
    if (mpz_cmp_ui(b, 0) == 0)
        mpz_set(result, a);
    else
    {
        mpz_t r;
        mpz_init(r);
        mpz_mod(r, a, b);
        gcd(result, b, r);
        mpz_clear(r);
    }
}

// Computes the modular inverse of a (modulo m) using the Extended Euclidean algorithm
// Alternatively use mpz_invert()
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

// Checks if the n is prime by dividing up to sqrt(n)
bool isprime_naive(const mpz_t n)
{
    if (mpz_cmp_ui(n, 1) <= 0) return false;
    if (mpz_cmp_ui(n, 2) == 0) return true;
    mpz_t r; mpz_init(r);
    mpz_mod_ui(r, n, 2);
    if (mpz_cmp_ui(r, 0) == 0) {mpz_clear(r); return false;}

    mpz_t check;
    mpz_init(check);
    mpz_sqrt(check, n);
    for (int i = 3; mpz_cmp_ui(check, i) > 0; i += 2) {
        mpz_mod_ui(r, n, (unsigned long)i);
        if (mpz_cmp_ui(r, 0) == 0) {mpz_clears(r, check, NULL); return false;}
    }

    mpz_clears(r, check, NULL);
    return true;
}

// Checks if the n is prime by dividing by primes + Miller-Rabin 25 rounds
bool isprime_mr(const mpz_t n)
{
    int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                   59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int i = 0; i < 25; i++){
        if (mpz_fdiv_ui(n, primes[i])==0) return false;
    }
    
    bool ret = true;
    mpz_t a, d, x, s, j, y;
    mpz_inits(a, d, x, s, y, NULL);

    mpz_sub_ui(d, n, 1);
    while (mpz_fdiv_ui(d, 2)==0) {
        mpz_add_ui(s, s, 1);
        mpz_fdiv_q_2exp(d, d, 1);
    }

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (size_t)time(NULL));
    for (int i = 0; i < 25; i++){
        mpz_urandomm(a, state, n);
        mpz_powm(x, a, d, n);

        for (mpz_init(j); mpz_cmp(j, s) < 0; mpz_add_ui(j, j, 1)){
            mpz_powm_ui(y, x, 2, n);
            if (mpz_cmp_ui(y, 1) == 0 && mpz_cmp_ui(x, 1) != 0) {
                ret = false;
                break;
            }
            mpz_set(x, y);
        }
        if (mpz_cmp_ui(y, 1) != 0) {ret = false; break;}
    }


    mpz_clears(a, d, x, j, y, NULL);
    return ret;
}

// Generates a random prime number with the given number of bits.
void randprime(mpz_t p, int bits, size_t seed) {
    // Generate a random number with the given number of bits.
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed+(size_t)time(NULL));
    mpz_urandomb(p, state, bits);

    if (mpz_fdiv_ui(p,2)==0) mpz_add_ui(p, p, 1);

    mpz_t t;
    mpz_init(t);
    mpz_tdiv_q_2exp(t, p, bits-1);
    if (mpz_cmp_ui(t, 1) < 0) randprime(p, bits, seed+1);
    
    // Keep incrementing p until it is prime
    // This is not super secure but probably good enough
    //while (!isprime_mr(p)) mpz_add_ui(p,p,2);
    mpz_nextprime(p, p);

    gmp_randclear(state);
    mpz_clear(t);
}

// Generates an RSA key pair with the specified bit length
void generate_key_pair(int nbits, mpz_t n, mpz_t e, mpz_t d, mpz_t p, mpz_t q)
{
    // Generate two prime numbers p and q
    randprime(p, nbits, 0xbc590c83fd);
    randprime(q, nbits, 0xa0530ab1f1);
    //assert(isprime(p) && isprime(q));

    // Compute the modulus n = p * q
    mpz_mul(n, p, q);

    // Compute the totient of n, which is equal to (p-1) * (q-1)
    mpz_t totient;
    mpz_init(totient);
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_mul(totient, p, q);

    // Reset p and q
    mpz_add_ui(p, p, 1);
    mpz_add_ui(q, q, 1);

    mpz_set_ui(e, 3);
    mpz_t res;
    mpz_init(res);
    //gcd(res, e, totient);
    mpz_gcd(res, e, totient);
    while (mpz_cmp_ui(res, 1) != 0){
        mpz_add_ui(e,e,2);
        //gcd(res, e, totient);
        mpz_gcd(res, e, totient);
    }

    //modinv(d, e, totient);
    mpz_invert(d, e, totient);

    mpz_clears(res, totient, NULL);
}

void encryptm(const mpz_t n, const mpz_t e, const char *msg, mpz_t enc)
{
    mpz_t m;
    mpz_inits(m, enc, NULL);
    mpz_set_ui(m, 0);
    size_t len = strlen(msg);

    // Convert the message to a numerical representation using ASCII
    // for (size_t i = 0; i < len; i++){
    //     mpz_mul_ui(m,m,256);
    //     mpz_add_ui(m,m, (uint8_t)msg[i]);
    // }
    mpz_import(m, len, 1, 1, 0, 0, msg);

    // Compute the encrypted message C = M^e (mod n)
    mpz_powm(enc, m, e, n);

    mpz_clear(m);
}

// Decrypts the specified message using the RSA algorithm
void decrypt(mpz_t n, mpz_t d, const mpz_t enc, char *dec)
{
    // Compute the decrypted message M = C^d (mod n)
    mpz_t m;
    mpz_init(m);
    mpz_powm(m, enc, d, n);

    // Convert the decrypted message back to a string representation using ASCII
    // size_t len = (mpz_sizeinbase(m, 2)+7)/8;
    // for (int i = len-1; i>=0; i--){
    //     dec[i] = (uint8_t)mpz_tdiv_q_ui(m, m, 256);
    // }

    size_t len;
    //*dec = calloc(1 + (mpz_sizeinbase(m, 2)+7)/8, sizeof(char));
    mpz_export(dec, &len, 1, 1, 0, 0, m);

    dec[len] = '\0';

    mpz_clear(m);
}

int generate_mode(char *msg, size_t p_bits, verb_t verb)
{
    // Generate an RSA key pair with a bit length of 1024
    mpz_t p, q, n, e, d, enc;
    mpz_inits(p, q, n, e, d, NULL);
    generate_key_pair(p_bits, n, e, d, p, q);

    // Encrypt message using the RSA algorithm
    encryptm(n, e, msg, enc);

    // Decrypt the encrypted message using the RSA algorithm
    char *dec=calloc(p_bits/4, sizeof(char));
    decrypt(n, d, enc, dec);

    char *enc_str = mpz_get_str(NULL, 16, enc);

    char *p_str = mpz_get_str(NULL, 16, p);
    char *q_str = mpz_get_str(NULL, 16, q);
    char *n_str = mpz_get_str(NULL, 16, n);
    char *e_str = mpz_get_str(NULL, 16, e);
    char *d_str = mpz_get_str(NULL, 16, d);

    // Print primes and keys
    if (verb == VERBOSE_MODE) printf("p: %s\n", p_str);
    if (verb == VERBOSE_MODE) printf("q: %s\n", q_str);
    printf("n: %s\n", n_str);
    printf("e: %s\n", e_str);
    printf("d: %s\n", d_str);

    // Print the original message, the encrypted message, and the decrypted message
    if (strlen(msg) > 0){
    if (verb == VERBOSE_MODE) printf("Original message: %s\n", msg);
    printf("Encrypted message: %s\n", enc_str);
    if (verb == VERBOSE_MODE) printf("Decrypted message: %s\n", dec);
    }

    // Safely erase private key values
    mpz_set_ui(p, 0);
    mpz_set_ui(q, 0);
    mpz_set_ui(d, 0);
    memset(d_str, 0, strlen(d_str));
    memset(p_str, 0, strlen(p_str));
    memset(q_str, 0, strlen(q_str));

    // Clean up
    mpz_clears(n, e, d, enc, p, q, NULL);
    free(msg);
    free(enc_str);
    free(p_str);
    free(q_str);
    free(n_str);
    free(e_str);
    free(d_str);
    free(dec);
    return 0;
}

int encrypt_mode(char *msg, char *n_tmp, char *e_tmp, verb_t verb)
{
    mpz_t n, e, enc;
    mpz_inits(n, e, NULL);
    
    mpz_set_str(n, n_tmp, 16);
    mpz_set_str(e, e_tmp, 16);

    // Encrypt message using the RSA algorithm
    encryptm(n, e, msg, enc);

    char *enc_str = mpz_get_str(NULL, 16, enc);
    printf("%s\n", enc_str);

    // Clean up
    mpz_clears(n, e, enc, NULL);
    free(msg);
    free(enc_str);
    return 0;
}

int decrypt_mode(char *msg, char *n_tmp, char *d_tmp, verb_t verb)
{
    mpz_t n, d, enc;
    mpz_inits(n, d, enc, NULL);
    
    mpz_set_str(n, n_tmp, 16);
    mpz_set_str(d, d_tmp, 16);
    mpz_set_str(enc, msg, 16);

    // Decrypt the encrypted message using the RSA algorithm
    char *dec = calloc(1+(mpz_sizeinbase(n, 2)+7)/8, sizeof(char));
    decrypt(n, d, enc, dec);

    printf("%s\n", dec);

    // Clean up
    mpz_clears(n, d, enc, NULL);
    free(msg);
    free(dec);
    return 0;
}


int main(int argc, char *argv[])
{
    int opt;
    verb_t verb = VERBOSE_MODE;
    enum { GENERATE_MODE, ENCRYPT_MODE, DECRYPT_MODE } mode = GENERATE_MODE;

    size_t p_bits = 512;

    char *e_tmp;
    char *n_tmp;
    char *d_tmp;
    
    while ((opt = getopt(argc, argv, "m:s:qe:n:d:h")) != -1) {
        switch (opt) {
        case 'm': //char *p = optarg;
            //for ( ; *p; ++p) *p = tolower(*p);
            if (strcmp(optarg, "decrypt") <= 0) mode = DECRYPT_MODE;
            else if (strcmp(optarg, "encrypt") <= 0) mode = ENCRYPT_MODE;
            break;
        case 'e': e_tmp = optarg; break;
        case 'n': n_tmp = optarg; break;
        case 'd': d_tmp = optarg; break;
        case 's': p_bits = (size_t)atoi(optarg); break;
        case 'q': verb = QUIET_MODE; break;
        case 'h':
            fprintf(stderr, "Usage: %s -m mode [-s level] [-q] secret phrase...\n", argv[0]);
            exit(1);
        }
    }

    char *msg = calloc(p_bits/4, sizeof(char));

    // Construct input string
    int i = optind;
    while (i < argc){
        //assert(strlen(msg)+strlen(argv[i])+1 < p_bits/4);
        if (!(strlen(msg)+strlen(argv[i])+1 < p_bits/4)){
            fprintf(stderr, "Secret phrase should be less than %lu chars\n", p_bits/4);
            exit(1);
        }
        if (i > optind) strcat(msg, " ");
        strcat(msg, argv[i]);
        i++;
    }
    

    switch (mode) {
    case GENERATE_MODE:
        //assert(p_bits % 64 == 0 && p_bits > 0);
        if (!(p_bits % 64 == 0 && p_bits > 0)){
            fprintf(stderr, "Security level (%lu) should be a positive multiple of 64\n", p_bits);
            exit(1);
        }

        if (verb == VERBOSE_MODE) printf("Encrypting at %lu bit security\n", p_bits);
        return generate_mode(msg, p_bits, verb);

    case ENCRYPT_MODE: return encrypt_mode(msg, n_tmp, e_tmp, verb);
    case DECRYPT_MODE: return decrypt_mode(msg, n_tmp, d_tmp, verb);
    }
}
