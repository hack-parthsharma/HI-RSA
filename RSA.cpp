/** 
 * HI-RSA
 * Highly Insecure RSA implementation
 * For a quick summary of RSA, wikipedia has a good basic explanation
 * https://en.wikipedia.org/wiki/RSA_%28algorithm%29
 *  
 * See the README file for more details
 **/

#include <stdio.h>
#include <ctype.h>
#include <gmp.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

//Global constants
  mpz_t ZERO;
  mpz_t ONE;
  mpz_t TWO;

void PrimeTest(mpz_t candidate) {
  // input an odd integer, output a number likely to be prime
  // Start by assuming that q is not prime, and only change this after last test
  bool isPrime = false;
  bool composite = false;
  // Setup all variables
  mpz_t subcandidate; // prime - 1
  mpz_init(subcandidate);
  mpz_t qodd; // odd cofactor of subcandidate
  mpz_init(qodd);
  mpz_t remainder;
  mpz_init(remainder);
  mpz_t divisor;
  mpz_init(divisor);
  // Need to try for first 20 small primes as bases in Miller's Test
  // Hand coded the smallest 20 primes for slight speed boost
  int smallprimes[20] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};

  if (mpz_even_p(candidate) != 0) {
      mpz_add(candidate, candidate, ONE);
  }

  while(!isPrime) {
    composite = false;
  // Ceiling for factoring is set arbitrarily. High is better
    for (int i = 3; i < 5000000; i+=2) { 
      mpz_set_d(divisor, i);
      mpz_mod(remainder, candidate, divisor);
      double rmd = mpz_get_d(remainder);
  // Optional debugging announce if candidate prime didn't work
      if (rmd == 0) {
  //      cout << "Remainder of prime divided by " << i << " = " << rmd << endl;
        i = 1;
        mpz_add(candidate, candidate, TWO);
      }
    }

  // Miller-Rabin primality test
  // Step 1: divide until an odd cofactor of (prime-1) is found
    int k = 0;
  // Last minute disaster check - don't want an even candidate
    if(mpz_even_p(candidate) != 0) {
cout<<"Candidate is even!!"<<endl;
    }
    mpz_sub(subcandidate, candidate, ONE);
    while (mpz_even_p(subcandidate) != 0) {
      mpz_cdiv_q(subcandidate, subcandidate, TWO);
      k++;
    }
    // Because subcandidate should always be even, k != 0
    mpz_set(qodd, subcandidate); 
    mpz_sub(subcandidate, candidate, ONE);

  // Step 2: Set i = 0, r = residue of b^q (mod n)
  
    mpz_t base, ivar, rvar, kvar; // base, counter, remainder, (2^k)*qodd
    mpz_init(base);
    mpz_init(ivar);
    mpz_init(rvar);
    mpz_init(kvar);
    mpz_set_d(kvar, k);
    bool MTest = true;

    for(int f = 0; f<10; f++) {
    mpz_set(ivar, ZERO);
    mpz_set_d(base, smallprimes[f]);
    while(MTest) {
    // Step 3: Check if prime (comparisons)
      mpz_powm(rvar, base, qodd, candidate); 
      if(((mpz_cmp(rvar, ONE) == 0) & (mpz_cmp(ivar, ZERO) == 0)) || ((mpz_cmp(ivar, ZERO) >= 0) & (mpz_cmp(rvar, subcandidate) == 0)) & !composite) {
        isPrime = true;
        MTest = false;
      }

    // Step 4: Setup for next iteration
      mpz_add(ivar, ivar, ONE);
      mpz_powm(rvar, rvar, TWO, candidate); // Test base^(2^i)*qodd
      if((mpz_cmp(ivar, kvar) >= 0) & !isPrime) {
        composite = true;
        MTest = false;
      }
    }
    }
 
    if(composite) {
      isPrime = false;
      mpz_add(candidate, candidate, TWO);
    }
  // End of Miller's test
}
}
// End of primality test


void SmallPrimes() {
  // Since this is just a quick demo that the RSA mathematics works
  //   the integers have been preselected. If you wish to change the
  //   given values, remember that e is coprime to tot(n), that p and
  //   q are prime, and that the plaintext should match the initial
  //   message if everything works properly.
  mpz_t q; // prime q
  mpz_t subq;
  mpz_t p; // prime p
  mpz_t subp;
  mpz_t e; // encryption exponent
  mpz_t d; // decryption exponent
  mpz_t n; // n = p*q
  mpz_t totn; // totient function of n = (p-1)(q-1)
  mpz_t message; // incoming message to encrypt
  mpz_t ct; // ciphertest
  mpz_t pt; // plaintext
  mpz_t gcdetotn;
  mpz_init(q);
  mpz_init(subq);
  mpz_init(p);
  mpz_init(subp);
  mpz_init(e);
  mpz_init(d);
  mpz_init(n);
  mpz_init(totn);
  mpz_init(message);
  mpz_init(ct);
  mpz_init(pt);
  mpz_init(gcdetotn);
  mpz_set_d(q, 47);
  mpz_sub_ui(subq, q, 1);
  mpz_set_d(p, 31);
  mpz_sub_ui(subp, p, 1);
  mpz_set_d(e, 7);
  mpz_set_d(ONE, 1);
  mpz_mul(totn, subq, subp);
  mpz_mul(n, p, q);

  // Need gcd(e, tot(n)) to be 1
  mpz_gcd(gcdetotn, e, totn);
  while(mpz_cmp(gcdetotn, ONE) > 0) {
  // Increment e to avoid gcd(e, tot(n)) > 1
    mpz_add_ui(e, e, 2);
    mpz_gcd(gcdetotn, e, totn);
  } 
  // Catch any last errors
  if(mpz_cmp(gcdetotn, ZERO) == 0) {
    cout << "ERROR! gcd(e, tot(n)) = 0 !" << endl;
  }
  gmp_printf("q = %Zd\n", q);
  gmp_printf("p = %Zd\n", p);
  gmp_printf("e = %Zd\n", e);
  gmp_printf("d = %Zd\n", d);
  gmp_printf("%Zd is the value of n\n", n);
  cout << "Enter an integer less than n to serve as a message: " << endl;
  int input = 0;
  cin >> input;
  mpz_set_d(message, input);
  if(mpz_cmp(message, n) > 0) {
    cout << "Error! Message is too big!\n";
  } else {
    mpz_invert(d, e, totn);
    mpz_powm(ct, message, e, n);
    mpz_powm(pt, ct, d, n);
    gmp_printf("Encrypted message: %Zd\n", ct);
    gmp_printf("Decrypted ciphertext: %Zd\n", pt);
    
  }
}



void LargePrimes() {

  int strength = 256;
  // Don't want the difference bwteen p and q to be too small
  // It avoids factoring with Fermat's factorization algorithm
  int pstrength = strength+12;
  int qstrength = strength-12;

  // p and q are eventual primes, n = p*q
  mpz_t p; // prime p used for encryption
  mpz_init(p);
  mpz_t q; // prime q used for encryption
  mpz_init(q);
  mpz_t n; // n = p*q
  mpz_init(n);
  mpz_t totn; // totient function of n
  mpz_init(totn);
  mpz_t e; // encryption exponent (currently hardwired in)
  mpz_init(e);
  mpz_set_d(e, 65537); // recommended value for e
  mpz_t d; // decryption exponent
  mpz_init(d);
  mpz_init(ZERO);
  mpz_set_d(ZERO, 0);
  mpz_init(ONE);
  mpz_set_d(ONE, 1);
  mpz_init(TWO);
  mpz_set_d(TWO, 2);
  mpz_t candremainder; // remainder of candidate during even/odd test
  mpz_init(candremainder);
  gmp_randstate_t rand; // holds random data used to get p and q
  gmp_randinit_mt(rand);
  gmp_randinit_lc_2exp_size(rand, strength);
 
// Now create candidate primes P, Q
// Get random candidate primes
// Use /dev/urandom for a pseudorandom seed

  FILE * fp = fopen("/dev/urandom", "r");
  int seed;
  fread(&seed, sizeof(int), 1, fp);
  gmp_randseed_ui(rand, seed);
  mpz_urandomb(q, rand, qstrength);
  gmp_printf("Possible q value\n  %Zd\n", q);

  FILE * ft = fopen("/dev/urandom", "r"); // need to check whether this alters the seed compared to fp
  fread(&seed, sizeof(int), 1, ft);
  gmp_randseed_ui(rand, seed);
  mpz_urandomb(p, rand, pstrength);
  gmp_printf("Possible p value\n  %Zd\n", p);

  // Test if these candidate P, Q values are likely prime
  // Start by testing q
  // First check if q is even. If so, make q odd

  mpz_mod(candremainder, q, TWO);
  double rem = mpz_get_d(candremainder);
    if (rem == 0) {
      mpz_add(q, q, ONE);
  }

  // Transform random p, q to prime p, q
  PrimeTest(q);
  gmp_printf("Post-prime test q\n  %Zd\n", q);
  PrimeTest(p);
  gmp_printf("Post-prime test p\n  %Zd\n", p);

 /** 
 *    Assume that we now know that p and q are large primes
 *    Also assume that RSA is as hard as factoring large number and that NP != P
 *    Then it is time to begin encrypting/decrypting
 **/

  // Find n, totn, and d
  mpz_t subp;
  mpz_t subq;
  mpz_init(subq);
  mpz_init(subp);
  mpz_sub(subp, p, ONE);
  mpz_sub(subq, q, ONE);
  mpz_mul(n, p, q);
  mpz_mul(totn, subp, subq);

  // Need to check gcd(e, totn) = 1
  mpz_t gcdvalue;
  mpz_init(gcdvalue);
  mpz_gcd(gcdvalue, e, totn);
  while( mpz_cmp(gcdvalue, ONE) > 0) {
    mpz_add_ui(e, e, 2);
    cout << "added 2 to e value";
  }

  // calculate decryption exponent
  mpz_invert(d, e, totn);
  gmp_printf("e = %Zd\n", e);
  gmp_printf("d = %Zd\n", d);
  gmp_printf("n = %Zd\n", n);
  if( mpz_cmp(d, n) > 0) {
    cout << "d is greater than n!!" << endl;
  }

  string message = "";
  cout << "Enter a short (alphabetic only) message to encrypt: " << endl;
  cin >> message;
  int mlength = message.length();
  // Convert message from chars into a big int
  mpz_t TEN;
  mpz_init(TEN);
  mpz_set_d(TEN, 10);
  mpz_t block;
  mpz_init(block);
  mpz_t powten;
  mpz_init(powten);
  mpz_t add2block;
  mpz_init(add2block);
  
  // iterate through each char of message
  for(int l = 0; l < mlength; l++) {
    char w = message.at(l);
    // Need uppercase message to allow for easier double digit ascii encoding
    int element = toupper(w);
    // Allocates two spaces for each char in string
    mpz_pow_ui(powten, TEN, l*2);
    mpz_set_d(add2block, element);
    mpz_mul(add2block, add2block, powten);
    //gmp_printf("add2block is %Zd\n", add2block);
    mpz_add(block, block, add2block);
  }

  // Encryption take in (n, e) as public key pair
  
  mpz_t ciphertext;
  mpz_init(ciphertext);
  // Actual encryption step, which requires that the block be less than n
  if(mpz_cmp(block, n) < 0) {
  // c = block^e (mod n)
    mpz_powm(ciphertext, block, e, n);
  } else {
    cout << "Block is too big!!";
  }
  cout<<endl; 
  gmp_printf("Numericized input message: %Zd\n", block);
  gmp_printf("Encrypted input message: %Zd\n", ciphertext);

  // Decryption of ciphertext takes in (n, d) as private key pair
  mpz_t plaintext;
  mpz_init(plaintext);
  //Actual decryption
  mpz_powm(plaintext, ciphertext, d, n);
  gmp_printf("Decrypted ciphertext: %Zd\n", plaintext);

  // Want to convert output int to string
  mpz_t ptlength;
  mpz_init(ptlength);
  int numlength =  mpz_sizeinbase(plaintext, 10);
  mpz_t letter;
  mpz_init(letter);
  mpz_t modpointer;
  mpz_init(modpointer);
  mpz_t divpointer;
  mpz_init(divpointer);
  mpz_t asciicode;
  mpz_init(asciicode);
  mpz_t pwrten, textexp; // Used during text decryption
  mpz_init(pwrten);
  mpz_init(textexp);

  cout << "Decrypted text: ";  
  for(int y = 0; y < numlength/2; y++) {
    //First calculate division offset
    mpz_set_d(textexp, 2*y);
    mpz_powm(pwrten, TEN, textexp, n); // pwrten = 10^textexp (mod n)
    mpz_set(divpointer, pwrten);
//gmp_printf("divpointer - %Zd\n", divpointer);
    // Now calculate modulus offset
    mpz_add(textexp, textexp, TWO);
    mpz_powm(pwrten, TEN, textexp, n);
    mpz_set(modpointer, pwrten);
//gmp_printf("modpointer - %Zd\n", modpointer);
 
    mpz_mod(letter, plaintext, modpointer);
    mpz_fdiv_q(letter, letter, divpointer);
    int asciicode = mpz_get_d(letter);
    char l = asciicode;
    cout << l; // decoded plaintext integer
    // Allocates two spaces for each char in string
  }
  cout << endl;

 // Check whether the primality tests are working well
  cout << "If 0, q is composite; if 1, q is likely prime - " << mpz_probab_prime_p(p, 100) << endl;
  cout << "If 0, p is composite; if 1, p is likely prime - " <<  mpz_probab_prime_p(p, 100) << endl;

}

int main() {
  
  cout << endl << "--  HI-RSA initiated  --" << endl << endl;
  cout << "Please enter an input of 1 if you want to test small primes or an input of 2 for a large prime test: " << endl;
  int choice = 0;
  cin >> choice;
  if (choice == 1) {
    SmallPrimes();
  } else if (choice == 2) {
    LargePrimes();
  } else {
    cout << "You entered an invalid input!" << endl;
    return 1;
  }
}

