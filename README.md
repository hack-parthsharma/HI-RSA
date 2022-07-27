HI-RSA
Highly Insecure RSA implementation

PURPOSE

This program was written as a proof of concept of RSA.
This implementation is NOT to be used for any real world encryption.
This encryption program demonstrates the basic mathematical concepts that
 make up the core of the RSA encryption system.

PROCESS

The two options at runtime allow for smaller primes to be used, where it
is easier to see what is happening, or larger primes, demonstrating the 
time requirement of the primality testing process. To test for primality, 
Miller's test is used with the first 20 small primes, with the
Miller-Rabin test completing the primality check. Finally, when using
option 2 (large primes), the ASCII message will only be stored in
capital letters. Encrypting an ASCII message with option 1 (small primes)
would require breaking the message into blocks due to the large message
size, which is not yet implemented.

Information about the GMP library is available at http://gmplib.org
The program has been tested and is known to work as of GMP 5.1.3.

COMPILATION

This program requires the GMP library and can be compiled as follows:
                g++ -Wall RSA.cpp -lgmpxx -lgmp -o RSA
Tip: if you are altering the code and get a segmentation fault, most
  likely a GMP variable is being used that hasn't been initialized!
