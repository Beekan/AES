# Mathematical Description of the S-Box
# It consists of two blocks:
# The 1st one is a galois field inversion block
# The 2nd one is an affine mapping block
# We will implement the two blocks and construct the S-Box and save it for later use
# Note: It shouldn't be generated at each AES Encryption or Decryption round

# We need to implement the EEA Algorithm for inversion
# imports
import numpy as np
import json as js
from numpy.polynomial import Polynomial as P


# Functions
def eea(poly, mod):
    # this function implements the extended euclidean algorithm
    # the input poly is the polynomial we are required to find its inverse
    # The input mod is the modulus
    # The output of this function is the inverse of poly where the modulus is mod
    t = 0
    new_t = 1
    r = mod
    new_r = poly
    i = 1
    while sum(new_r.coef) != 0:
        q, _ = divmod(r, new_r)
        (r, new_r) = (new_r, r - q * new_r)
        (t, new_t) = (new_t, t - q * new_t)
        new_r.coef = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in new_r.coef]
        r.coef = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in r.coef]
        new_t.coef = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in new_t.coef]
        if i > 1:
            t.coef = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in t.coef]
        i = i + 1
    inverse, _ = divmod(t, r)
    inverse.coef = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in inverse.coef]
    return inverse


def aff_transform(byte):
    # This function implements the AES affine transformation which is used after the modular
    # inverse operation to construct the S-Box
    A = np.array([[1, 0, 0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, 1, 1, 1], [1, 1, 1, 0, 0, 0, 1, 1],
                  [1, 1, 1, 1, 0, 0, 0, 1], [1, 1, 1, 1, 1, 0, 0, 0], [0, 1, 1, 1, 1, 1, 0, 0],
                  [0, 0, 1, 1, 1, 1, 1, 0], [0, 0, 0, 1, 1, 1, 1, 1]])
    b = np.array([1, 1, 0, 0, 0, 1, 1, 0])
    result = np.matmul(A, byte) + b
    result = [0 * (bit % 2 == 0) + 1 * (bit % 2 != 0)for bit in result]
    return result


def inv_aff_transform(byte):
    A = np.array([[0, 0, 1, 0, 0, 1, 0, 1], [1, 0, 0, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0, 0, 1],
                  [1, 0, 1, 0, 0, 1, 0, 0], [0, 1, 0, 1, 0, 0, 1, 0], [0, 0, 1, 0, 1, 0, 0, 1],
                  [1, 0, 0, 1, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0, 1, 0]])
    b = np.array([1, 0, 1, 0, 0, 0, 0, 0])
    result = np.matmul(A, byte) + b
    result = [0 * (bit % 2 == 0) + 1 * (bit % 2 != 0) for bit in result]
    return result


# Generate the S-Box
polynomials = []
modulus = P([1, 1, 0, 1, 1, 0, 0, 0, 1])
for i in range(2**8):
    polynomials.append(np.binary_repr(i, width=8))
SBox = {polynomial: None for polynomial in polynomials}
for polynomial in polynomials:
    poly = P([int(coefficient) for coefficient in polynomial])
    inverse = eea(poly, modulus)
    byte = [coefficient for coefficient in inverse]
    for i in range(8-len(byte)):
        byte.append(0)
    byte = aff_transform(byte)
    SBox[polynomial[::-1]] = ''.join([str(elem) for elem in byte[::-1]])
with open("S-Box.json", "w") as write_file:
    js.dump(SBox, write_file, sort_keys=True,indent=4, separators=(',', ': '))

# Generate the Inverse S-Box
polynomials = []
modulus = P([1, 1, 0, 1, 1, 0, 0, 0, 1])
for i in range(2**8):
    polynomials.append(np.binary_repr(i, width=8))
InverseSBox = {polynomial: None for polynomial in polynomials}
for polynomial in polynomials:
    byte = [int(coefficient) for coefficient in polynomial]
    byte = inv_aff_transform(byte)
    poly = P(byte)
    inverse = eea(poly, modulus)
    byte = [coefficient for coefficient in inverse]
    for i in range(8-len(byte)):
        byte.append(0)
    SBox[polynomial[::-1]] = ''.join([str(elem) for elem in byte[::-1]])
with open("InverseS-Box.json", "w") as write_file:
    js.dump(SBox, write_file, sort_keys=True,indent=4, separators=(',', ': '))
