import numpy as np
from numpy.polynomial import Polynomial as P
import json as js

def reduce_galois(byte):
    byte = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in byte]
    return byte

bytes = []
modulus = P([1, 1, 0, 1, 1, 0, 0, 0, 1])
multipliers = ['00000001', '00000010', '00000011','00001001', '00001011', '00001101', '00001110']
for i in range(2**8):
    bytes.append(np.binary_repr(i, width=8))
MulTable = {byte: {multiplier: None for multiplier in multipliers} for byte in bytes}
for byte in bytes:
    for multiplier in multipliers:
        polynomial = P([int(bit) for bit in byte][::-1])
        mulPoly = P([int(bit) for bit in multiplier][::-1])
        MulTable[byte][multiplier] = int(''.join(list(map(str, reduce_galois(map(int, divmod(polynomial*mulPoly, modulus)[1]
                                                                             .coef[::-1]))))).zfill(8),2)
with open("MulTable.json", "w") as write_file:
    js.dump(MulTable, write_file, sort_keys=True,indent=4, separators=(',', ': '))
print(MulTable['00000001'])
