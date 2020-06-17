# In this script, we implement the key schedule specified in the AES standard

from numpy.polynomial import Polynomial as P
import json as js


def reduce_galois(byte):
    # this function is used to change all the polynomial coefficients in the equivalence
    # class of 0 or 1 to 0 or 1
    # Input: byte Output: Reduced byte
    new_byte = [0 * (coefficient % 2 == 0) + 1 * (coefficient % 2 != 0) for coefficient in byte]
    return new_byte


def read_key(key_string):
    # This function reads a binary key in string format and it is out is a polynomial
    # split into words if size 32 bits
    # Input: Key in string format Output: Key in Polynomial format split into words
    key = []
    if len(key_string)==128:
        ind_array = [0, 4, 8, 12]
    elif len(key_string)==192:
        ind_array = [0, 4, 8, 12, 16, 20]
    elif len(key_string)==256:
        ind_array = [0, 4, 8, 12, 16, 20, 24, 28]
    for j in ind_array:
        key.append(P([int(bit) for i in range(j, j+4) for bit in key_string[i*8:(i+1)*8][::-1]]))
    return key


def h(word, s_box):
    # This function implements the h-box which is in the key schedule of the AES-256
    # Input: 32-bit Word in polynomial format, SBox
    # Output: Result of the h-box operation (32-bit Word in polynomial format)
    first_byte = ''.join([str(int(elem)) for elem in word.coef[0:8][::-1]])
    second_byte = ''.join([str(int(elem)) for elem in word.coef[8:16][::-1]])
    third_byte = ''.join([str(int(elem)) for elem in word.coef[16:24][::-1]])
    fourth_byte = ''.join([str(int(elem)) for elem in word.coef[24:32][::-1]])
    (first_byte, second_byte, third_byte, fourth_byte) = \
        (s_box[first_byte], s_box[second_byte], s_box[third_byte], s_box[fourth_byte])
    new_word = P([int(bit) for bit in first_byte[::-1]] + [int(bit) for bit in second_byte[::-1]] +
    [int(bit) for bit in third_byte[::-1]] + [int(bit) for bit in fourth_byte[::-1]])
    return new_word


def g(word, rc, s_box):
    # This function implements the g-box which is in the key schedules for all key sizes
    # Input: 32-bit Word in polynomial format, current round coefficient ,SBox
    # Output: Result of the h-box operation (32-bit Word in polynomial format)
    first_byte = ''.join([str(int(elem)) for elem in word.coef[0:8][::-1]])
    second_byte = ''.join([str(int(elem)) for elem in word.coef[8:16][::-1]])
    third_byte = ''.join([str(int(elem)) for elem in word.coef[16:24][::-1]])
    fourth_byte = ''.join([str(int(elem)) for elem in word.coef[24:32][::-1]])
    (first_byte, second_byte, third_byte, fourth_byte) = \
        (s_box[second_byte], s_box[third_byte], s_box[fourth_byte], s_box[first_byte])
    new_word = P([int(bit) for bit in first_byte[::-1]] + [int(bit) for bit in second_byte[::-1]]
                 + [int(bit) for bit in third_byte[::-1]] + [int(bit) for bit in fourth_byte[::-1]])
    temp = P(new_word.coef[0:8]) + rc
    temp = reduce_galois(temp.coef)
    for i in range(8 - len(temp)):
        temp.append(0)
    new_word.coef[0:8] = temp
    return new_word


def genKeyState(key_string, base=16):
    # This function is used to generate the keys for all rounds
    # Input: Key in string format, base of the key (optional argument and the default is base = 16)
    # Output: All the round keys in hex stored in a list
    hex_keys = []
    SBox = js.load(open("S-Box.json"))
    modulus = P([1, 1, 0, 1, 1, 0, 0, 0, 1])
    if base == 16:
        if len(key_string) == 128/4:
            binaryData = bin(int(key_string, 16))[2:].zfill(128)
            key_size = 128
        elif len(key_string) == 192/4:
            binaryData = bin(int(key_string, 16))[2:].zfill(192)
            key_size = 192
        elif len(key_string) == 256/4:
            binaryData = bin(int(key_string, 16))[2:].zfill(256)
            key_size = 256
        else:
            raise Exception('The Key size is invalid')
    elif base == 2:
        binaryData = key_string
    elif base != 2:
        raise Exception('The input data should be in either base 16 or 2')

    W = read_key(binaryData)
    rc = P([1])
    if key_size == 128:
        rounds = 10
        for i in range(rounds):
            temp_list = []
            for j in range(4):
                if j == 0:
                    word = reduce_galois((W[len(W) - 4] + g(W[len(W) - 1], rc, SBox)).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
                else:
                    word = reduce_galois((W[len(W) - (4 - j)] + prev_word).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
            W.extend(temp_list)
            _, rc = divmod(rc * P([0, 1]), modulus)
    elif key_size == 192:
        rounds = 8
        for i in range(rounds):
            temp_list = []
            for j in range(6):
                if j == 0:
                    word = reduce_galois((W[len(W) - 6] + g(W[len(W) - 1], rc, SBox)).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
                else:
                    word = reduce_galois((W[len(W) - (6 - j)] + prev_word).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
            W.extend(temp_list)
            _, rc = divmod(rc * P([0, 1]), modulus)
    elif key_size == 256:
        rounds = 7
        for i in range(rounds):
            temp_list = []
            for j in range(8):
                if j == 0:
                    word = reduce_galois((W[len(W) - 8] + g(W[len(W) - 1], rc, SBox)).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
                elif j == 4:
                    word = reduce_galois((W[len(W) - (8 - j)] + h(prev_word, SBox)).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
                else:
                    word = reduce_galois((W[len(W) - (8 - j)] + prev_word).coef)
                    while len(word) < 32:
                        word.append(0)
                    temp_list.append(P(word))
                    prev_word = P(word)
            W.extend(temp_list)
            _, rc = divmod(rc * P([0, 1]), modulus)
    #The next part converts the polynomials to hexadecimal values
    list_key = []
    counter = 0
    key_counter = 0
    for word in W:
        counter += 1
        list_key.extend(word.coef[0:8][::-1])
        list_key.extend(word.coef[8:16][::-1])
        list_key.extend(word.coef[16:24][::-1])
        list_key.extend(word.coef[24:32][::-1])
        if counter == 4:
            expanded_key_binary = ''.join([str(int(bit)) for bit in list_key])
            hexKey = hex(int(expanded_key_binary, 2))[2:]
            while 32 > len(hexKey):
                hexKey = "0" + hexKey
            hex_keys.append(hexKey)
            counter = 0
            key_counter += 1
            list_key = []
        if key_counter==15:
            break
    return hex_keys
