# In this script, we implement the key schedule specified in the AES standard
# First, we implement the g-block

from numpy.polynomial import Polynomial as P
import json as js
import KeySchedule as KS
import numpy as np
import matplotlib.pyplot as plt



def pad(plaintext):
    padding_len=bin((128 - (len(plaintext) % 128))%128)[2:].zfill(8)
    plaintext=plaintext.ljust(len(plaintext)+int(padding_len,2), '0')
    return plaintext,padding_len


def unpad(plaintext,padding_len):
    if int(padding_len,2)>0:
        plaintext=plaintext[:-int(padding_len,2)]
    return plaintext



def sBox(inputData,SBox,base=2):
    """
    This function operates on multiples of 32 Hex Characters or 128 bits. Each 8 bits are swapped with their
    equivalent 8 bits in the dictonary passed to the function.
    
    A cautionary notice when using this funciton is that the data that does not form a complete state is 
    discarded.
    
    
    Inputs:
    inputData: A string containing the hex or binary
    base: The numerical base of the inputData whether it is in base 16 or 2
    
    Outputs:
    subHex: A string of hex charachters after byte substitution
    """
   
    #convert from hexadecimal to binary
    if base==16:
        inputData = bin(int(inputData, 16))[2:].zfill(len(inputData)*4)
    elif base !=2:
        raise Exception('The input data should be in either base 16 or 2')
    
    
    subBin=""
    for i in range(int(len(inputData)/128)):
        #loop over 16 elements (the number of elements in each state) where each element is 8 bits, 
        #then substitute those 8 bits using SBox
        subBin+=bin(int("".join([SBox[inputData[byte*8:(byte+1)*8]] for byte in range(i*16,(i+1)*16)]),2))[2:].zfill(128)
         
    return subBin

def createStates(inputData,base=2):
    """
    This function operates on multiples of 32 Hex Characters or 128 bits. Each byte is flipped and turned into a polynomial 
    as to facilitate the operations done on the polynomials such as EEA. An array of states containing polynomial elements is returned.
    
    A cautionary notice when using this funciton is that the data that does not form a complete state is 
    discarded.
    
    Inputs: 
    inputData: A string containing the hex or binary
    base: The numerical base of the inputData whether it is in base 16 or 2
    
    Outputs:
    states: An numpy array of states of shape (number of states, 4,4)
    
    
    """
    
    if base==16:
        inputData = bin(int(inputData, 16))[2:].zfill(len(inputData)*4)
    elif base !=2:
        raise Exception('The input data should be in either base 16 or 2')
    numStates=int(len(inputData)/128)
    states=[]
    for j in range(numStates):
        states.append(np.array([[inputData[byte*8:(byte+1)*8]] for byte in range(j*16,(j+1)*16)]).reshape((4,4)).T)
    return np.array(states)


def StateToHex(state):
    """
    This function operates on a single full state and return the hex value of that state in a string
    
    Inputs: 
    state: A polynmial np array of shape (4,4) with each element being a polynomial of 8 bits
    
    Outputs:
    hexValue: A string containing the hex value of that state.
    
    """
    
    return hex(int("".join(state.T.reshape(-1)),2))[2:].zfill(32)


def StateToBin(state):
    return "".join(state.T.reshape(-1))

def Addition(inputData, roundKey, inputDataBase=2,roundKeyBase=16):
    """
    This function takes on 2 Strings of either base 16 or 2 and adds them together using the xor operator.
    If the input data containes multiple of a single state length, then each state is added on its own to the roundKey.
    The result of these seperate additions is then concatintated together to form the result of the addition.
    
    A cautionary notice when using this funciton is that the data that does not form a complete state is 
    discarded.
    
    Inputs:
    inputData: A string of hex or binary charachters
    roundKey: A string of hex or binary characters
    inputDataBase: The base at which inputdata operates.
    roundKeyBase: The base at which roundKey operates
    
    Outputs:
    result: the result of the addition of the 2 strings has the same length as inputData
    
    """
    if(inputDataBase==2):
        stateLength=128
    else:
        stateLength=32

    result=""
    for i in range(int(len(inputData)/stateLength)):
        result+=bin(int(roundKey,roundKeyBase)^int(inputData[i*stateLength:(i+1)*stateLength],inputDataBase))[2:].zfill(128)
    if(result==""):
        raise Exception('The given data does not complete a state')
    return result

def mixColumns(state,mult):
    newstate=np.array(state)
    for i in range(4):
        newstate[0][i]=bin(mult[state[0][i]]["00000010"]^mult[state[1][i]]["00000011"]^mult[state[2][i]]["00000001"]^mult[state[3][i]]["00000001"])[2:].zfill(8)
        newstate[1][i]=bin(mult[state[0][i]]["00000001"]^mult[state[1][i]]["00000010"]^mult[state[2][i]]["00000011"]^mult[state[3][i]]["00000001"])[2:].zfill(8)
        newstate[2][i]=bin(mult[state[0][i]]["00000001"]^mult[state[1][i]]["00000001"]^mult[state[2][i]]["00000010"]^mult[state[3][i]]["00000011"])[2:].zfill(8)
        newstate[3][i]=bin(mult[state[0][i]]["00000011"]^mult[state[1][i]]["00000001"]^mult[state[2][i]]["00000001"]^mult[state[3][i]]["00000010"])[2:].zfill(8)
    return newstate


def InvmixColumns(state,mult):
    newstate=np.array(state)
    for i in range(4):
        newstate[0][i]=bin(mult[state[0][i]]["00001110"]^mult[state[1][i]]["00001011"]^mult[state[2][i]]["00001101"]^mult[state[3][i]]["00001001"])[2:].zfill(8)
        newstate[1][i]=bin(mult[state[0][i]]["00001001"]^mult[state[1][i]]["00001110"]^mult[state[2][i]]["00001011"]^mult[state[3][i]]["00001101"])[2:].zfill(8)
        newstate[2][i]=bin(mult[state[0][i]]["00001101"]^mult[state[1][i]]["00001001"]^mult[state[2][i]]["00001110"]^mult[state[3][i]]["00001011"])[2:].zfill(8)
        newstate[3][i]=bin(mult[state[0][i]]["00001011"]^mult[state[1][i]]["00001101"]^mult[state[2][i]]["00001001"]^mult[state[3][i]]["00001110"])[2:].zfill(8)
    return newstate

def shift_rows(state):
    # the function takes the state as an input and undergoes the shifts rows operation specified in the AES cipher
    (state[1][0], state[1][1], state[1][2], state[1][3]) = (state[1][1], state[1][2], state[1][3], state[1][0])
    (state[2][0], state[2][1], state[2][2], state[2][3]) = (state[2][2], state[2][3], state[2][0], state[2][1])
    (state[3][0], state[3][1], state[3][2], state[3][3]) = (state[3][3], state[3][0], state[3][1], state[3][2])

    return state


def inv_shift_rows(state):
    (state[1][0], state[1][1], state[1][2], state[1][3]) = (state[1][3], state[1][0], state[1][1], state[1][2])
    (state[2][0], state[2][1], state[2][2], state[2][3]) = (state[2][2], state[2][3], state[2][0], state[2][1])
    (state[3][0], state[3][1], state[3][2], state[3][3]) = (state[3][1], state[3][2], state[3][3], state[3][0])
    return state

def encrypt_state(inputData,roundKeys,mult,SBox,showRounds=True,rounds=10,inputDataBase=16):
    """
    inputs:
    inputData:  Either a hex string of length 32 or binary string of length 128
    roundkeys:  An array containg all the subkeys needed.(this is already computed in a higher level function to decrease the number of computations)
    mult:       A dictionary with format {element{element1:element2}} which is used to the remove the use of multiplications for increased performance
    SBox:       A dictionary with format {element1:element2} used for byte substitution
    showRounds: Used to show the output of each layer
    rounds:     the number of rounds needed for decryption according to the size of the initial key
    inputDataBase: either 2 or 16
    
    outputs:
    keyAdded:   A string containg binary characters of length 128
    """
    
    
    #convert to hex
    if (inputDataBase==16):
        inputData = bin(int(inputData, 16))[2:].zfill(int(len(inputData)*4))
    elif (inputDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')
    
    
        
    #first key Addition
    keyAdded=Addition(inputData,roundKeys[0])
    
    #AES encryption rounds
    for i in range(rounds):
        if(showRounds):
            print(f"Round {i+1}")
        subBin=sBox(keyAdded,SBox)
        states=createStates(subBin)
        if(showRounds):
            print(f"SBox output:         {StateToHex(states[0])}")
            
        states[0]=shift_rows(states[0])
        if(showRounds):
            print(f"shift rows output:   {StateToHex(states[0])}")
            
        if(i<rounds-1):
            states[0]=mixColumns(states[0],mult)
            if(showRounds):
                print(f"mix columns output:  {StateToHex(states[0])}")
                
        keyAdded=Addition(StateToBin(states[0]),roundKeys[i+1])
        if(showRounds):
            print(f"Key Addition output: {hex(int(keyAdded,2))[2:].zfill(32)}")
            print("")
            
    return keyAdded


def decrypt_state(cipherData,roundKeys,mult,ISBox,showRounds=True,rounds=10,cipherDataBase=16):
    """
    inputs:
    cipherData: Either a hex string of length 32 or binary string of length 128
    roundkeys:  An array containg all the subkeys needed.(this is already computed in a higher level function to decrease the number of computations)
    mult:       A dictionary with format {element{element1:element2}} which is used to the remove the use of multiplications for increased performance
    ISBox:      A dictionary with format {element1:element2} used for byte substitution
    showRounds: Used to show the output of each layer
    rounds:     the number of rounds needed for decryption according to the size of the initial key
    cipherDataBase: either 2 or 16
    
    outputs:
    keyAdded:   A string containg binary characters of length 128
    """
    
    if (cipherDataBase==16):
        cipherData = bin(int(cipherData, 2))[2:].zfill(int(len(cipherData)*4))
    elif (cipherDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')
    
    
    for i in range(rounds):
        if(showRounds):
            print(f"Round {i+1}")
        cipherData=Addition(cipherData,roundKeys[-(i+1)])
        if (showRounds):
            print(f"Key Addition output: {hex(int(cipherData,2))[2:].zfill(32)}")
        states=createStates(cipherData)
        
        if(i>0):
            states[0]=InvmixColumns(states[0],mult)
            if (showRounds):
                print(f"mix columns output:  {StateToHex(states[0])}")

        states[0]=inv_shift_rows(states[0])
        if (showRounds):
            print(f"shift rows output:   {StateToHex(states[0])}")
            
        cipherData=sBox(StateToBin(states[0]),ISBox)
        if (showRounds):
            print(f"SBox output:         {hex(int(cipherData,2))[2:].zfill(32)}")
            print("")
            
    keyAdded=Addition(cipherData,roundKeys[0])
    return keyAdded




def decrypt_cbc(inputData,initialKey,IV,showRounds=False,inputDataBase=16,initialKeyBase=16):
    """
    This function deals with either hex or binary strings resulting from encrypt_ecb.
    
    inputs:
    
    inputData:  Hex or binary string arrays
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    IV:         A hex string of length 32 charcters
    showRounds: show the output of each layer
    inputDataBase:  either 2 or 16
    initialKeyBase: either 2 or 16
    
    outputs:
    plaintext: this represents the decrypted data.
    """
    
    
    if (inputDataBase==16):
        inputData = bin(int(inputData, 16))[2:].zfill(int(len(inputData)*4))
    elif (inputDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')

    if (initialKeyBase==2):
        initialKey = hex(int(initialKey, 2))[2:].zfill(int(len(initialKey)/4))
    elif (initialKeyBase!=16):
        raise Exception('The input data should be in either base 16 or 2')
        
    assert (len(IV)==32)
    prev=bin(int(IV,16))[2:].zfill(128)
    
    #key Generation for all rounds
    roundKeys = KS.genKeyState(initialKey, initialKeyBase)
    rounds=len(roundKeys)-1
        
    #implementing dictionaries instead of calculation and opening them in a high level function to reduce compuation time
    mult = js.load(open("MulTable.json"))
    ISBox = js.load(open("InverseS-Box.json"))
    
    #cutting the input into elements with each of 128 bit representing a state after removing the last 8 bits which
    #represent the number of bits padded
    paddedLen=inputData[-8:]
    inputData=inputData[:-8]
    numStates=int(len(inputData)/128)
    states=[inputData[i*128:(i+1)*128] for i in range(numStates)]
    
    
    plaintext=""
    counter=1
    
    for i in range(len(states)):
        if(showRounds):
            print(f"State {counter}")
            
        plaintext+=Addition(decrypt_state(states[i],roundKeys,mult,ISBox,cipherDataBase=2,showRounds=showRounds,rounds=rounds)
                                              ,prev,inputDataBase=2,roundKeyBase=2)
        prev=states[i]
        if(showRounds):
            counter+=1
            print("")
    plaintext=unpad(plaintext,paddedLen)
 
    return hex(int(plaintext,2))[2:].zfill(int(len(plaintext)/4))


def decrypt_ecb(inputData,initialKey,showRounds=False,inputDataBase=16,initialKeyBase=16):
    """
    This function deals with either hex or binary strings resulting from encrypt_ecb.
    
    inputs:
    
    inputData:  Hex or binary string arrays
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    showRounds: show the output of each layer
    inputDataBase:  either 2 or 16
    initialKeyBase: either 2 or 16
    
    outputs:
    plaintext: this represents the decrypted data.
    """
    
    
    if (inputDataBase==16):
        inputData = bin(int(inputData, 16))[2:].zfill(int(len(inputData)*4))
    elif (inputDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')

    if (initialKeyBase==2):
        initialKey = hex(int(initialKey, 2))[2:].zfill(int(len(initialKey)/4))
    elif (initialKeyBase!=16):
        raise Exception('The input data should be in either base 16 or 2')
        
   

    #key Generation for all rounds
    roundKeys = KS.genKeyState(initialKey, initialKeyBase)
    rounds=len(roundKeys)-1
    
    #implementing dictionaries instead of calculation and opening them in a high level function to reduce compuation time
    mult = js.load(open("MulTable.json"))
    ISBox = js.load(open("InverseS-Box.json"))
    
    #cutting the input into elements with each of 128 bit representing a state after removing the last 8 bits which
    #represent the number of bits padded
    paddedLen=inputData[-8:]
    inputData=inputData[:-8]
    numStates=int(len(inputData)/128)
    states=[inputData[i*128:(i+1)*128] for i in range(numStates)]
    
    plaintext=""
    counter=1
    
    for state in states:
        if(showRounds):
            print(f"State {counter}")
            counter+=1
        plaintext+=decrypt_state(state,roundKeys,mult,ISBox,cipherDataBase=2,showRounds=showRounds,rounds=rounds)        
        if(showRounds):
            print("")
            
    plaintext=unpad(plaintext,paddedLen)
    
    return hex(int(plaintext,2))[2:].zfill(int(len(plaintext)/4))


def encrypt_cbc(inputData,initialKey,IV,showRounds=False,inputDataBase=16,initialKeyBase=16):
    """
    This function deals with either hex or binary strings. The string does not need to form a complete state and will be padded if needed.
    
    inputs:
    
    IV:         A hex string of length 32 charcters
    inputData:  Hex or binary string arrays
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    showRounds: show the output of each layer
    inputDataBase:  either 2 or 16
    initialKeyBase: either 2 or 16
    
    outputs:
    ciphertext: This represents the encrypted data. 
    
    
    
    """
    
    if (inputDataBase==16):
        inputData = bin(int(inputData, 16))[2:].zfill(int(len(inputData)*4))
    elif (inputDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')
    
    if (initialKeyBase==2):
        initialKey = hex(int(initialKey, 2))[2:].zfill(int(len(initialKey)/4))
    elif (initialKeyBase!=16):
        raise Exception('The input data should be in either base 16 or 2')
        
    assert (len(IV)==32)
    prev=bin(int(IV,16))[2:].zfill(128)
    
    #keyGeneration for all rounds
    roundKeys = KS.genKeyState(initialKey, initialKeyBase)
    rounds=len(roundKeys)-1
    
    #implementing dictionaries instead of calculation and opening them in a high level function to reduce compuation time
    mult = js.load(open("MulTable.json"))
    SBox = js.load(open("S-Box.json"))
    
    #cutting the input into elements with each of 128 bit representing a state.
    inputData,paddedLen=pad(inputData)
    numStates=int(len(inputData)/128)
    states=[inputData[i*128:(i+1)*128] for i in range(numStates)]
    
    ciphertext=""
    counter=1
    
    for state in states:
        if(showRounds):
            print(f"State {counter}")
            counter+=1
        prev=encrypt_state(Addition(state,prev,inputDataBase=2,roundKeyBase=2),roundKeys,mult,SBox,inputDataBase=2,showRounds=showRounds,rounds=rounds)
        ciphertext+=prev
        if(showRounds):
            print("")
            
    return hex(int(ciphertext+paddedLen,2))[2:].zfill(int(len(ciphertext+paddedLen)/4))

def encrypt_ecb(inputData,initialKey,showRounds=False,inputDataBase=16,initialKeyBase=16):
    """
    This function deals with either hex or binary strings. The string does not need to form a complete state and will be padded if needed.
    
    inputs:
    
    inputData:  Hex or binary string arrays
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    showRounds: show the output of each layer
    inputDataBase:  either 2 or 16
    initialKeyBase: either 2 or 16
    
    outputs:
    ciphertext: This represents the encrypted data. 
    
    
    
    """
    
    if (inputDataBase==16):
        inputData = bin(int(inputData, 16))[2:].zfill(int(len(inputData)*4))
    elif (inputDataBase!=2):
        raise Exception('The input data should be in either base 16 or 2')
    
    if (initialKeyBase==2):
        initialKey = hex(int(initialKey, 2))[2:].zfill(int(len(initialKey)/4))
    elif (initialKeyBase!=16):
        raise Exception('The Key should be in either base 16 or 2')
        
    #key generation for all rounds
    roundKeys = KS.genKeyState(initialKey, initialKeyBase)
    rounds=len(roundKeys)-1
    
    #implementing dictionaries instead of calculation and opening them in a high level function to reduce compuation time
    mult = js.load(open("MulTable.json"))
    SBox = js.load(open("S-Box.json"))
    
    #cutting the input into elements with each of 128 bit representing a state.
    inputData,paddedLen=pad(inputData)
    numStates=int(len(inputData)/128)
    states=[inputData[i*128:(i+1)*128] for i in range(numStates)]
    
    ciphertext=""
    counter=1
    
    for state in states:
        if(showRounds):
            print(f"State {counter}")
            counter+=1
        ciphertext+=encrypt_state(state,roundKeys,mult,SBox,inputDataBase=2,showRounds=showRounds,rounds=rounds)
        if(showRounds):
            print("")
            
    return hex(int(ciphertext+paddedLen,2))[2:].zfill(int(len(ciphertext+paddedLen)/4))



def encrypt_img(img,initialKey,IV="00000000000000000000000000000000",mode="ECB",showImage=False,showRounds=False,initialKeyBase=16):
    """
    This function deals with grey scale image arrays. The array is flattened and converted into hex values. The array is then passed to the respective 
    function according to the mode. The image is also padded if there are incomplete states and the last byte contains the number of padded characters.
    
    inputs:
    img: Grey scale matrix
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    IV:         initialization vector that must be 32 hex characters
    mode:       either ECB or CBC
    showImage:  Show the image after encryption, this cannot be done outside this function as the array is padded and a byte is added at the end and is
                not in a matrix format
    showRounds: show the output of each layer
    initialKeyBase: either 2 or 16 
    
    outputs:
    
    cipherImg:  An array containing the encrypted hex values 
    row:        The number of rows in the original image
    col:        The number of cols in the original image
    """
    
    
    assert (len(IV)==32)
    
    row,col=img.shape
    imgArray=img.reshape(-1)
    imgHex=""
    for i in range(imgArray.shape[0]):
        imgHex+=hex(imgArray[i])[2:].zfill(2)

    if(mode=="ECB"):
        cipherImg=encrypt_ecb(imgHex,initialKey,showRounds=showRounds,inputDataBase=16,initialKeyBase=initialKeyBase)
    elif mode=="CBC":
        cipherImg=encrypt_cbc(imgHex,initialKey,IV,showRounds=showRounds,inputDataBase=16,initialKeyBase=initialKeyBase)
    else:
        raise Exception('Valid modes are "ECB" and "CBC" only') 
    
    if(showImage):
        
        # the added padding is removed to be able to show the image
        paddedLen=bin(int(cipherImg[-2:],16))[2:].zfill(8)
        ImgShown=bin(int(cipherImg[:-2],16))[2:].zfill(len(cipherImg[:-2])*4)
        ImgShown=unpad(ImgShown,paddedLen)
        ImgShown=np.array([int(ImgShown[i*8:(i+1)*8],2) for i in range(int(len(ImgShown)/8))]).reshape((row,col))
        plt.imshow(ImgShown)
        
    return cipherImg,row,col


def decrypt_img(cipherImg,initialKey,row,col,IV="00000000000000000000000000000000",mode="CBC",showImage=False,showRounds=False,initialKeyBase=16):
    """
    This function deals with encrypted arrays that results from encrypt_img functions.
    
    inputs:
    cipherImg:  Hex array
    initialKey: either a base 2 or 16 string of characters with length varying between 128,192,and 256 bits
    IV:         initialization vector that must be 32 hex characters
    mode:       either ECB or CBC
    showImage:  Show the image after encryption, this cannot be done outside this function as the array is padded and a byte is added at the end and is
                not in a matrix format
    showRounds: show the output of each layer
    initialKeyBase: either 2 or 16 
    
    outputs:
    
    clearImg:  The original image matrix with correct dimensions 
    """
    assert (len(IV)==32)
    if(mode=="ECB"):
        clearImg=decrypt_ecb(cipherImg,initialKey,showRounds=showRounds,inputDataBase=16,initialKeyBase=initialKeyBase)
    elif(mode=="CBC"):
        clearImg=decrypt_cbc(cipherImg,initialKey,IV,showRounds=showRounds,inputDataBase=16,initialKeyBase=initialKeyBase)
    else:
        raise Exception('Valid modes are "ECB" and "CBC" only') 
    
    clearImg=np.array([int(clearImg[i*2:(i+1)*2],16) for i in range(int(len(clearImg)/2))]).reshape((row,col))
    if(showImage):
        plt.imshow(clearImg)
    return clearImg