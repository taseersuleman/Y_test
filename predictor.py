import numpy as np
import pandas as pd
import streamlit as st
from Bio import SeqIO
import math
import csv
import pickle
from PIL import Image
from sklearn.preprocessing import StandardScaler

try:
    str22 = "CGGGGGACAGGAUAAACUGUUACGCAUAUAUGACUUGAACA"
    st.subheader("Kindly click the EXAMPLE button for sample RNA sequence")
    if st.button('Example'):
        st.write(str22)


    def seqToMat(seq):
        encoder = ['a', 'c', 'u', 'g']

        lent = len(seq)
        n = int(math.ceil(math.sqrt(lent)))
        seqMat = [[0 for x in range(n)] for y in range(n)]
        i = 0
        seqiter = 0
        for i in range(n):
            j = 0
            for j in range(n):
                if seqiter < lent:
                    try:
                        aa = int(encoder.index(seq[seqiter]))
                    except ValueError:
                        #                    print(seq)
                        #                    print(seqiter)
                        exit(0)
                    else:
                        seqMat[i][j] = aa
                    seqiter += 1
        return seqMat

        ###frequencyVector###


    def frequencyVec4(seq):
        encoder = ['a', 'c', 'u', 'g']
        fv = [0 for x in range(4)]
        i = 0
        for i in range(4):
            fv[i] = seq.count(encoder[i])
        #    print('FV')
        # print(fv)
        return fv


    def frequencyVec16(encoding1, seq):
        fv = [0 for x in range(16)]
        i = 0
        for i in range(16):
            fv[i] = encoding1.count(i)
        # print('FV')
        # print(fv)
        return fv


    def frequencyVec64(encoding2, seq):
        fv = [0 for x in range(64)]
        i = 0
        for i in range(64):
            fv[i] = encoding2.count(i)
        # print('FV')
        # print(fv)
        return fv
        ###AAPIV Matrices###


    def AAPIV4(seq):
        encoder = ['a', 'c', 'u', 'g']
        apv = [0 for x in range(4)]
        i = 1
        sum = 0
        for i in range(4):
            j = 0
            for j in range(len(seq)):
                if seq[j] == encoder[i]:
                    sum = sum + j + 1
            apv[i] = sum
            sum = 0
        # print('AAPIV')
        # print(apv[1:] + apv[0:1])
        return apv[1:] + apv[0:1]


    def AAPIV16(encoding1, seq):
        apv = [0 for x in range(16)]
        i = 1
        sum = 0
        for i in range(16):
            j = 0
            for j in range(len(encoding1)):
                if encoding1[j] == i:
                    sum = sum + j + 1
            apv[i] = sum
            sum = 0
        # print('AAPIV')
        # print(apv[1:] + apv[0:1])
        return apv[1:] + apv[0:1]


    def AAPIV64(encoding2, seq):
        apv = [0 for x in range(64)]
        i = 1
        sum = 0
        for i in range(64):
            j = 0
            for j in range(len(encoding2)):
                if encoding2[j] == i:
                    sum = sum + j + 1
            apv[i] = sum
            sum = 0
        # print('AAPIV')
        # print(apv[1:] + apv[0:1])
        return apv[1:] + apv[0:1]


    def print2Dmat(mat):
        n = len(mat)
        i = 0
        strOut = ''
        for i in range(n):
            strOut = strOut + str(mat[i]) + '<br>'
        return strOut
        ####Prim matrices###


    def PRIM4(seq):
        encoder = ['a', 'c', 'u', 'g']
        prim = [[0 for x in range(4)] for y in range(4)]
        i = 0
        for i in range(4):
            aa1 = encoder[i]
            aa1index = -1
            for x in range(len(seq)):
                if seq[x] == aa1:
                    aa1index = x + 1
                    break
            if aa1index != -1:
                j = 0
                for j in range(4):
                    if j != i:
                        aa2 = encoder[j]
                        aa2index = 0
                        for y in range(len(seq)):
                            if seq[y] == aa2:
                                aa2index = aa2index + ((y + 1) - aa1index)
                        prim[i][j] = int(aa2index)
        # print('prim')
        # print(prim)
        return prim


    def cal_dibase_index(seq):
        dibase_index = {'aa': 1, 'ac': 2, 'au': 3, 'ag': 4, 'ca': 5, 'cu': 6, 'cc': 7, 'cg': 8, 'gg': 9, 'ga': 10, 'gc': 11,
                        'gu': 12, 'uu': 13, 'ua': 14, 'uc': 15, 'ug': 16}
        return dibase_index[seq]


    def PRIM16(encoding1):
        # print(encoding1)
        prim = [[0 for x in range(16)] for y in range(16)]
        i = 0
        # print(np.array(prim).shape)
        for i in range(len(prim)):
            # print('--------------------------------------')
            # print(encoding1[i],(i+1))
            # print('--------------------------------------')
            f = -1
            for j in range(len(encoding1)):
                # print(encoding1[j],(j+1))
                if (i + 1) == encoding1[j] and f == -1:
                    f = j
                    #  print('i',f,encoding1[f])
                    break
            for j in range(len(encoding1)):
                if f != -1 and f != j:
                    #                if i == 4 and encoding1[j] == 5:
                    #                    print(j+1,encoding1[j],f+1,(j+1)-(f+1)-0)
                    #   print(encoding1[f])
                    if encoding1[f] == (i + 1):
                        prim[i][encoding1[j] - 1] += (j) - (f)
        #               break
        # print('prim')
        # print(prim)
        return prim


    def first_index(encoding1):
        # print(encoding1)
        prim = [[0 for x in range(16)] for y in range(16)]
        i = 0
        print(np.array(prim).shape)
        results = np.array([0] * len(prim))
        for i in range(len(prim)):
            # print('--------------------------------------')
            # print(encoding1[i],(i+1))
            # print('--------------------------------------')
            for j in range(len(encoding1)):
                # print(encoding1[j],(j+1))
                if (i + 1) == encoding1[j]:
                    results[i] = j + 1
                    # print('i',f,encoding1[f])
                    break

        return results


    def cal_tribase_index(seq):
        if seq[0] == 'a':
            x = cal_dibase_index(seq[1] + seq[2])
            return (x + 0)
        elif seq[0] == 'c':
            x = cal_dibase_index(seq[1] + seq[2])
            return (x + 16)
        elif seq[0] == 'g':
            x = cal_dibase_index(seq[1] + seq[2])
            return (x + 32)
        elif seq[0] == 'u':
            x = cal_dibase_index(seq[1] + seq[2])
            return (x + 48)


    def PRIM64(encoding2):
        # print(encoding2)
        prim = [[0 for x in range(64)] for y in range(64)]
        i = 0
        # print(np.array(prim).shape)
        for i in range(len(prim)):
            # print('--------------------------------------')
            # print(encoding2[i],(i+1))
            # print('--------------------------------------')
            f = -1
            for j in range(len(encoding2)):
                # print(encoding2[j],(j+1))
                if (i + 1) == encoding2[j] and f == -1:
                    f = j
                    #  print('i',f,encoding2[f])
                    break
            for j in range(len(encoding2)):
                if f != -1 and f != j:
                    #                if i == 4 and encoding2[j] == 5:
                    #                    print(j+1,encoding2[j],f+1,(j+1)-(f+1)-0)
                    #   print(encoding2[f])
                    if encoding2[f] == (i + 1):
                        prim[i][encoding2[j] - 1] += (j) - (f)
        #               break
        # print('prim')
        # print(prim)
        return prim


    def first_index(encoding2):
        # print(encoding2)
        prim = [[0 for x in range(64)] for y in range(64)]
        i = 0
        # print(np.array(prim).shape)
        results = np.array([0] * len(prim))
        for i in range(len(prim)):
            # print('--------------------------------------')
            # print(encoding2[i],(i+1))
            # print('--------------------------------------')
            for j in range(len(encoding2)):
                # print(encoding[j],(j+1))
                if (i + 1) == encoding2[j]:
                    results[i] = j + 1
                    # print('i',f,encoding2[f])
                    break
        return results


    def rawMoments(mat, order):
        n = len(mat)
        rawM = []
        sum = 0
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    p = 0
                    for p in range(n):
                        q = 0
                        for q in range(n):
                            sum = sum + (((p + 1) ** i) * ((q + 1) ** j) * int(mat[p][q]))
                    rawM.append(sum)
                    sum = 0
        # print('rawmmoments')
        # print(rawM)
        return rawM


    def centralMoments(mat, order, xbar, ybar):
        n = len(mat)
        centM = []
        sum = 0
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    p = 0
                    for p in range(n):
                        q = 0
                        for q in range(n):
                            sum = sum + ((((p + 1) - xbar) ** i) * (((q + 1) - ybar) ** j) * mat[p][q])
                    centM.append(sum)
                    sum = 0
        # print('central moments')
        # print(centM)
        return centM


    def hahnMoments(mat, order):
        N = len(mat)
        hahnM = []
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    answer = hahnMoment(i, j, N, mat)
                    hahnM.append(answer)
        # print('hahn moments')
        # print(hahnM)
        return hahnM


    def hahnMoment(m, n, N, mat):
        value = 0.0
        x = 0
        for x in range(N):
            y = 0
            for y in range(N):
                value = value + (
                        mat[x][y] * (hahnProcessor(x, m, N)) * (hahnProcessor(x, n, N)))
        #    #print('value')
        #    #print(value)
        return value


    def hahnProcessor(x, n, N):
        return hahnPol(x, n, N) * math.sqrt(roho(x, n, N))


    def hahnPol(x, n, N):
        answer = 0.0
        ans1 = pochHammer(N - 1.0, n) * pochHammer(N - 1.0, n)
        ans2 = 0.0
        k = 0
        for k in range(n + 1):
            ans2 = ans2 + math.pow(-1.0, k) * ((pochHammer(-n, k) * pochHammer(-x, k) *
                                                pochHammer(2 * N - n - 1.0, k)))
        answer = ans1 + ans2
        return answer


    def roho(x, n, N):
        return gamma(n + 1.0) * gamma(n + 1.0) * pochHammer((n + 1.0), N)


    def gamma(x):
        return math.exp(logGamma(x))


    def logGamma(x):
        temp = (x - 0.5) * math.log(x + 4.5) - (x + 4.5)
        ser = 101.19539853003
        return temp + math.log(ser * math.sqrt(2 * math.pi))


    def pochHammer(a, k):
        answer = 1.0
        i = 0
        for i in range(k):
            answer = answer * (a + i)
        return answer


    def calcFV(seq):
        encoding1 = list()
        for i in range(0, len(seq) - 1, 2):
            encoding1.append(cal_dibase_index(seq[i] + seq[i + 1]))
        encoding2 = list()
        for i in range(0, len(seq) - 2, 3):
            encoding2.append(cal_tribase_index(seq[i] + seq[i + 1] + seq[i + 2]))

        el1 = len(encoding1)
        el12 = int(math.ceil(math.sqrt(el1)))
        new_enc1 = []
        n = 0
        for i in range(el12):
            temp = []
            for j in range(el12):
                try:
                    temp.append(encoding1[n])
                except:
                    temp.append(0)
                n += 1
            new_enc1.append(temp)

        encoding1_2 = new_enc1

        el1 = len(encoding2)
        el12 = int(math.ceil(math.sqrt(el1)))
        new_enc1 = []
        n = 0
        for i in range(el12):
            temp = []
            for j in range(el12):
                try:
                    temp.append(encoding2[n])
                except:
                    temp.append(0)
            new_enc1.append(temp)
            n += 1

        encoding2_2 = new_enc1

        fv = [0 for x in range(522)]

        fvIter = 0
        myMat = seqToMat(seq)

        # print('seq r m')
        myRawMoments = rawMoments(myMat, 3)
        for ele in myRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        xbar = myRawMoments[4]
        ybar = myRawMoments[1]

        # print('seq C m')
        myCentralMoments = centralMoments(myMat, 3, xbar, ybar)

        for ele in myCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq h m')
        myHahnMoments = hahnMoments(myMat, 3)
        for ele in myHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq f 4')
        myFrequencyVec4 = frequencyVec4(seq)
        for ele in myFrequencyVec4:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 r m')
        myRawMoments = rawMoments(encoding1_2, 3)
        for ele in myRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        xbar = myRawMoments[4]
        ybar = myRawMoments[1]

        # print('e1 c m')
        myCentralMoments = centralMoments(encoding1_2, 3, xbar, ybar)
        for ele in myCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 h m')
        myHahnMoments = hahnMoments(encoding1_2, 3)
        for ele in myHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 f 16')
        myFrequencyVec16 = frequencyVec16(encoding1, seq)
        for ele in myFrequencyVec16:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 r m')
        myRawMoments = rawMoments(encoding2_2, 3)
        for ele in myRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        xbar = myRawMoments[4]
        ybar = myRawMoments[1]

        # print('e2 c m')
        myCentralMoments = centralMoments(encoding2_2, 3, xbar, ybar)
        for ele in myCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 h m')
        myHahnMoments = hahnMoments(encoding2_2, 3)
        for ele in myHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 f 64')
        myFrequencyVec64 = frequencyVec64(encoding2, seq)
        for ele in myFrequencyVec64:
            fv[fvIter] = ele
            fvIter = fvIter + 1

            # print('seq p 4')
        myPRIM = PRIM4(seq)
        # print('seq r m')
        myPRIMRawMoments = rawMoments(myPRIM, 3)
        xbar2 = myPRIMRawMoments[4]
        ybar2 = myPRIMRawMoments[1]

        # print('seq c m')
        myPRIMCentralMoments = centralMoments(myPRIM, 3, xbar2, ybar2)
        for ele in myPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        for ele in myPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq h m')
        myPRIMHahnMoments = hahnMoments(myPRIM, 3)
        for ele in myPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 p 16')
        myPRIM = PRIM16(encoding1)

        # print('e1 r m')
        myPRIMRawMoments = rawMoments(myPRIM, 3)
        xbar2 = myPRIMRawMoments[4]
        ybar2 = myPRIMRawMoments[1]

        # print('e1 c m')
        myPRIMCentralMoments = centralMoments(myPRIM, 3, xbar2, ybar2)
        for ele in myPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 h m')
        myPRIMHahnMoments = hahnMoments(myPRIM, 3)
        for ele in myPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 p 64')
        myPRIM = PRIM64(encoding2)
        # print('e2 r m')
        myPRIMRawMoments = rawMoments(myPRIM, 3)
        xbar2 = myPRIMRawMoments[4]
        ybar2 = myPRIMRawMoments[1]

        # print('e2 c m')
        myPRIMCentralMoments = centralMoments(myPRIM, 3, xbar2, ybar2)
        for ele in myPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 h m')
        myPRIMHahnMoments = hahnMoments(myPRIM, 3)
        for ele in myPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq a 4')
        myAAPIV = AAPIV4(seq)
        for ele in myAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 a 16')
        myAAPIV = AAPIV16(encoding1, seq)
        for ele in myAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 a 64')
        myAAPIV = AAPIV64(encoding2, seq)
        for ele in myAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq p 4')
        myRPRIM = PRIM4(seq[::-1])

        # print('seq r m')
        myRPRIMRawMoments = rawMoments(myRPRIM, 3)
        xbar3 = myRPRIMRawMoments[4]
        ybar3 = myRPRIMRawMoments[1]

        # print('seq C m')
        myRPRIMCentralMoments = centralMoments(myRPRIM, 3, xbar3, ybar3)
        for ele in myRPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myRPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('seq h m')
        myRPRIMHahnMoments = hahnMoments(myRPRIM, 3)
        for ele in myRPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 p 16')
        myRPRIM = PRIM16(encoding1[::-1])

        # print('e1 r m')
        myRPRIMRawMoments = rawMoments(myRPRIM, 3)
        xbar3 = myRPRIMRawMoments[4]
        ybar3 = myRPRIMRawMoments[1]

        # print('e1 c m')
        myRPRIMCentralMoments = centralMoments(myRPRIM, 3, xbar3, ybar3)
        for ele in myRPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myRPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e1 h m')
        myRPRIMHahnMoments = hahnMoments(myRPRIM, 3)
        for ele in myRPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 p 64')
        myRPRIM = PRIM64(encoding2[::-1])

        # print('e2 r m')
        myRPRIMRawMoments = rawMoments(myRPRIM, 3)
        xbar3 = myRPRIMRawMoments[4]
        ybar3 = myRPRIMRawMoments[1]

        # print('e2 c m')
        myRPRIMCentralMoments = centralMoments(myRPRIM, 3, xbar3, ybar3)
        for ele in myRPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myRPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        # print('e2 h m')
        myRPRIMHahnMoments = hahnMoments(myRPRIM, 3)
        for ele in myRPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        myRAAPIV = AAPIV4(seq[::-1])
        for ele in myRAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        myRAAPIV = AAPIV16(encoding1, seq[::-1])
        for ele in myRAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        myRAAPIV = AAPIV64(encoding2, seq[::-1])
        for ele in myRAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1

        return fv


    def processAllStrings(fname):
        seqs = []
        allFVs = []
        with open(fname, 'r') as filehandle:
            for line in filehandle:
                currentPlace = line[:-1]
                seqs.append(currentPlace)
        allowed_chars = set('aguc')
        i = 0
        for seq in seqs:
            seq = seq.lower()
            print(str(i) + ': ' + seq)
            if seq != '':
                if set(seq).issubset(allowed_chars):
                    # print(len(seq))
                    allFVs.append(calcFV(seq))
                    i = i + 1
                else:
                    print('Invalid Sequence\n' + str(i))
                    i = i + 1
        return allFVs


    i = 0
    allFVs = []


    def input_seq():
        st.subheader("Input Sequence of any length")
        sequence1 = st.text_area("Sequence Input", height=20)
        abc = []
        if st.button("Submit"):

            abc = sequence1
            count = []
            keeper = []
            len1 = len(abc)
            len2 = len1 - 1

            for i in range(len1):
                if abc[i] == "U":
                    count.append(i)

            len3 = len(count)
            for i in range(len3):
                s = count[i]

                if s <= 20 and s <= len2:
                    n = len2 - s
                    m = 20 - s
                    # n = 20 - m
                    str1 = ("A" * m)
                    str2 = abc[s - s:s]

                    if n <= 20:
                        o = 20 - n
                        str4 = abc[s:s + n + 1]
                        str5 = ("A" * o)
                        str6 = "".join((str1, str2, str4, str5))
                        keeper.append(str6)
                    elif n > 20:
                        str4A = abc[s:s + 21]
                        str6A = "".join((str1, str2, str4A))

                        keeper.append(str6A)

                elif s > 20:
                    n1 = len1 - s
                    str7 = abc[s - 20:s + 1]
                    if n1 <= 20:
                        o1 = 20 - n1
                        o1 = o1 + 1
                        str9 = ("A" * o1)
                        str8 = abc[s:s + n1 - 1]

                        str10 = "".join((str7, str8, str9))

                        keeper.append(str10)
                    elif n1 > 20:
                        str8A = abc[s:s + 20]
                        str10A = "".join((str7, str8A))

                        keeper.append(str10A)

            # print("##############################################################################")
            klen = len(keeper)

            for i in range(klen):
                seq = keeper[i]
                allFVs.append(calcFV(seq.lower()))

            with open('./IISequence_FVs_for_test.csv', mode='w') as fvFile:
                fvWriter = csv.writer(fvFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                for fv in allFVs:
                    fvWriter.writerow(fv)

            np.random.seed(5)
            inputSize = 522
            outputcol = inputSize + 1
            # dataset = np.genfromtxt("IISequence_FVs_for_test.csv", delimiter=",", dtype=float)
            df = pd.read_csv("IISequence_FVs_for_test.csv", sep=',', header=None)
            W = df.iloc[:, :].values
            Y = df.iloc[:, -1].values
            # W = dataset[:, 0:inputSize]
            # Y = dataset[:, -1]
            std_scale = StandardScaler().fit(W)
            W = std_scale.transform(W)
            load_model = pickle.load(open('Y_site_model.pkl', 'rb'))
            pred = load_model.predict(W)
            output_proba = load_model.predict_proba(W)[:, 1]

            lno = len(output_proba)

            for i in range(lno):

                st.subheader("Site Number = ")
                st.write(count[i])
                st.subheader("Sequence")
                st.write(keeper[i])

                if output_proba[i] > 0.5:
                    st.info("Output =  ψ-site")
                else:
                    st.info("Output = Non-ψ-site")


    input_seq()


except:
  print("An exception occurred")
