import numpy as np
import operator
import itertools

from Util.util import *

def encode(raw_matrix, bitwise):
    row_numbers = len(raw_matrix)
    col_numbers = int(len(raw_matrix[0]) / bitwise)
    encode_matrix = np.zeros((row_numbers, col_numbers), dtype=np.int)
    for row in range(row_numbers):
        new_column = 0
        for col in range(0, col_numbers * bitwise, bitwise):
            temp_list = []
            for t in range(bitwise):
                temp_list.append(raw_matrix[row][col + t])
            value = 0
            for i in range(bitwise):
                value = value * 2 + temp_list[i]
            encode_matrix[row][new_column] = int(value)
            new_column += 1

    print("Now, Encoding has completed !")

    return encode_matrix

def get_pattern(row1, row2, bitwise):
    p_res = []
    col_pattern = []
    for index in range(len(row1)):
        num1 = int(row1[index])
        num2 = int(row2[index])
        temp = num1 & num2
        p_res.append(temp)
        res = bin(temp).replace("0b", "")
        while (len(res) < bitwise):                   #不足bitwise位元数，则用0填充至足量
            res = "0" + res
        for i in range(len(res)):
            if (res[i] == '1'):
                col_pattern.append(index * bitwise + i)
    return p_res, col_pattern

def bibit(raw_matrix, nr, nc, bitwise):

    encode_matrix = encode(raw_matrix, bitwise)
    print(encode_matrix.shape)

    biclusters = []
    bicluster_number = 1
    patterns = []

    rows_pairs = itertools.combinations([i for i in range(len(encode_matrix))], 2)

    for row_pattern in rows_pairs:
        row1 = encode_matrix[row_pattern[0]]
        row2 = encode_matrix[row_pattern[1]]

        p_res, col_pattern = get_pattern(row1, row2, bitwise)

        if p_res in patterns:
            continue
        patterns.append(p_res)
        # 列数不满足最低限度值或者p_res不为新，则退出这次工作，从下一个pair中继续寻找双簇
        if len(col_pattern) < nc:
            continue

        bic_candidate = [list(row_pattern), col_pattern]  # 候选的一个双簇

        for index in range(len(encode_matrix)):
            if index != row_pattern[0] and index != row_pattern[1]:
                p_temp = get_pattern(encode_matrix[index], p_res, bitwise)[0]
                if operator.eq(p_temp, p_res):
                    bic_candidate[0].append(index)

        if len(bic_candidate[0]) >= nr:
            print("Now, we find the %dth bicluster" % bicluster_number)
            # print(bic_candidate[0])
            bicluster_number += 1 # 895个簇
            biclusters.append(bic_candidate)

    biclusters = adjust_bic(biclusters)

    return biclusters



if __name__ == '__main__':
    raw_matrix = np.loadtxt("/data/miRNA-mRNA matrix.txt")
    biclusters = bibit(raw_matrix, 4, 6, 4)
    print_bic(biclusters, "/data/4X6/Bibit biclusters(4X6).txt")