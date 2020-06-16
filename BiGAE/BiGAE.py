from Util.util import *
import numpy as np
import itertools


def col_euclid_distance(measure_matrix, col_pattern, c):

    sum_dist = 0
    count = 0
    for col in col_pattern:
        sum_dist += np.linalg.norm(measure_matrix[:, col] - measure_matrix[:, c])
        count += 1

    return sum_dist / count

def bigae(raw_matrix, measure_matrix, nr, nc):

    number = 1

    # 构建双聚类种子
    # 先对列进行分组，缩减遍历规模
    # 分为nc-1组，可以保证不会遗漏可能的双聚类结果
    group_len = len(raw_matrix[0]) // (nc - 1)
    left_len = len(raw_matrix[0]) - group_len * (nc - 1)
    all_col_pairs = []
    for k in range(nc - 1):
        all_col_pairs += list(itertools.combinations([i for i in range(k * group_len, (k + 1) * group_len)], 2))
    if left_len > 0:
        all_col_pairs += list(itertools.combinations([i for i in range((nc - 1) * group_len, len(raw_matrix[0]))], 2))

    bic_seeds = []
    for col_pair in all_col_pairs:
        col1 = col_pair[0]
        col2 = col_pair[1]
        col_pattern = [col1, col2]
        row_pattern = []
        for row in range(len(raw_matrix)):
            if raw_matrix[row][col1] == 1 and raw_matrix[row][col2] == 1:
                row_pattern.append(row)
        if len(row_pattern) >= nr and \
                np.linalg.norm(measure_matrix[:, col_pair[0]] - measure_matrix[:, col_pair[1]]) <= 0.05:
            bic_seeds.append([row_pattern, col_pattern])

    # 逐行添加表现好的列，保证非0覆盖率一直在85%以上
    biclusters = []
    temp_count = 0
    for bic_seed in bic_seeds:
        temp_count += 1

        row_pattern = bic_seed[0]
        col_pattern = bic_seed[1]

        # 添加列
        for c in range(len(raw_matrix[0])):
            if c in col_pattern:
                continue
            temp_dist = col_euclid_distance(measure_matrix, col_pattern, c)
            temp_rate = compute_one_rate(row_pattern, col_pattern + [c], raw_matrix)
            if temp_dist <= 0.05 and temp_rate >= 0.85:
                col_pattern.append(c)

        if len(col_pattern) >= nc:
            biclusters.append([row_pattern, col_pattern])
            print("Now, we find %dth bicluster" % number)
            number += 1
            print("\tNow bic_seed number" + str(temp_count) + "  /  " + "\tAll bic_seeds numbers" + str(len(bic_seeds)))

    biclusters = remove_duplicate_bic(biclusters)
    biclusters = adjust_bic(biclusters)

    return biclusters

if __name__ == '__main__':
    raw_matrix = np.loadtxt("/data/miRNA-mRNA matrix.txt")
    measure_matrix = np.loadtxt("/data/mRNA embedding.txt")
    measure_matrix = np.transpose(convert(measure_matrix))
    biclusters = bigae(raw_matrix, measure_matrix, 4, 6)
    print_bic(biclusters, "/data/4X6/BiGAE biclusters(0.05)(4X6).txt")




