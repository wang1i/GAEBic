import numpy as np

np.seterr(divide='ignore',invalid='ignore')

def print_bic(biclusters, bic_path):
    count = 1
    fw = open(bic_path, "w")
    for bic in biclusters:
        fw.write("The %dth bicluster:\n" % count)
        fw.write("\trow_pattern:\t")
        for r in bic[0]:
            fw.write(str(r) + "\t")
        fw.write("\n")
        fw.write("\tcol_pattern:\t")
        for c in bic[1]:
            fw.write(str(c) + "\t")
        fw.write("\n")
        count += 1
    fw.close()

def read_bic(bic_path):
    row_pattern = []
    col_pattern = []
    biclusters = []
    fr = open(bic_path, "r")
    i = 0
    while True:
        line = fr.readline()
        if not line:
            break
        if i == 0:
            i += 1
            continue
        temp = line.strip().split("\t")
        if i == 1:
            row_pattern = [int(temp[r]) for r in range(1, len(temp))]
        if i == 2:
            col_pattern = [int(temp[c]) for c in range(1, len(temp))]
            biclusters.append([row_pattern, col_pattern])
        i += 1
        if i == 3:
            i = 0
    fr.close()

    return biclusters

def adjust_bic(biclusters):

    def get_col_len(elem):
        return len(elem[1])

    biclusters.sort(key=get_col_len)
    biclusters.reverse()

    return biclusters

def remove_duplicate_bic(biclusters):

    remove_index = set([])
    for i in range(len(biclusters) - 1):
        bic1 = biclusters[i]
        rows1 = bic1[0]
        cols1 = bic1[1]
        for j in range(i + 1, len(biclusters)):
            bic2 = biclusters[j]
            rows2 = bic2[0]
            cols2 = bic2[1]
            if sorted(rows1) == sorted(rows2) and sorted(cols1) == sorted(cols2):
                remove_index.add(i)

    remove_index = list(remove_index)
    new_biclusters = []
    for i in range(len(biclusters)):
        if i not in remove_index:
            new_biclusters.append(biclusters[i])

    return new_biclusters


def convert(measure_matrix):

    min_val = measure_matrix.min(1)
    max_val = measure_matrix.max(1)
    ranges = max_val - min_val
    num_col = len(measure_matrix[0])

    new_matrix = measure_matrix - np.transpose(np.tile(min_val, (num_col, 1)))
    new_matrix = new_matrix / np.transpose(np.tile(ranges, (num_col, 1)))

    return new_matrix

def compute_one_rate(row_pattern, col_pattern, raw_matrix):

    # 先计算非0覆盖率
    one_count = 0
    for r in row_pattern:
        for c in col_pattern:
            if raw_matrix[r][c] == 1:
                one_count += 1
    one_rate = one_count / (len(row_pattern) * len(col_pattern))

    return one_rate

def get_name_list(name_path):

    name_list = []
    with open(name_path) as fr:
        for line in fr.readlines():
            if line != "":
                name_list.append(line.strip().replace("LYMA_", "lyma."))

    return name_list

