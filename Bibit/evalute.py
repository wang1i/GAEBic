from Bibit.bibit import *
from Util.util import *

def compute_mean_euclid_distance(bic_path, measure_matrix_path):

    sum_dist = 0
    count = 0
    biclusters = read_bic(bic_path)
    measure_matrix = np.loadtxt(measure_matrix_path)
    measure_matrix = np.transpose(convert(measure_matrix))

    for bic in biclusters:
        col_pattern = bic[1]
        col_pairs = itertools.combinations(col_pattern, 2)
        for col_pair in col_pairs:
            temp_dist = np.linalg.norm(measure_matrix[:, col_pair[0]] - measure_matrix[:, col_pair[1]])
            if temp_dist < 0.1:
                sum_dist += temp_dist
                count += 1

    # temp_test = np.linalg.norm(measure_matrix[biclusters[599][0][0]] - new_matrix[biclusters[599][0][2]])

    print("sum_dist = " + str(sum_dist))
    print("count = " + str(count))
    # print("某一个双聚类的两行距离" + str(temp_test))

    return sum_dist / count


def result_without_merge():
    raw_matrix = np.loadtxt("/data/miRNA-4X4 matrix.txt")
    biclusters = bibit(raw_matrix, 4, 4, 4)
    print_bic(biclusters, "/data/Bibit biclusters.txt")


if __name__ == '__main__':
    bic_path = "/data/Bibit biclusters.txt"
    measure_matrix_path = "/data/4X4 embedding.txt"
    # Bibit的双聚类的平均欧式距离是：	0.06636658216603161
    print("Bibit的双聚类的平均欧式距离是：\t" + str(compute_mean_euclid_distance(bic_path, measure_matrix_path)))

