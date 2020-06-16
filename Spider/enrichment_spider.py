import requests
import re
import random
import math
from Util.util import *


def read_genes(bic_path, name_path):

    biclusters = read_bic(bic_path)
    genes_with_index = [item[1] for item in biclusters]

    name_list = get_name_list(name_path)
    genes_with_name = []
    for item in genes_with_index:
        new_item = ""
        for i in range(len(item)):
            new_item += name_list[item[i]]
            if i != len(item) - 1:
                new_item += "\n"

        genes_with_name.append(new_item)

    return genes_with_name

def get_random_ip():

    proxies_list = ['39.137.69.9:80', '117.87.177.247:9000', '47.106.96.233:80', '123.57.84.116:8118', '61.143.38.53:8118',
               '117.90.5.173:9000', '58.87.73.43:8118', '47.98.198.125:8118', '117.87.177.135:9000',
               '61.183.176.122:57210',
               '117.28.97.121:808', '60.6.241.72:808', '117.114.149.10:45801', '103.242.202.178:39015',
               '123.127.93.188:44399',
               '222.180.162.245:9999', '183.166.129.53:8080', '163.125.112.83:8118', '112.16.169.194:32576']
    proxies_list = ['http://' + item for item in proxies_list]
    proxy_ip = random.choice(proxies_list)
    proxies = {'http': proxy_ip}

    return proxies


def open_url(url, featurelist):

    # 使用代理
    proxies = get_random_ip()
    data = {"featurelist": featurelist}
    headers = {"user-agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.108 Safari/537.36"}


    result = requests.post(url, headers = headers, data = data, proxies = proxies, verify = False)


    return result

def is_enrichment(result, length):

    def switch_number(number):

        return number // 4 + 1

    form = re.search("GO_id	Genome_GO_count	Expressed_GO	Expected_expression	Status	Corrected_P	GO_desc"
               + "((.*\n)*)", result.text).group(1)
    infos = form.split("\n")
    infos = infos[1:-1]
    for info in infos:
        info_split = info.split("\t")
        if len(info_split) == 7 and float(info_split[5]) < 0.05 and float(info_split[2]) >= switch_number(length):
                return True
    return False

if __name__ == '__main__':
    bic_path = "C:\\Users\\王利\\PycharmProjects\\GAEBic\\data\\4X6\\BiGAE biclusters(0.03)(4X6).txt"
    name_path = "C:\\Users\\王利\\PycharmProjects\\GAEBic\\data\\mRNA nodes.name"

    url = "https://www.soybase.org/Enrichment_v2/GO_Enrichment.php?enrichmentdl=true"
    requests.packages.urllib3.disable_warnings()

    count = 0
    num = 1

    genes_list = read_genes(bic_path, name_path)

    for featurelist in genes_list:
        res = open_url(url, featurelist)
        feature_len = len(featurelist.split("\n"))
        if is_enrichment(res, feature_len):
            count += 1

        print("%dth genes is fin" % num)
        num += 1

    with open("C:\\Users\\王利\\PycharmProjects\\GAEBic\\data\\4X6\\BiGAE(0.03) enrichment rate(4X6).txt", "w") as fw:
        fw.write("enrichment numbers = %d" % count)
        fw.write("\n")
        fw.write("enrichment rate = %f" % (count / len(genes_list)))

    print(count)
    print(count / len(genes_list))

