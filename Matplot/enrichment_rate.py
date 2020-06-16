import matplotlib.pyplot as plt
import numpy as np

algorithm_name = ['Spectral Biclustering', 'Bibit', 'Bimax', 'GAEBic(0.03)', 'GAEBic(0.04)', 'GAEBic(0.05)']
rate_44 = [0.0116, 0.4302, 0.4608, 0.7114,  0.5637, 0.5036]
rate_46 = [0.0244, 0.3679, 0.4672, 0.8305,  0.5526, 0.4855]

#准备X轴，宽带为0.3
bar_width = 0.3
x_44 = list(range(len(algorithm_name)))
x_46 = [i + bar_width for i in x_44]
x_ticks = [i - 0.15 for i in x_46]

#绘图
rect_1 = plt.bar(x_44,rate_44,width=bar_width,label = '4X4',color = 'Green')
rect_2 = plt.bar(x_46,rate_46,width=bar_width,label = '4X6',color = 'Orange')

#图例
plt.legend(loc = "upper left")

#X轴的刻度
plt.xticks(x_ticks,algorithm_name)

#坐标
for rect_1s in rect_1:
    height1 = rect_1s.get_height()
    plt.text(rect_1s.get_x() + rect_1s.get_width()/2, 1.01*height1, "%2.0f" % (height1 * 100) + "%", ha='center')

for rect_2s in rect_2:
    height2 = rect_2s.get_height()
    plt.text(rect_2s.get_x() + rect_2s.get_width()/2, 1.01*height2, "%2.0f" % (height2 * 100) + "%", ha='center')

plt.show()



