# -*- coding: utf-8 -*-
from CIKM_TOOLS import *
import matplotlib
import matplotlib.pyplot as plt


data_folder = '../data/'

#set_list = ['train','testA','testB']
#size_list = [10000,2000,2000]

set_list = ['testA']
size_list = [2000]

for set_name,set_size in zip(set_list,size_list):
    output_file = data_folder + set_name +  '_ubyte.txt'
    f = open(output_file, "w")
    f.close()
    Img_ind = 0
    input_file = data_folder + set_name +'.txt'
    with open(input_file) as f:
        for content in f:
            Img_ind = Img_ind +1
            line = content.split(',')
            title = line[0] + '    '+line[1]
            data_write = np.asarray(line[2].strip().split(' ')).astype(np.ubyte)
            data_write = (data_write + 1).astype(np.ubyte)
            if data_write.max()>255:
                print('too large')
            if data_write.min()<0:
                print('too small')
            f = open(output_file, "a")
            f.write(data_write.tobytes())
            f.close()

targetDir = '../figure'
for set_name, N_pic in zip(set_list, size_list):
    input_file = data_folder + set_name + '_ubyte.txt'
    fig = plt.figure(frameon = False)
    for pic_id1 in range(1, N_pic + 1):
        print('transforming ' + set_name + ': ' + str(pic_id1).zfill(5))
        for T_id in range(1, 16):
            for H_id in range(1, 5):
                FAIL_CORNER = 0
                data_mat1 = read_data(input_file, pic_id1, T_id, H_id)

                d_min = np.min(data_mat1)
                d_max = np.max(data_mat1)

                for i in range(len(data_mat1)):
                    for j in range(len(data_mat1[i])):
                        data_mat1[i][j] = 1.0 * (data_mat1[i][j] - d_min) / (d_max - d_min + 0.01) * 255
                        #data_mat1[i][j] = data_mat1[i][j]
                    #print len(data_mat1), len(data_mat1[0])
                fig.clf()
                ax = fig.gca()
                ax.set_frame_on(False)
                ax.margins(0, 0)
                #pm = ax.pcolormesh(range(0, 101), range(0, 101), data_mat1)
                pm = ax.pcolormesh(range(0, 101), range(0, 101), data_mat1, cmap = 'gray')
                #fig.colorbar(pm)
                #plt.plot(data_mat1)
                plt.xticks([])
                plt.yticks([])
                #plt.show()
                fig.savefig(targetDir + '/' + str(H_id - 0.5) + '/' + set_name + str(pic_id1).zfill(5)+ '.T' + str(T_id).zfill(2) + '.png', bbox_inches = 'tight', pad_inches = 0)