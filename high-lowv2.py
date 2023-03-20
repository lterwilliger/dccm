import pandas as pd
import numpy as np
from tkinter import filedialog
import os
import random
import matplotlib.pyplot as plt
#### CURRENT PREREQS BEFORE RUNNING!!!!!! ####
####
# 1 You must convert all tab files -> csv in excel and also split the text->columns
# 2 You must have two folders: Wild type Folder and Variant type folder
#   2a These folders have nothing but the DCCM csv files in them
# 3 You must know the number of variants in the gene
####
###############################################

#### GLOBAL CONFIGURATION SETTINGS #### 
####                               ####
 # set this to True to bring up prompt to enter a filename manually, else defaults to filename + _TOP_COR.csv
rename_outfile = False
# set this to however many variants you want in the top gain and top loss
num_top = 10
# set this to True if you want an output file of all the variants, False to turn off
will_write_all_var = True
# set this to True if you have not already averaged and subtracted the matrices
will_subtract = True
# set this to True to produce three heatmaps: wild avg, var avg, and subtracted 
will_plot = True
####                               ####
#######################################

def main():
    gene_name = input("Enter the name of the gene: ")
    num_var = int(input(f"Enter the number of amino acids in {gene_name}: "))
    if will_subtract:
        print("Wild File Directory")
        wild_files = get_filenames_in_direct()
        wild_file_count = len(wild_files)
        print("Variant File Directory")
        var_files = get_filenames_in_direct()
        var_file_count = len(var_files)
        
        wild_array = add_matrices(wild_files,num_var)
        wild_avg_array = average_matrices(wild_array,wild_file_count)
        var_array = add_matrices(var_files,num_var)
        var_avg_array= average_matrices(var_array,var_file_count)
        sub_array = subtract_matrices(var_avg_array,wild_avg_array)
    else:
        sub_filename = get_filename_pick()
        sub_array = sort_subtracted_file(sub_filename)
        
    if rename_outfile:
        out_name = input("Enter the name of the output file (end it with .csv): ")
    else:
        out_name = gene_name + "_TOP_COR.csv"

    if will_plot:
        heatmap(wild_avg_array, "Wild Type DCCM")
        heatmap(var_avg_array, "Variant Type DCCM")
        heatmap(sub_array, "Variant - Wild Type DCCM")
    sorted_index_array, np_dccm_sub = sort_subtracted(sub_array)
    outfile = open(out_name,"w")
    write_out_top(sorted_index_array,np_dccm_sub,outfile,False)
    write_out_top(sorted_index_array,np_dccm_sub,outfile,True)
    print(f"{num_top} variants have been written to {out_name}")
    all_out_filename = gene_name + "_ALL_COR.csv"
    if will_write_all_var:
        write_out_all(all_out_filename,sorted_index_array,num_var,np_dccm_sub) 
    
#######################################
def write_out_top(sorted_index_array,np_dccm_sub,outfile,is_gain):
    # write the top ten in pairs
    var_range = (num_top * 4) - 1
    if is_gain:
        change = "Loss"
        start_range = 0
        step_range = 4
    else:
        change = "Gain"
        start_range = -2
        var_range = (var_range * -1) - 1
        step_range = -4
    outfile.write(f"Greatest {change} of correlation\n")
    outfile.write("AA1,AA2,Correlation Change\n")
    for i in range(start_range,var_range,step_range):
        rounded_val = int((np_dccm_sub[int(i/2)]) * 1000) / 1000
        val = f"{sorted_index_array[i] + 1},{sorted_index_array[i+1] + 1},{rounded_val}\n"
        outfile.write(val)
    
def write_out_all(all_out_filename,sorted_index_array,num_var,np_dccm_sub):
    all_outfile = open(all_out_filename,"w")
    all_outfile.write("AA1,AA2,Correlation Change\n")
    i = 0
    while i < len(sorted_index_array) - 1:
        rounded_val = int((np_dccm_sub[int(i/2)]) * 1000) / 1000
        val = f"{sorted_index_array[i] + 1},{sorted_index_array[i+1] + 1},{rounded_val}\n"
        all_outfile.write(val)
        i = i + 2
    print(f"{num_var} variants have been written to {all_out_filename}")
 
def add_matrices(file_array,num_var):
    head = ""
    for x in range(0,num_var):
        head += f"{x}\t"
    head_list = head.split("\t")
    head_list.pop(num_var)
    sum_data_frame = pd.read_csv(file_array[0],sep=",",names = head_list, skiprows=1,skipinitialspace=True)
    print("Adding arrays..")
    for i in range(1,len(file_array)):
        temp_data_frame = pd.read_csv(file_array[i],sep="," ,names = head_list,skiprows=1,skipinitialspace=True)
        res_data_frame = sum_data_frame.add(temp_data_frame, fill_value=0)
        sum_data_frame = res_data_frame
    #print(sum_data_frame)
    return sum_data_frame
    # num_var = len(sum_data_frame)

def average_matrices(df,file_count):
    df = df / file_count
    df = df * 10000
    df = df.round(1)
    df = df.floordiv(10)
    df = df / 1000
    print("Averaging Matrices...")
    #print(df)
    return df

def subtract_matrices(df1,df2): 
    df_result = df1 - df2
    print("Subtracting Matrices...")
    print(df_result)
    return df_result

def get_file_path():
    path = filedialog.askdirectory()
    return path

def get_filenames_in_direct():
    path = get_file_path()
    print(path)
    dir_list = os.listdir(path)
    print("Files and directories in '", path, "' :")
    for i in range(0,len(dir_list)):
        dir_list[i] = os.path.join(path,dir_list[i])
    return dir_list

def make_unique_indices(np_arr):
    # this seems weird but the sorting algo needs non duplicates and this is the method i went with
    # I want when x = y, y = x to be the same random num
    for x in range (0,len(np_arr)):
        val = (random.random() / 1000)
        np_arr[x] += val
        #print(np_arr[x])

def sort_subtracted(data_frame):
    #np_dccm_sub = data_frame.iloc[:, :].values # [:, :] => [rows, columns]
    #print(len(data_frame))
    np_dccm_sub = data_frame.to_numpy()
    make_unique_indices(np_dccm_sub)
    sorted_index_array = np.unravel_index(np.argsort(np_dccm_sub, axis=None, kind="stable"),np_dccm_sub.shape)
    np_dccm_sub = np_dccm_sub[sorted_index_array]
    sorted_index_array = np.ravel(sorted_index_array, order="F")
    #print(sorted_index_array)
    
    print(len(sorted_index_array))
    print(len(np_dccm_sub))
    return sorted_index_array,np_dccm_sub 

def sort_subtracted_file(sub_filename):
     ###### Reading data, argsorting ######
    data_frame = pd.read_csv(sub_filename, index_col=0) 
    return data_frame
    #####################################

def get_filename_pick():
    filename = filedialog.askopenfilename(initialdir = "/",
                                        title = "Select a File",
                                        filetypes = (("CSV files",
                                                    "*.csv*"),
                                                    ("all files",
                                                    "*.*")))
    return filename

def heatmap(df,title):
    plt.imshow(df , cmap = 'autumn' , interpolation = 'nearest')
    plt.title( title )
    plt.show()
if __name__ == "__main__":
    main()
