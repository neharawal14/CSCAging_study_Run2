import csv
import os
from datetime import datetime
debug = False
year = 2016
def split_csv(path_csv, path_run_list):
    f = open(path_csv)
    print("file opening" , path_csv)
    output_string = path_run_list+f"pressurecsc_{year}.h"
    output_file = open(output_string,"w") 
    lines = f.readlines()[3:]
    first_line = "double getpressure"+str(year)+"(UInt_t time){\n"
    output_file.write(first_line)
    print(" output file ", first_line)
   # output_file.write("run\tstart_time\tend_time\n") 
    for line in lines:
            #print(line)
            line = line.strip()
            # Check your delimiter
            splitted_line = line.split(";")
            time_value = splitted_line[0]
            pressure_value = splitted_line[1]
            #print("time value ",time_value) 
            time_value = time_value.replace('"','')
            time_value = int(time_value) 
            pressure_value = pressure_value.replace('"','')
            pressure_value = round(float(pressure_value),3)


            if(pressure_value ==-100) :
                 print("time : pressure : ",time_value, " : ",pressure_value)
                 continue
            #print(time_value, " ", pressure_value)            
            #print("time ", int(time_value))
            #print( " pressure ",float(pressure_value))
            new_string = "else if(time < "+str(time_value/1000)+") return "+str(pressure_value)+";\n"
            output_file.write(new_string)
    
    output_file.write("}")
 
if __name__ == '__main__':
    curr_dir = os.getcwd()
    path_csv = f"{curr_dir}/Barometric_pressures_2016.csv"
    #path_run_list = f"{curr_dir}/../Src/"
    path_run_list = f"{curr_dir}/"
    print(" output path ", path_run_list)
    split_csv(path_csv,path_run_list)
