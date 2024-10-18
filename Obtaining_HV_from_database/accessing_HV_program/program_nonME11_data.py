# importing module
import cx_Oracle
import argparse
#from nominal_HV_function_ME11 import nominal_HV_function_ME11 
from datetime import datetime
#import db_config
# Create a table in Oracle database
# It takes input dpid names and find the time during which there was a trip corresponding to all dpid
# Furthe it uses the function nominal_HV_function to evaluate the nominal time period corresponding to a dpid value
import time

parser = argparse.ArgumentParser()
parser.add_argument("chamber_name", help="chamber for which we need to derived HV")
parser.add_argument("sign_name", help="sign for the chamber")
parser.add_argument("year_name", help="year of the data")
args = parser.parse_args()
chamber = args.chamber_name
sign = args.sign_name
year = args.year_name
st_time = time.time()

# we should write a try and exception in case if the connection to server fails
try:
  con = cx_Oracle.connect("cms_csc_pvss_cond_r", "DQM_dcsread_pvss1", "cms_omds_adg")
  print(con.version)
 
  # Now execute the sqlquery
  # To get inside sqlserver
  cur = con.cursor()

  # We will look for all the CMS global runs ,and make arrays of run number, start time and end time
  #f = open("/afs/cern.ch/work/n/nrawal/sql_access/2022_trips/duration_list_ME11.txt","r")
  input_file_name = "../input_golden_lumi_files/duration_"+year+"_goldenjson.txt" 
  f = open(input_file_name,"r")
  lines = f.readlines()[1:]
  run_nb_list = []
  start_time_list = []
  end_time_list = []
  for line in lines:
    line_list = line.split("\t")
    run_nb = line_list[0].strip()
    start_time = line_list[1].strip()
    end_time = line_list[2].strip()
    run_nb_list.append(run_nb)
    start_time_list.append(start_time)
    end_time_list.append(end_time) 

  # Output file in which we will write the HV information for with start and end time
  input_text_file = "../final_HV_info_files/HV_info_golden_"+chamber+"_"+sign+"_"+year+"_asc.txt"
  with open(input_text_file, "w") as f:
    f.write("dpid\trhid\tRunNb\tStart time\tEnd time\tHVvalue\n")

  # Opening the file containing dpid and rhid of channels and read them into arrays
  
  #open_file_string = "/afs/cern.ch/work/n/nrawal/sql_access/DPID_project/dpid_name_ME11_plus_run3.txt"
  open_file_string = "../DPID_project/dpid_name_"+chamber+"_"+sign+".txt"
  f_dpid = open(open_file_string, "r")
  f_dpid_values = f_dpid.readlines()[1:]

  # for each channel we look over all the runs and find the HV corresponding to individual runs
  for value_dpid in f_dpid_values:
    dpid_value_string = value_dpid.split("\t")
    dpid_value = dpid_value_string[0]
    rhid_value = dpid_value_string[1].strip()
    for i in range(len(run_nb_list)):
      date_start_list = start_time_list[i]
      date_end_list = end_time_list[i]
      query_start = "SELECT VALUE,CHANGE_DATE FROM CMS_CSC_PVSS_COND.CSC_HV_V_DATA WHERE CHANGE_DATE <  '"+str(date_start_list) + "'  AND DPID="+dpid_value+" AND VALUE IS NOT NULL ORDER BY CHANGE_DATE DESC FETCH FIRST 1 ROWS ONLY" ; 
      query_middle = "SELECT VALUE,CHANGE_DATE FROM CMS_CSC_PVSS_COND.CSC_HV_V_DATA WHERE DPID="+dpid_value+" AND CHANGE_DATE between '"+date_start_list+"' AND '"+str(date_end_list)+"' AND VALUE IS NOT NULL ORDER BY CHANGE_DATE ASC"; 
      #    query = "SELECT ACTUAL_VMON,CHANGE_DATE FROM CMS_CSC_PVSS_COND.FWCAENCHANNEL WHERE DPID="+dpid_value +"AND CHANGE_DATE between '"+ start_time_list[i]+"' and '"+ end_time_list[i]+"'"
      
      # query = "SELECT DPID, CHANGE_DATE,  VALUE FROM CMS_CSC_PVSS_TEST.CSC_HV_V_DATA WHERE DPID=118018 AND CHANGE_DATE LIKE '31-JUL-17%'"
      
      # the obtained row is a list and itse element can be accessed similar to list elements

      #print(" start ", date_start_list, type(date_start_list))
      #print(" end ", date_end_list, type(date_end_list))
      dt_start = datetime.strptime(date_start_list, '%d-%b-%y %I.%M.%S.%f %p')
      dt_end = datetime.strptime(date_end_list, '%d-%b-%y %I.%M.%S.%f %p')

      # Convert the datetime object back to the desired string format
      m_start = dt_start.strftime('%Y-%m-%d %H:%M:%S')
      m_end = dt_end.strftime('%Y-%m-%d %H:%M:%S')

      value_list = []
      date_list = []
      num_medium = 0
      for row in cur.execute(query_start):
        value_list.append(row[0]) 
        date_list.append(row[1])
        #print(" start time and HV ",row[0], " time ",row[1])
      for row in cur.execute(query_middle):
        value_list.append(row[0])
        date_value_str = str(row[1])
        try:
          date_value = datetime.strptime(date_value_str,'%Y-%m-%d %H:%M:%S.%f')
        except:
          date_value = datetime.strptime(date_value_str,'%Y-%m-%d %H:%M:%S')
        m_value = datetime.strftime(date_value, '%Y-%m-%d %H:%M:%S')
        date_list.append(m_value)
        num_medium=+1
        #print(" between time and HV ",row[0], " time ",row[1])
     
       # Prepare the data
      data = []
      if(len(date_list)==1):
        data.append((dpid_value, rhid_value, run_nb_list[i], str(m_start), str(m_end), str(value_list[0])))
      if(len(date_list)==2):
        data.append((dpid_value, rhid_value,  run_nb_list[i], str(m_start), str(date_list[1]), str(value_list[0])))
        data.append((dpid_value, rhid_value, run_nb_list[i], str(date_list[1]), str(m_end), str(value_list[1])))
      if(len(date_list)>2):
        data.append((dpid_value, rhid_value, run_nb_list[i], str(m_start), str(date_list[1]), str(value_list[0])))
        for j in range(1,len(date_list) - 1):
            data.append((dpid_value, rhid_value, run_nb_list[i], str(date_list[j]), str(date_list[j+1]), str(value_list[j])))
        data.append((dpid_value, rhid_value, run_nb_list[i], str(date_list[-1]), str(m_end), str(value_list[-1])))

      #print("data from here \n")
      #print(data)

# Write data to a text file
      with open(input_text_file, "a") as f:
        for d in data:
         f.write("\t".join(d) + "\n")

  print("Data written to output.txt!")  

  print("Table Created successfully") 
except cx_Oracle.DatabaseError as e: 
  print("There is a problem with Oracle", e) # by writing finally if any error occurs # then also we can close the all database operation finally: 
if con: 
  con.close() 
  en_time = time.time() 
print("start time of program", st_time) 
print("end time of program", en_time)
