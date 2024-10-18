# importing module
import cx_Oracle
from nominal_HV_function_ME11 import nominal_HV_function_ME11 
from datetime import datetime
#import db_config
# Create a table in Oracle database
# It takes input dpid names and find the time during which there was a trip corresponding to all dpid
# Furthe it uses the function nominal_HV_function to evaluate the nominal time period corresponding to a dpid value
import time

st_time = time.time()
try:
  con = cx_Oracle.connect("cms_csc_pvss_cond_r", "DQM_dcsread_pvss1", "cms_omds_adg")
  print(con.version)
 
  # Now execute the sqlquery
  cur = con.cursor()
 
  date_start_list = "12-AUG-18 12.22.13.000000 PM";
  date_end_list = "13-AUG-18 11.31.56.000000 AM";
  query_start = "SELECT ACTUAL_VMON, CHANGE_DATE FROM CMS_CSC_PVSS_COND.FWCAENCHANNEL WHERE CHANGE_DATE <  '"+str(date_start_list) + "'  AND DPID=753147  AND ACTUAL_VMON IS NOT NULL ORDER BY CHANGE_DATE DESC FETCH FIRST 1 ROWS ONLY" ; 
  query_middle = "SELECT ACTUAL_VMON, CHANGE_DATE FROM CMS_CSC_PVSS_COND.FWCAENCHANNEL WHERE DPID=753147 AND CHANGE_DATE between '"+date_start_list+"' AND '"+str(date_end_list)+"'   AND ACTUAL_VMON IS NOT NULL ORDER BY CHANGE_DATE ASC"; 
  #    query = "SELECT ACTUAL_VMON,CHANGE_DATE FROM CMS_CSC_PVSS_COND.FWCAENCHANNEL WHERE DPID="+dpid_value +"AND CHANGE_DATE between '"+ start_time_list[i]+"' and '"+ end_time_list[i]+"'"
  
  # query = "SELECT DPID, CHANGE_DATE,  VALUE FROM CMS_CSC_PVSS_TEST.CSC_HV_V_DATA WHERE DPID=118018 AND CHANGE_DATE LIKE '31-JUL-17%'"
  
  # the obtained row is a list and itse element can be accessed similar to list elements

  dt_start = datetime.strptime(date_start_list, '%d-%b-%y %I.%M.%S.%f %p')
  dt_end = datetime.strptime(date_end_list, '%d-%b-%y %I.%M.%S.%f %p')

  # Convert the datetime object back to the desired string format
  m_start = dt_start.strftime('%Y-%m-%d %H:%M:%S.%f')
  m_end = dt_end.strftime('%Y-%m-%d %H:%M:%S.%f')

  value_list = []
  date_list = []
  num_medium = 0
  for row in cur.execute(query_start):
    value_list.append(row[0]) 
    date_list.append(row[1])
    print(" start time and HV ",row[0], " time ",row[1])
  for row in cur.execute(query_middle):
    value_list.append(row[0])
    date_list.append(row[1])
    num_medium=+1
    print(" between time and HV ",row[0], " time ",row[1])
 
   # Prepare the data
  data = []
  if(len(date_list)==1):
    data.append((str(m_start), str(m_end), str(value_list[0])))
  if(len(date_list)==2):
    data.append((str(m_start), str(date_list[1]), str(value_list[0])))
    data.append((str(date_list[1]), str(m_end), str(value_list[1])))
  if(len(date_list)>2):
    data.append((str(m_start), str(date_list[1]), str(value_list[0])))
    for i in range(1,len(date_list) - 1):
        data.append((str(date_list[i]), str(date_list[i+1]), str(value_list[i])))
    data.append((str(date_list[-1]), str(m_end), str(value_list[-1])))

  print("data from here \n")
  print(data)

# Write data to a text file
  with open("output.txt", "w") as f:
    f.write("Start time\tEnd time\tHV value\n")
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
