# importing module
import cx_Oracle
#import db_config
# Create a table in Oracle database
try:
  con = cx_Oracle.connect("cms_csc_pvss_cond_r", "DQM_dcsread_pvss1", "cms_omds_adg")
  
  print(con.version)
 
  # Now execute the sqlquery
  cur = con.cursor()
 
  # Open the text file and find start and end of the run 
  # output is stored in the lists
  #dpid names are read from file    
  #chambers_name = ["ME21","ME22", "ME31","ME32", "ME41", "ME42"]
  #chamber_name_list = ["ME22", "ME32", "ME42"]
  #query_str = ["ME_P22","ME_P32", "ME_P42"]
  input_file = "/afs/cern.ch/work/n/nrawal/sql_access/DPID_project/ME11_HV_Mapping_CAEN.csv"
  f = open(input_file,"r")
  lines = f.readlines()[1:]
  open_string = "dpid_name_ME11_Run3.txt"
  f = open(open_string,"w")
  statement_first = "ME11\tdpid\trhid\n"
  f.write(statement_first)

  station = 1
  ring = 1
  for line in lines:
       split_line = line.split(";")
       dpid_name = split_line[0]
       chamber_name = split_line[1]
       chamber_endcap = chamber_name[2:3]
       chamber_number = chamber_name[7:9]
       hv_channel = chamber_name[10]
       
       if(chamber_endcap=="+") :
           endcap = 1
       if(chamber_endcap=="-") :
           endcap =2 
       layer = hv_channel
       rhid = 1000000*endcap + 100000 * station + 10000*ring + 100 *int(chamber_number) + 10*int(layer) + int(hv_channel)
       new_rhid = int((rhid - rhid%10)/10)
       #print("chamber name ",chamber_name, " rhid ",rhid, " new rhid ",new_rhid) 
       dp_name_string = dpid_name
       query = "SELECT ID, DPNAME FROM CMS_CSC_PVSS_COND.DP_NAME2ID WHERE DPNAME='"+dp_name_string+"'"
       print(query)
       for row in cur.execute(query):
          print(row)
          print(row[0],"name : ",row [1])
          string_write_one =str(row[0])+"\t"+str(new_rhid) 
          f.write(string_write_one)
          f.write("\n")

  print("Table Created successfully")

except cx_Oracle.DatabaseError as e:
  print("There is a problem with Oracle", e)

# by writing finally if any error occurs
# then also we can close the all database operation

finally:
  if cur:
    cur.close()
  if con:
    con.close()

