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
  #DPID NAMES I DO have

  #For ME12 minus
  chamber_length = 18
  HV_channels = 18
  #ME 21 everything 
  endcap = 1 #plus ; we are not sure of this
  station = 2 #second station
  ring = 1 
  f = open("dpid_name_ME21_plus_again.txt","w")
  f.write("ME21 plus : dpid \t rhid \n")
  for i in range(1,chamber_length+1) : 
    for j in range(1,HV_channels+1):
      chamber =i
      if(j==6 or j==12 or j==18) :
        layer = 6
      else :
        layer = j%6
      hvsegment = int((j-1)/6)+1
    
      rhid = 1000000*endcap + 100000 * station + 10000*ring + 100 *chamber + 10*layer + hvsegment
      # assuming we are taking correct values of endcap for plus
      dp_name_string = "cms_csc_dcs_3:CSC_ME_P21_C"+(str(i)).zfill(2)+"_HV_V"+(str(j)).zfill(2)+"_VMON"
      query = "SELECT ID, DPNAME FROM CMS_CSC_PVSS_COND.DP_NAME2ID WHERE DPNAME='"+dp_name_string+"'"
      print(query)
      for row in cur.execute(query):
        print(row)
        print(row[0],"name : ",row [1])
        string_write_one =str(row[0])+"\t"+str(rhid) 
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

