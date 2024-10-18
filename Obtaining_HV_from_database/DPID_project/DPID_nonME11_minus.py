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
  #ME11 details

  #chambers_name = ["ME21","ME22", "ME31","ME32", "ME41", "ME42"]
  #chamber_name_list = ["ME21", "ME31", "ME41"]
  #chamber_name_list = ["ME21", "ME31", "ME41"]
  #chamber_name_list = ["ME12", "ME13","ME21","ME31","ME41","ME22","ME32","ME42"]
  chamber_name_list = ["ME22"]
  #chamber_number_list = [36,36,18,18,18,36,36,36]
  chamber_number_list = [36]
  #HV_channels_list = [18,18,18,18,18,18,30,30,30]
  HV_channels_list = [30]
  #endcap_list = [2,2,2,2,2,2,2,2]
  endcap_list = [2]
  #station_list = [1,1,2,3,4,2,3,4]
  station_list = [2]
  #ring_list = [2,3,1,1,1,2,2,2]
  ring_list = [2]
  #query_str = ["ME_P22","ME_P32", "ME_P42"]
  #query_str = ["ME_M12","ME_M13","ME_M21","ME_M31","ME_M41","ME_M22","ME_M32","ME_M42"]
  query_str = ["ME_M22"]
  for num in range(len(chamber_name_list)):
    chamber_name = chamber_name_list[num]
    HV_channel = HV_channels_list[num]
    endcap = endcap_list[num]
    station = station_list[num]
    ring = ring_list[num]
    chamber_number = chamber_number_list[num]
   # chamber_length = 18
    #HV_channels = 18
  #ME 21 everything 
#  endcap = 1 #plus ; we are not sure of this
#  station = 2 #second station
#  ring = 1 
    open_string = "dpid_name_"+chamber_name+"_minus_new.txt"
    f = open(open_string,"w")
    statement_first = chamber_name +"minus\tdpid\trhid\n"
    f.write(statement_first)
    for i in range(1,chamber_number+1) : 
      for j in range(1,HV_channel+1):
        chamber =i
        if(j==6 or j==12 or j==18 or j==24 or j==30) :
        #if(j==6 or j==12 or j==18) :
          layer = 6
        else :
          layer = j%6
        hvsegment = int((j-1)/6)+1
      
        rhid = 1000000*endcap + 100000 * station + 10000*ring + 100 *chamber + 10*layer + hvsegment
        # assuming we are taking correct values of endcap for plus
        dp_name_string = "cms_csc_dcs_2:CSC_"+query_str[num]+"_C"+(str(i)).zfill(2)+"_HV_V"+(str(j)).zfill(2)+"_VMON"
        query = "SELECT ID, DPNAME FROM CMS_CSC_PVSS_COND.DP_NAME2ID WHERE DPNAME='"+dp_name_string+"'"
        print(query)
        for row in cur.execute(query):
          print(row)
          print(row[0],"name : ",row [1])
          string_write_one =str(row[0])+"\t"+str(rhid) 
          f.write(string_write_one)
          f.write("\n")
    f.close()
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

