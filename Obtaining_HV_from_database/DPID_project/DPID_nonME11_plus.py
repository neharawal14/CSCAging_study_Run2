# importing module
import cx_Oracle
#import db_config
# Create a table in Oracle database
try:
  con = cx_Oracle.connect("cms_csc_pvss_cond_r", "DQM_dcsread_pvss1", "cms_omds_adg")
  
  print(con.version)
 
  # Now execute the sqlquery
  cur = con.cursor()
# Chamber configurations for different chambers (ME12, ME13, etc.)
chamber_configurations = [
    {"name": "ME12", "chamber_count": 36, "hv_channels": 18, "endcap": 2, "station": 1, "ring": 2, "query_str": "ME_P12", "system" : "cms_csc_dcs_3"},
    {"name": "ME13", "chamber_count": 36, "hv_channels": 18, "endcap": 2, "station": 1, "ring": 3, "query_str": "ME_P13","system" : "cms_csc_dcs_3"},
    {"name": "ME21", "chamber_count": 18, "hv_channels": 18, "endcap": 2, "station": 2, "ring": 1, "query_str": "ME_P21","system" : "cms_csc_dcs_3"},
    {"name": "ME31", "chamber_count": 18, "hv_channels": 18, "endcap": 2, "station": 3, "ring": 1, "query_str": "ME_P31","system" : "cms_csc_dcs_3"},
    {"name": "ME41", "chamber_count": 18, "hv_channels": 18, "endcap": 2, "station": 4, "ring": 1, "query_str": "ME_P41","system" : "cms_csc_dcs_3"},
    {"name": "ME22", "chamber_count": 36, "hv_channels": 30, "endcap": 2, "station": 2, "ring": 2, "query_str": "ME_P22","system" : "cms_csc_dcs_3"},
    {"name": "ME32", "chamber_count": 36, "hv_channels": 30, "endcap": 2, "station": 3, "ring": 2, "query_str": "ME_P32","system" : "cms_csc_dcs_3"},
    {"name": "ME42", "chamber_count": 36, "hv_channels": 30, "endcap": 2, "station": 4, "ring": 2, "query_str": "ME_P42","system" : "cms_csc_dcs_3"}
#    {"name": "ME12", "chamber_count": 36, "hv_channels": 18, "endcap": 1, "station": 1, "ring": 2, "query_str": "ME_M12", "system" : "cms_csc_dcs_2"}, 
#    {"name": "ME21", "chamber_count": 36, "hv_channels": 30, "endcap": 1, "station": 2, "ring": 1, "query_str": "ME_P21","system" : "cms_csc_dcs_3"},
    ]

for config in chamber_configurations:
        chamber_name = config["name"]
        HV_channel = config["hv_channels"]
        endcap = config["endcap"]
        station = config["station"]
        ring = config["ring"]
        chamber_number = config["chamber_count"]
        query_str = config["query_str"]
        system_str = config["system"]

        if(endcap==1):
            sign = "plus"
        elif(endcap==2):
            sign = "minus"


        # Open the file to write results for this chamber
        output_filename = f"dpid_name_{chamber_name}_{sign}_newOne.txt"
 
	    f = open(open_string,"w")
		statement_first = chamber_name +f"{sign}\tdpid\trhid\n"
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
    	    dp_name_string = "cms_csc_dcs_3:CSC_"+query_str[num]+"_C"+(str(i)).zfill(2)+"_HV_V"+(str(j)).zfill(2)+"_VMON"
    	    query = "SELECT ID, DPNAME FROM CMS_CSC_PVSS_COND.DP_NAME2ID WHERE DPNAME='"+dp_name_string+"'"
    	    print(query)
    	    for row in cur.execute(query):
    	      print(row)
    	      print(row[0],"name : ",row [1])
    	      string_write_one =f"{row[0]}\t{rhid}\t{row[1]}" 
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

