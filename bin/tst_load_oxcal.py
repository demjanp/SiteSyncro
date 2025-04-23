import os
from sitesyncro.utils.fnc_oxcal import load_oxcal_data

if __name__ == '__main__':
	
	fname = "c:\\documents_synced\\SiteSyncro\\sitesyncro\\stage_2\\model.js"
	data = load_oxcal_data(fname)
	print()
	print(data)
	print()