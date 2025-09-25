# Autofill Google Form
# based on https://github.com/jamesshah/GoogleForm-AutoFill
# Bram VDB
# 14/10/2022

from requests import post
from time import sleep
from sys import argv

def send_response(url, data):
    """It takes google form url which is to be submitted and also data which is a list of data to be submitted in the form iteratively."""
    for d in data:
        try:
            post(url, data=d)
            print("Form Submitted.")
            sleep(5)
        except:
            print("Error Occured!")
# URL to the form you want to fill. formResponse should be used instead of viewform
url = 'https://docs.google.com/forms/d/1OdyJcOKBngdj8xOWg3s7w6065TfF6JeE4PiWuIrhLA0/formResponse'

# assign values based on arguments in command line
groupname  = argv[1]
name       = argv[2]
dist       = argv[3]
vel        = argv[4]
llen       = argv[5]
mass       = argv[6]
shoulders  = argv[7]
arm_upper  = argv[8]
arm_lower  = argv[9]
leg_upper  = argv[10]
leg_lower  = argv[11]
foot       = argv[12]
strength   = argv[13]

# define values for google form fiels
values = {
		 # groupname
		"entry.104868765": str(groupname),
         # name
        "entry.378928282": str(name),
         # distance
        "entry.948065567": str(dist),
         # velocity
        "entry.1829825032": str(vel),
         # body length
        "entry.1615981409":str(llen),
         # mass
        "entry.862288878": str(mass),
         # shoulder width
        "entry.265262110": str(shoulders),
         # upper arm
        "entry.328724874": str(arm_upper),
         # lower arm 
        "entry.1769227683": str(arm_lower),
         # upper leg 
        "entry.1309586461": str(leg_upper),
         # lower leg 
        "entry.253816664": str(leg_lower),
         # foot length 
        "entry.516261500": str(foot),
         # muscle force
        "entry.1463940260": str(strength)
        }

values_list = []  
values_list.append(values)

send_response(url, values_list)
