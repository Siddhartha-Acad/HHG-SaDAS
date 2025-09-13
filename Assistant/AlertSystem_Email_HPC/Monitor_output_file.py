"""
File: Monitor_output_file.py
Project: HHG-SaDAS
Code Description:


Original Author: Gaurav Bhutani
Affiliation: Indian Institute of Technology (IIT) Mandi
Repository: https://github.com/gbhutani/hpsc_2025/blob/main/codes/hpc/python_script_email_hpc_job_status.py

modified by: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
-
- This file is part of the HHG-SaDAS package, modified during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

import time
import smtplib
import argparse
import os, subprocess
from email.mime.text import MIMEText

# Default values (modifiable)
DEFAULT_HPC_PATH = "/home/s23092/colloc_pt"
DEFAULT_JOBID = "568181"


parser = argparse.ArgumentParser(description="Monitor HPC job and send email notifications.")
parser.add_argument("--jobid", type=str, default=DEFAULT_JOBID, help="HPC job ID to monitor.")
parser.add_argument("--hpc_path", type=str, default=DEFAULT_HPC_PATH, help="Path to the HPC job directory.")
args = parser.parse_args()

hpc_path = args.hpc_path        # Uses command-line argument if provided, otherwise keeps the default
jobid = args.jobid              # Uses command-line argument if provided, otherwise keeps the default

hpc_userid = "s23092"
print(f"Monitoring job {jobid} in {hpc_path}")



completed_file_name = 'out_N200.o'
local_copy_path = r"/Assistant/AlertSystem_Email_HPC/data_from_cluster"
completed_local_file_path = os.path.join(local_copy_path, completed_file_name)

completed_flag = False
cmd_copy_completed_file = f'scp {hpc_userid}@10.8.1.19:{hpc_path}/{completed_file_name} {local_copy_path}'
def send_email():
    myEmail = 'your_email_id@gmail.com'
    smtp_server = 'smtp.gmail.com'
    smtp_port = 587
    server = smtplib.SMTP(smtp_server, smtp_port)
    server.starttls()

    with open(r'path\to\password\random_letters.txt', 'r') as file:
        password = file.read().strip()
    server.login(myEmail, password)
    print('logged in Successfully!  Sending mail...')

    with open(completed_local_file_path, 'r') as file:
        contents = file.read()

    msg = MIMEText(contents)
    msg['From'] = myEmail
    msg['To'] = myEmail
    msg['Subject'] = jobid + ' completed'
    server.send_message(msg)
    server.quit()



attempt_no = 0
while True:

    if not completed_flag:
        print(f'[Attempt-{attempt_no}]: looking for: {completed_file_name}')
        subprocess.call(cmd_copy_completed_file, shell=True, stderr=subprocess.DEVNULL)
        if os.path.exists(completed_local_file_path):
            print("File found and copied locally!")
            send_email()
            completed_flag = True
        attempt_no += 1

    else: break
    time.sleep(5)
