import os
import sys
import time
from ClusterUtils import *

QSUBCommands = """
qsub jA07_4297_LepOGE_F1_c.0.sh
qsub jA07_4299_LepOGE_F3_c.0.sh
qsub jA07_4301_LepOGE_F5_c.0.sh
qsub jA07_4304_LepOGE_F7_c.0.sh
qsub jA07_4306_LepOGE_F9_c.0.sh
qsub jA07_6292_c.0.sh
qsub jA07_6294_c.0.sh
qsub jA07_6296_c.0.sh
qsub jA07_6298_c.0.sh
qsub jA07_6300_c.0.sh
qsub jA07_6389_OGE2_3_7_11_c.0.sh
qsub jA07_6391_OGE2_3_7_13_c.0.sh
qsub jA07_6393_OGE2_3_7_15_c.0.sh
qsub jA07_6395_OGE2_3_7_17_c.0.sh
qsub jA07_6397_OGE2_3_7_19_c.0.sh
qsub jA07_7072_OGE2_3_7_21_c.0.sh
qsub jA07_7074_OGE2_3_7_23_c.0.sh
qsub jA07_7076_OGE2_3_10_1_c.0.sh
qsub jA07_7078_OGE2_3_10_3_c.0.sh
qsub jA07_7080_OGE2_3_10_5_c.0.sh
qsub jA07_7082_OGE2_3_10_7_c.0.sh
qsub jA07_7119_OGE2_3_10_9_c.0.sh
qsub jA07_7121_OGE2_3_10_11_c.0.sh
qsub jA07_7123_OGE2_3_10_13_c.0.sh
qsub jA07_7125_OGE2_3_10_15_c.0.sh
qsub jA07_7243_OGE2_3_10_17_c.0.sh
qsub jA07_7245_OGE2_3_10_19_c.0.sh
qsub jA07_7247_OGE2_3_10_21_c.0.sh
qsub jA07_7249_OGE2_3_10_23_c.0.sh
"""

JobList = []
for JobChunk in QSUBCommands.split("\n"):
    String = JobChunk.strip()
    if String:
        JobList.append(String)

while 1:
    JobCount = GetRunningJobCount()
    while JobCount <= MAXIMUM_JOBS:
        JobCommand = JobList.pop(0)
        print JobCommand
        os.system(JobCommand)
        JobCount += 1
        if len(JobList) <= 0:
            break
    print "(Running %s jobs; %s awaiting submission...)"%(JobCount, len(JobList))
    time.sleep(59)
    continue
    
