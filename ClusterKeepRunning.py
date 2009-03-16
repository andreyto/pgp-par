import os
import sys
import time
from ClusterUtils import *

QSUBCommands = """
qsub j3_inf_sol_01_itms.0.sh
qsub j3_inf_sol_03_itms.0.sh
qsub j3_inf_sol_05_itms.0.sh
qsub j3_inf_sol_07_itms.0.sh
qsub j3_inf_sol_09_itms.0.sh
qsub j3_inf_sol_11_itms.0.sh
qsub j3_inf_sol_1_01_itms.0.sh
qsub j3_inf_sol_1_03_itms.0.sh
qsub j3_inf_sol_1_05_itms.0.sh
qsub j3_inf_sol_1_07_itms.0.sh
qsub j3_inf_sol_1_09_itms.0.sh
qsub j3_inf_sol_1_11_itms.0.sh
qsub j3_ni_III_01_itms.0.sh
qsub j3_ni_III_03_itms.0.sh
qsub j3_ni_III_05_itms.0.sh
qsub j3_ni_III_07_itms.0.sh
qsub j3_ni_III_09_itms.0.sh
qsub j3_ni_III_11_itms.0.sh
qsub j3_ni_II_01_itms.0.sh
qsub j3_ni_II_03_itms.0.sh
qsub j3_ni_II_05_itms.0.sh
qsub j3_ni_II_07_itms.0.sh
qsub j3_ni_II_09_itms.0.sh
qsub j3_ni_II_11_itms.0.sh
qsub j3_ni_I_01_080630195258_itms.0.sh
qsub j3_ni_I_03_itms.0.sh
qsub j3_ni_I_05_itms.0.sh
qsub j3_ni_I_07_itms.0.sh
qsub j3_ni_I_09_itms.0.sh
qsub j3_ni_I_11_itms.0.sh
qsub j3is_01_itms.0.sh
qsub j3is_02_itms.0.sh
qsub j3is_04_itms.0.sh
qsub j3is_06_itms.0.sh
qsub j3is_08_itms.0.sh
qsub j3is_10_itms.0.sh
qsub j3is_12_itms.0.sh
qsub j3ni_02_itms.0.sh
qsub j3ni_04_itms.0.sh
qsub j3ni_06_itms.0.sh
qsub j3ni_08_itms.0.sh
qsub j3ni_10_itms.0.sh
qsub j3ni_12_itms.0.sh
qsub j3ni_sol_3_02_itms.0.sh
qsub j3ni_sol_3_04_itms.0.sh
qsub j3ni_sol_3_06_itms.0.sh
qsub j3ni_sol_3_08_itms.0.sh
qsub j3ni_sol_3_10_itms.0.sh
qsub j3ni_sol_3_12_itms.0.sh
qsub j3p1_02_itms.0.sh
qsub j3p1_04_itms.0.sh
qsub j3p1_06_itms.0.sh
qsub j3p1_08_itms.0.sh
qsub j3p1_10_itms.0.sh
qsub j3p1_12_itms.0.sh
qsub j3p1_1_02_itms.0.sh
qsub j3p1_1_04_itms.0.sh
qsub j3p1_1_06_itms.0.sh
qsub j3p1_1_08_itms.0.sh
qsub j3p1_1_10_itms.0.sh
qsub j3p1_1_12_itms.0.sh
qsub j6_inf_sol_02_itms.0.sh
qsub j6_inf_sol_04_itms.0.sh
qsub j6_inf_sol_06_itms.0.sh
qsub j6_inf_sol_08_itms.0.sh
qsub j6_inf_sol_10_itms.0.sh
qsub j6_inf_sol_12_itms.0.sh
qsub j6ni_1_02_itms.0.sh
qsub j6ni_1_04_itms.0.sh
qsub j6ni_1_06_itms.0.sh
qsub j6ni_1_08_itms.0.sh
qsub j6ni_1_10_itms.0.sh
qsub j6ni_1_12_itms.0.sh
qsub j6ni_2_02_itms.0.sh
qsub j6ni_2_04_itms.0.sh
qsub j6ni_2_06_itms.0.sh
qsub j6ni_2_08_itms.0.sh
qsub j6ni_2_10_itms.0.sh
qsub j6ni_2_12_itms.0.sh
qsub j6pi_02_itms.0.sh
qsub j6pi_04_itms.0.sh
qsub j6pi_06_itms.0.sh
qsub j6pi_08_itms.0.sh
qsub j6pi_10_itms.0.sh
qsub j6pi_12_itms.0.sh
qsub jPARC_3II_rep3_0208_02_itms.0.sh
qsub jPARC_3II_rep3_0208_04_itms.0.sh
qsub jPARC_3II_rep3_0208_06_itms.0.sh
qsub jPARC_3II_rep3_0208_08_itms.0.sh
qsub jPARC_3II_rep3_0208_10_itms.0.sh
qsub jPARC_3II_rep3_0208_12_itms.0.sh
qsub jPARC_6II_rep1_0707_02_itms.0.sh
qsub jPARC_6II_rep1_0707_04_itms.0.sh
qsub jPARC_6II_rep1_0707_06_itms.0.sh
qsub jPARC_6II_rep1_0707_08_itms.0.sh
qsub jPARC_6II_rep1_0707_10_itms.0.sh
qsub jPARC_6II_rep1_0707_12_itms.0.sh
qsub jPARC_6II_rep4_0208_02_itms.0.sh
qsub jPARC_6II_rep4_0208_04_itms.0.sh
qsub jPARC_6II_rep4_0208_06_itms.0.sh
qsub jPARC_6II_rep4_0208_08_itms.0.sh
qsub jPARC_6II_rep4_0208_10_itms.0.sh
qsub jPARC_6II_rep4_0208_12_itms.0.sh
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
    