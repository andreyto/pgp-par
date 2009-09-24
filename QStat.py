"""
Wrapper for qstat that reveals the ACTUAL NAME of the damn job script.
"""
import os
import re

class InfoBag:
    def __init__(self):
        self.ScriptName = None
    def __cmp__(self, Other):
        return cmp(self.ScriptName, getattr(Other, "ScriptName", None))

# Assume that the environment variable USER is set properly.
USER_NAME = os.environ["USER"]
Command = "qstat -u %s"%USER_NAME
Pipe = os.popen(Command, "r")
Jobs = []
# Read the job IDs and submission times from qstat:
for FileLine in Pipe.xreadlines():
    Bits = FileLine.split()
    try:
        JobID = int(Bits[0])
        JobInfo = InfoBag()
        JobInfo.ID = JobID
        JobInfo.Status = Bits[4]
        JobInfo.SubmitTime = "%s %s"%(Bits[5], Bits[6])
    except:
        continue
    Jobs.append(JobInfo)
Pipe.close()

Jobs.sort()
# Run qstat -j to get the script name for each job of interest:
for Job in Jobs:
    Command = "qstat -j %s"%Job.ID
    Pipe = os.popen(Command)
    for FileLine in Pipe.xreadlines():
        if FileLine[:9] == "job_name:":
            Job.ScriptName = FileLine[9:].strip()
            break
    Pipe.close()

Jobs.sort()
for Job in Jobs:
    print "%s\t%s\t%s\t%s\t"%(Job.ID, Job.ScriptName, Job.Status, Job.SubmitTime)
print
print "%s jobs running."%len(Jobs)