import os
import sys
import time
from ClusterUtils import *

QSUBCommands = """
qsub jA415A001.0.sh
qsub jA415A003.0.sh
qsub jA415A005.0.sh
qsub jA415A007.0.sh
qsub jA415A009.0.sh
qsub jA415A011.0.sh
qsub jA415A013.0.sh
qsub jA415A015.0.sh
qsub jA415A017.0.sh
qsub jA415A019.0.sh
qsub jA415A021.0.sh
qsub jA415A023.0.sh
qsub jA415A025.0.sh
qsub jA415A027.0.sh
qsub jA415A029.0.sh
qsub jA415A031.0.sh
qsub jA415A033.0.sh
qsub jA415A035.0.sh
qsub jA415A037.0.sh
qsub jA415A039.0.sh
qsub jA415A041.0.sh
qsub jA415A043.0.sh
qsub jA415A045.0.sh
qsub jA415A047.0.sh
qsub jA415A049.0.sh
qsub jA415A051.0.sh
qsub jA415A053.0.sh
qsub jA415A055.0.sh
qsub jA415A057.0.sh
qsub jA415A059.0.sh
qsub jA415A061.0.sh
qsub jA415A063.0.sh
qsub jA415A065.0.sh
qsub jA415A067.0.sh
qsub jA415A069.0.sh
qsub jA415A071.0.sh
qsub jA415A073.0.sh
qsub jA415A075.0.sh
qsub jA415A077.0.sh
qsub jA415A079.0.sh
qsub jA415A082.0.sh
qsub jA415A084.0.sh
qsub jA415A086.0.sh
qsub jA415A088.0.sh
qsub jA415A090.0.sh
qsub jA415A092.0.sh
qsub jA415A094.0.sh
qsub jA415A096.0.sh
qsub jA415B002.0.sh
qsub jA415B004.0.sh
qsub jA415B006.0.sh
qsub jA415B008.0.sh
qsub jA415B010.0.sh
qsub jA415B013.0.sh
qsub jA415B015.0.sh
qsub jA415B017.0.sh
qsub jA415B019.0.sh
qsub jA415B021.0.sh
qsub jA415B024.0.sh
qsub jA415B026.0.sh
qsub jA415B028.0.sh
qsub jA415B030.0.sh
qsub jA415B032.0.sh
qsub jA415B035.0.sh
qsub jA415B037.0.sh
qsub jA415B039.0.sh
qsub jA415B041.0.sh
qsub jA415B043.0.sh
qsub jA415B045.0.sh
qsub jA415B048.0.sh
qsub jA415B050.0.sh
qsub jA415B052.0.sh
qsub jA415B054.0.sh
qsub jA415B057.0.sh
qsub jA415B059.0.sh
qsub jA415B061.0.sh
qsub jA415B063.0.sh
qsub jA415B065.0.sh
qsub jA415B067.0.sh
qsub jA415B069.0.sh
qsub jA415B071.0.sh
qsub jA415B073.0.sh
qsub jA415B075.0.sh
qsub jA415B077.0.sh
qsub jA415B079.0.sh
qsub jA415B081.0.sh
qsub jA415B083.0.sh
qsub jA415B085.0.sh
qsub jA415B087.0.sh
qsub jA415B090.0.sh
qsub jA415B092.0.sh
qsub jA415B094.0.sh
qsub jA415B096.0.sh
qsub jA415C002.0.sh
qsub jA415C004.0.sh
qsub jA415C007.0.sh
qsub jA415C009.0.sh
qsub jA415C011.0.sh
qsub jA415C013.0.sh
qsub jA415C015.0.sh
qsub jA415C017.0.sh
qsub jA415C019.0.sh
qsub jA415C022.0.sh
qsub jA415C024.0.sh
qsub jA415C026.0.sh
qsub jA415C028.0.sh
qsub jA415C030.0.sh
qsub jA415C032.0.sh
qsub jA415C034.0.sh
qsub jA415C036.0.sh
qsub jA415C038.0.sh
qsub jA415C040.0.sh
qsub jA415C042.0.sh
qsub jA415C044.0.sh
qsub jA415C046.0.sh
qsub jA415C048.0.sh
qsub jA415C050.0.sh
qsub jA415C052.0.sh
qsub jA415C054.0.sh
qsub jA415C056.0.sh
qsub jA415C059.0.sh
qsub jA415C061.0.sh
qsub jA415C063.0.sh
qsub jA415C065.0.sh
qsub jA417A002.0.sh
qsub jA417A004.0.sh
qsub jA417A006.0.sh
qsub jA417A008.0.sh
qsub jA417A010.0.sh
qsub jA417A012.0.sh
qsub jA417A014.0.sh
qsub jA417A018.0.sh
qsub jA417A020.0.sh
qsub jA417A022.0.sh
qsub jA417A024.0.sh
qsub jA417A026.0.sh
qsub jA417A028.0.sh
qsub jA417A030.0.sh
qsub jA417A032.0.sh
qsub jA417A034.0.sh
qsub jA417A036.0.sh
qsub jA417A038.0.sh
qsub jA417A040.0.sh
qsub jA417A042.0.sh
qsub jA417A044.0.sh
qsub jA417A046.0.sh
qsub jA417A048.0.sh
qsub jA417A050.0.sh
qsub jA417A052.0.sh
qsub jA417A054.0.sh
qsub jA417A056.0.sh
qsub jA417A058.0.sh
qsub jA417A060.0.sh
qsub jA417A062.0.sh
qsub jA417A064.0.sh
qsub jA417A066.0.sh
qsub jA417A068.0.sh
qsub jA417A070.0.sh
qsub jA417A072.0.sh
qsub jA417A074.0.sh
qsub jA417A076.0.sh
qsub jA417A078.0.sh
qsub jA417A081.0.sh
qsub jA417A083.0.sh
qsub jA417A085.0.sh
qsub jA417A087.0.sh
qsub jA417A089.0.sh
qsub jA417A091.0.sh
qsub jA417A093.0.sh
qsub jA417A095.0.sh
qsub jA420A001.0.sh
qsub jA420A004.0.sh
qsub jA420A009.0.sh
qsub jA420A016.0.sh
qsub jA420A021.0.sh
qsub jA420A024.0.sh
qsub jA420A027.0.sh
qsub jA420A029.0.sh
qsub jA420A041.0.sh
qsub jA420A055.0.sh
qsub jA420A059.0.sh
qsub jA420A062.0.sh
qsub jA420A068.0.sh
qsub jA420A072.0.sh
qsub jA420A074.0.sh
qsub jA420A080.0.sh
qsub jA420A087.0.sh
qsub jA420A093.0.sh
qsub jA424A001.0.sh
qsub jA424A003.0.sh
qsub jA424A005.0.sh
qsub jA424A007.0.sh
qsub jA424A009.0.sh
qsub jA424A011.0.sh
qsub jA424A013.0.sh
qsub jA424A015.0.sh
qsub jA424A017.0.sh
qsub jA424A019.0.sh
qsub jA424A021.0.sh
qsub jA424A023.0.sh
qsub jA424A025.0.sh
qsub jA424A027.0.sh
qsub jA424A029.0.sh
qsub jA424A031.0.sh
qsub jA424A033.0.sh
qsub jA424A036.0.sh
qsub jA424A038.0.sh
qsub jA424A040.0.sh
qsub jA424A042.0.sh
qsub jA424A045.0.sh
qsub jA424A047.0.sh
qsub jA424A049.0.sh
qsub jA424A051.0.sh
qsub jA424A053.0.sh
qsub jA424A055.0.sh
qsub jA424A058.0.sh
qsub jA424A061.0.sh
qsub jA424A065.0.sh
qsub jA424A067.0.sh
qsub jA424A069.0.sh
qsub jA424A072.0.sh
qsub jA424A074.0.sh
qsub jA424A076.0.sh
qsub jA424A078.0.sh
qsub jA424A080.0.sh
qsub jA424A083.0.sh
qsub jA424A085.0.sh
qsub jA424A087.0.sh
qsub jA424A089.0.sh
qsub jA424A092.0.sh
qsub jA424A094.0.sh
qsub jA424A096.0.sh
qsub jA425A004.0.sh
qsub jA425A006.0.sh
qsub jA425A009.0.sh
qsub jA425A011.0.sh
qsub jA425A013.0.sh
qsub jA425A015.0.sh
qsub jA425A017.0.sh
qsub jA425A019.0.sh
qsub jA425A021.0.sh
qsub jA425A023.0.sh
qsub jA425A025.0.sh
qsub jA425A027.0.sh
qsub jA425A029.0.sh
qsub jA425A031.0.sh
qsub jA425A033.0.sh
qsub jA425A035.0.sh
qsub jA425A037.0.sh
qsub jA425A039.0.sh
qsub jA425A041.0.sh
qsub jA425A043.0.sh
qsub jA425A045.0.sh
qsub jA425A047.0.sh
qsub jA425B002.0.sh
qsub jA425B004.0.sh
qsub jA425B006.0.sh
qsub jA425B008.0.sh
qsub jA425B010.0.sh
qsub jA425B012.0.sh
qsub jA425B014.0.sh
qsub jA425B016.0.sh
qsub jA425B018.0.sh
qsub jA425B020.0.sh
qsub jA425B022.0.sh
qsub jA425B024.0.sh
qsub jA425B026.0.sh
qsub jA425B028.0.sh
qsub jA425B030.0.sh
qsub jA425B032.0.sh
qsub jA425B034.0.sh
qsub jA425B036.0.sh
qsub jA425B038.0.sh
qsub jA425B040.0.sh
qsub jA425B042.0.sh
qsub jA425B044.0.sh
qsub jA425B046.0.sh
qsub jA425B048.0.sh
qsub jA443A003.0.sh
qsub jA443A005.0.sh
qsub jA443A007.0.sh
qsub jA443A010.0.sh
qsub jA443A012.0.sh
qsub jA443A015.0.sh
qsub jA443A018.0.sh
qsub jA443A020.0.sh
qsub jA443A022.0.sh
qsub jA443A024.0.sh
qsub jA443A029.0.sh
qsub jA443A031.0.sh
qsub jA443A033.0.sh
qsub jA443A035.0.sh
qsub jA443A037.0.sh
qsub jA443A039.0.sh
qsub jA443A041.0.sh
qsub jA443A043.0.sh
qsub jA443A045.0.sh
qsub jA443A047.0.sh
qsub jA443A049.0.sh
qsub jA443A051.0.sh
qsub jA443A053.0.sh
qsub jA443A056.0.sh
qsub jA443A058.0.sh
qsub jA443A060.0.sh
qsub jA443A062.0.sh
qsub jA443A064.0.sh
qsub jA443A066.0.sh
qsub jA443A068.0.sh
qsub jA443A070.0.sh
qsub jA443A072.0.sh
qsub jA443A074.0.sh
qsub jA443A076.0.sh
qsub jA443A078.0.sh
qsub jA443A081.0.sh
qsub jA443A083.0.sh
qsub jA443A085.0.sh
qsub jA443A087.0.sh
qsub jA443A089.0.sh
qsub jA443A092.0.sh
qsub jA443A094.0.sh
qsub jA443A096.0.sh
qsub jA443B002.0.sh
qsub jA443B006.0.sh
qsub jA443B008.0.sh
qsub jA443B016.0.sh
qsub jA443B019.0.sh
qsub jA443B024.0.sh
qsub jA443B030.0.sh
qsub jA443B032.0.sh
qsub jA443B034.0.sh
qsub jA443B038.0.sh
qsub jA443B043.0.sh
qsub jA443B045.0.sh
qsub jA443B049.0.sh
qsub jA443B056.0.sh
qsub jA443B063.0.sh
qsub jA443B069.0.sh
qsub jA443B073.0.sh
qsub jA444A001.0.sh
qsub jA444A003.0.sh
qsub jA444A005.0.sh
qsub jA444A007.0.sh
qsub jA444A009.0.sh
qsub jA444A011.0.sh
qsub jA444A013.0.sh
qsub jA444A015.0.sh
qsub jA444A018.0.sh
qsub jA444A020.0.sh
qsub jA444A022.0.sh
qsub jA444A024.0.sh
qsub jA444A026.0.sh
qsub jA444A028.0.sh
qsub jA444A030.0.sh
qsub jA444A032.0.sh
qsub jA444A034.0.sh
qsub jA444A036.0.sh
qsub jA444A038.0.sh
qsub jA444A040.0.sh
qsub jA444A042.0.sh
qsub jA444A044.0.sh
qsub jA444A046.0.sh
qsub jA444A049.0.sh
qsub jA444A051.0.sh
qsub jA444A053.0.sh
qsub jA444A055.0.sh
qsub jA444A059.0.sh
qsub jA444A061.0.sh
qsub jA444A063.0.sh
qsub jA444A066.0.sh
qsub jA444A068.0.sh
qsub jA444A070.0.sh
qsub jA444A072.0.sh
qsub jA444A074.0.sh
qsub jA444A076.0.sh
qsub jA444A078.0.sh
qsub jA444A080.0.sh
qsub jA444A082.0.sh
qsub jA444A084.0.sh
qsub jA444A086.0.sh
qsub jA444A088.0.sh
qsub jA444A091.0.sh
qsub jA444A093.0.sh
qsub jA444B001.0.sh
qsub jA444B003.0.sh
qsub jA444B010.0.sh
qsub jA444B014.0.sh
qsub jA444B020.0.sh
qsub jA444B023.0.sh
qsub jA444B026.0.sh
qsub jA444B030.0.sh
qsub jA444B032.0.sh
qsub jA444B035.0.sh
qsub jA444B039.0.sh
qsub jA444B044.0.sh
qsub jA444B046.0.sh
qsub jA444B054.0.sh
qsub jA444B057.0.sh
qsub jA444B060.0.sh
qsub jA444B063.0.sh
qsub jA444B067.0.sh
qsub jA444B069.0.sh
qsub jA444B072.0.sh
qsub jA444B074.0.sh
qsub jA444B083.0.sh
qsub jA444B085.0.sh
qsub jA444B088.0.sh
qsub jA444B090.0.sh
qsub jA444B093.0.sh
qsub jA444C002.0.sh
qsub jA444C018.0.sh
qsub jA444C031.0.sh
qsub jA444C039.0.sh
qsub jA444C058.0.sh
qsub jA470A001.0.sh
qsub jA470A003.0.sh
qsub jA470A005.0.sh
qsub jA470A007.0.sh
qsub jA470A010.0.sh
qsub jA470A012.0.sh
qsub jA470A014.0.sh
qsub jA470A017.0.sh
qsub jA470A019.0.sh
qsub jA470A021.0.sh
qsub jA470A023.0.sh
qsub jA470A026.0.sh
qsub jA470A028.0.sh
qsub jA470A030.0.sh
qsub jA470A032.0.sh
qsub jA470A034.0.sh
qsub jA470A036.0.sh
qsub jA470A038.0.sh
qsub jA470A042.0.sh
qsub jA470A044.0.sh
qsub jA470A046.0.sh
qsub jA470A048.0.sh
qsub jA470A051.0.sh
qsub jA470A056.0.sh
qsub jA474A005.0.sh
qsub jA474A007.0.sh
qsub jA474A011.0.sh
qsub jA474A013.0.sh
qsub jA474A016.0.sh
qsub jA474A018.0.sh
qsub jA474A024.0.sh
qsub jA474A026.0.sh
qsub jA474A028.0.sh
qsub jA474A031.0.sh
qsub jA474A034.0.sh
qsub jA474A036.0.sh
qsub jA474A042.0.sh
qsub jA474A044.0.sh
qsub jA474A047.0.sh
qsub jA474A051.0.sh
qsub jA474A054.0.sh
qsub jA474A057.0.sh
qsub jA474A061.0.sh
qsub jA474A065.0.sh
qsub jA474A071.0.sh
qsub jA474A076.0.sh
qsub jA474A081.0.sh
qsub jA474A091.0.sh
qsub jA476A002.0.sh
qsub jA476A004.0.sh
qsub jA476A006.0.sh
qsub jA476A008.0.sh
qsub jA476A010.0.sh
qsub jA476A012.0.sh
qsub jA476A014.0.sh
qsub jA476A016.0.sh
qsub jA476A018.0.sh
qsub jA476A020.0.sh
qsub jA476A022.0.sh
qsub jA476A024.0.sh
qsub jA476A026.0.sh
qsub jA476A028.0.sh
qsub jA476A030.0.sh
qsub jA476A032.0.sh
qsub jA476A034.0.sh
qsub jA476A036.0.sh
qsub jA476A038.0.sh
qsub jA476A040.0.sh
qsub jA476A042.0.sh
qsub jA476A044.0.sh
qsub jA476A047.0.sh
qsub jA476A049.0.sh
qsub jA476A051.0.sh
qsub jA476A053.0.sh
qsub jA476A055.0.sh
qsub jA476A057.0.sh
qsub jA476A059.0.sh
qsub jA476A062.0.sh
qsub jA476A064.0.sh
qsub jA478A002.0.sh
qsub jA478A004.0.sh
qsub jA478A006.0.sh
qsub jA478A008.0.sh
qsub jA478A010.0.sh
qsub jA478A013.0.sh
qsub jA478A015.0.sh
qsub jA478A017.0.sh
qsub jA478A019.0.sh
qsub jA478A021.0.sh
qsub jA478A023.0.sh
qsub jA478A025.0.sh
qsub jA478A029.0.sh
qsub jA478A031.0.sh
qsub jA478A033.0.sh
qsub jA478A035.0.sh
qsub jA478A037.0.sh
qsub jA478A039.0.sh
qsub jA478A041.0.sh
qsub jA478A044.0.sh
qsub jA478A046.0.sh
qsub jA478A048.0.sh
qsub jA478A050.0.sh
qsub jA478A052.0.sh
qsub jA478A054.0.sh
qsub jA478A056.0.sh
qsub jA478A058.0.sh
qsub jA478A060.0.sh
qsub jA478A062.0.sh
qsub jA478A064.0.sh
qsub jA478A066.0.sh
qsub jA478A068.0.sh
qsub jA478A070.0.sh
qsub jA478A072.0.sh
qsub jA478A074.0.sh
qsub jA478A076.0.sh
qsub jA478A078.0.sh
qsub jA478A081.0.sh
qsub jA478A083.0.sh
qsub jA478A086.0.sh
qsub jA478A089.0.sh
qsub jA478A091.0.sh
qsub jA478A093.0.sh
qsub jA478A095.0.sh
qsub jA479A002.0.sh
qsub jA479A008.0.sh
qsub jA479A013.0.sh
qsub jA479A020.0.sh
qsub jA479A024.0.sh
qsub jA479A029.0.sh
qsub jA479A033.0.sh
qsub jA479A035.0.sh
qsub jA479A038.0.sh
qsub jA479A043.0.sh
qsub jA479A048.0.sh
qsub jA479A050.0.sh
qsub jA479A059.0.sh
qsub jA479A061.0.sh
qsub jA479A073.0.sh
qsub jA479A085.0.sh
qsub jA479A089.0.sh
qsub jA479A094.0.sh
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
    