#!/usr/bin/env python3
import os
import time
import numpy

matlabFlagName1 = "/mnt/hgfs/share/notDone.txt"
matlabFlagName2 = "/mnt/hgfs/share/Done.txt"
grloraFlagscript = "touch /mnt/hgfs/share/grloraDone.txt"
cliper_script = "/usr/bin/python3 /home/zkw/gr-lora/apps/lora_realtime_clip_signal.py --num "
while True:
    # 等待matlab处理
    print("wait matlab...")
    while not os.path.exists(matlabFlagName1) and not os.path.exists(matlabFlagName2):
        time.sleep(1)
    # matlab处理完成，读取数据并删除标志位
    print("matlab Done!")
    if os.path.exists(matlabFlagName1):
        transArray = numpy.loadtxt(matlabFlagName1, dtype='int', delimiter=',')
        os.remove(matlabFlagName1)
    elif os.path.exists(matlabFlagName2):
        transArray = numpy.loadtxt(matlabFlagName2, dtype='int', delimiter=',')
        os.remove(matlabFlagName2)
    # 采集信号
    transArray = transArray.tolist();
    if isinstance(transArray, list):
        cliper_num = len(transArray)
    else:
        cliper_num = 1
    
    cliper_script_all = cliper_script + str(cliper_num)
    os.system(cliper_script_all)
    # 采集完成写入标志位
    os.system(grloraFlagscript)
