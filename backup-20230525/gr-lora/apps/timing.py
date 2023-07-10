#!/usr/bin/env python3
import os
import time

count = 0;
while True:
	os.system("/usr/bin/python3 /home/zkw/gr-lora/apps/lora_realtime_clip_signal.py")
	count = count + 1
	print("Done\n")
	if count == 1:
		break
	time.sleep(3*60)
