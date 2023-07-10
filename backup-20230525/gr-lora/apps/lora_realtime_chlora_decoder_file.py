#!/usr/bin/env python3
# Gets a USRP capture trace from my research page and decodes it using gr-lora.
# Author: Pieter Robyns

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
from lora.loraconfig import LoRaConfig
from time import sleep
import lora
import argparse
import os
import requests
import json
import osmosdr
import ctypes
import sys

class lora_receive_file_nogui(gr.top_block):
    def __init__(self, pkg_num_all, data_file_path):
        gr.top_block.__init__(self, "Lora Receive File, No GUI")

        ##################################################
        # Variables
        ##################################################
        self.sf = sf = 10
        self.bw = bw = 125000
        self.target_freq = target_freq = 433.5e6
        # self.target_freq = target_freq = 470.0e6
        self.symbols_per_sec = symbols_per_sec = float(bw) / (2**sf)
        self.sample_rate = sample_rate = 2e6
        # self.sample_rate = sample_rate = 1e6
        self.decimation = decimation = 1
        self.capture_freq = capture_freq = 433.5e6
        # self.capture_freq = capture_freq = 470.0e6
        self.bitrate = bitrate = sf * (1 / (2**sf / float(bw)))
        self.decimation = 1
        cr = "4/5";

        ##################################################
        # Blocks
        ##################################################
        # self.osmosdr_source_0 = osmosdr.source(
        #     args="numchan=" + str(1) + " " + ''
        # )
        # self.osmosdr_source_0.set_sample_rate(sample_rate)
        # self.osmosdr_source_0.set_center_freq(capture_freq, 0)
        # self.osmosdr_source_0.set_freq_corr(0, 0)
        # self.osmosdr_source_0.set_gain(10, 0)
        # self.osmosdr_source_0.set_if_gain(20, 0)
        # self.osmosdr_source_0.set_bb_gain(20, 0)
        # self.osmosdr_source_0.set_antenna('', 0)
        # self.osmosdr_source_0.set_bandwidth(0, 0)
        self.lora_message_socket_sink_0 = lora.message_socket_sink('127.0.0.1', 40868, 0)
        lora_config = LoRaConfig(capture_freq, sf, cr, bw, prlen=8, crc=True, implicit=False)
        self.lora_receiver = lora.chlora_decoder(sample_rate, capture_freq, ([lora_config.freq]), lora_config.bw, lora_config.sf, lora_config.implicit, lora_config.cr_num, lora_config.crc, pkg_num_all, reduced_rate=False, decimation=self.decimation)
        self.blocks_throttle = blocks.throttle(gr.sizeof_gr_complex, sample_rate, True)
        self.blocks_file_source = blocks.file_source(gr.sizeof_gr_complex, data_file_path, False)

        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.lora_receiver, 'frames'), (self.lora_message_socket_sink_0, 'in'))
        # self.connect((self.osmosdr_source_0, 0), (self.lora_receiver, 0))
        self.connect((self.blocks_file_source, 0), (self.blocks_throttle, 0))
        self.connect((self.blocks_throttle, 0), (self.lora_receiver, 0))
        # self.connect((self.blocks_file_source, 0), (self.lora_receiver, 0))

def download_file(source, destination):
    print("[+] Downloading %s -> %s" % (source, destination)),
    try:
        response = requests.get(source, stream=True)

        with open(destination, 'wb') as f:
            for chunk in response.iter_content(1000*1000):
                f.write(chunk)
                print("."),
            print(".")
    except Exception as e:
        print("[-] " + str(e))
        exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Gets a USRP capture trace from my research page and decodes it using gr-lora.")
    parser.add_argument('--num', type=int, default=128, help='The num of pkg')
    # parser.add_argument('--file', type=str, default="example-trace", help='File name of the example trace.')
    args = parser.parse_args()
    # get the file name
    # list to store files
    dir_path = '/mnt/hgfs/share/samples/'
    res = []
    for path in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            res.append(path);
    # append file suffix
    data_file_path = os.path.join('/mnt/hgfs/share/samples/', res[0])
    # meta_file_path = os.path.join('./', args.file + '.sigmf-meta')
    # check if exist 
    # print(data_file_path)
    if not os.path.exists(data_file_path):
        raise SystemExit('file not exist!')
    
    tb = lora_receive_file_nogui(args.num, data_file_path)
    tb.start()
    tb.wait()
    
    
    
    
    
    
    
    
