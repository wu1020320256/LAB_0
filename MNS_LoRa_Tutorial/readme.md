# 整体介绍
--------------------
&ensp;&ensp;该教程主要负责Mobinets LoRa组信号处理部分的内容，旨在教会大家如何使用matlab、GNURadio配合rtl-sdr、USRP硬件进行信号的采集和处理，达到实验所需要求

## 背景
--------------------
&ensp;&ensp;当实验需要对LoRa信号进行调制解调的分析时，我们可以通过python、matlab等其他编程语言，通过公式进行信号仿真，来非常便捷的原理验证，如下：  

    cmx = 1+1*1i;                                         % 调整相位
    pre_dir = 2*pi;                                       % 2pi参数
    f0 = lora_set.bw/2;                                   % 设置理想upchirp和downchirp的初始频率
    d_symbols_per_second = lora_set.bw / lora_set.fft_x;  % 计算每秒symbol数目 
    T = -0.5 * lora_set.bw * d_symbols_per_second; 
    d_samples_per_second = 1000000;                       % sdr-rtl的采样率
    d_dt = 1/d_samples_per_second;                        % 采样点间间隔的时间
    t = d_dt*(0:1:lora_set.dine-1);

    % 计算理想downchirp和upchirp存入d_downchirp和d_upchirp数组中（复数形式）
    d_downchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t)) + sin(pre_dir .* t .* (f0 + T * t))*1i);
    d_upchirp = cmx * (cos(pre_dir .* t .* (f0 + T * t) * -1) + sin(pre_dir .* t .* (f0 + T * t) * -1)*1i);

&ensp;&ensp;将上述代码在matlab中直接运行，即可生成一个理想的ideachirp存入在d_downchirp和d_upchirp变量中，用于后续处理。  

&ensp;&ensp;当在matlab或者python进行理论验证后，我们需要实际生成这样一个信号，并且采集在实际环境中的信号进行解调和处理，与理论相印证。此时我们就需要使用实验室的LoRa testbed或者USRP发送一个实际的LoRa信号，然后通过硬件+软件的方式采集testbed发出的信号，并再通过软件的方式处理信号，完成实验步骤。

    发送信号                   采集信号                   信号处理
    testbed  ------------->   rtl-sdr  ----------->  gr-lora/Matlab

&ensp;&ensp;其中需要大家对信号处理流程，通信原理和信号与系统有一个基本的了解，比如需要清楚一个信号发出后是以模拟电磁波的方式在空中散射，但在经过rtl-sdr或者其他硬件设备的采样后变成计算机能够读懂的数字信号：  

                模拟信号                采样/数字信号
    testbed  ------------->   rtl-sdr  ----------->  gr-lora/Matlab

&ensp;&ensp;对于LoRa基础知识可以看[Known and Unknown Facts of LoRa: Experiences from a Large-scale Measurement Study](https://dl.acm.org/doi/10.1145/3293534)，通信原理、信号与系统和数字信号处理可以通过b站和书籍的方式学习，这里不过多介绍了（自己学的也不好）

## 组件介绍
--------------------
### [GNURadio](https://www.gnuradio.org/)
&ensp;&ensp;GNURadio是GNU项目下的软件定义无线电工具包，其中包含了通用的信号采集模块，FFT模块，滤波器模块等用于处理数字信号的软件模块，通过GNURadio平台，我们可以很方便在在此基础之上构建我们所需的信号处理模块（只是一个平台，只能完成一些最基本的功能比如实时查看一个频段下的FFT，信号采集（无信号检测的采集）等基本功能）

### [gr-lora](https://github.com/rpp0/gr-lora)
&ensp;&ensp;基于GNURadio平台下的一个开源lora解码模块，其提供了基于LoRa信号的detect，time offset compensation, align, demodulate, decode功能，我们后续的实验都是通过修改gr-lora模块的代码来实现信号采集和实时信号处理的功能。

### [matlab](https://ww2.mathworks.cn/?s_tid=gn_logo)
&ensp;&ensp;目前实验室LoRa组最常用的实验仿真和信号处理软件，其提供了大量的数字信号处理工具包和便于调试的功能界面，方便快速进行offline的实验和debug。同时他的simulink工具包也可以实现GNURadio的作用，但是由于账号已经网上资源的限制，目前实验室还未进行simulink的使用。

<font color='red'>**PS: 其实可以直接使用GNURadio完成整个信号处理流程和实验，但是因为调试（信号处理的debug非常麻烦）方便，所以使用matlab来完成整个实验**</font>

## 安装/部署教程
--------------------
[GNURadio安装教程]()  
[gr-lora安装教程]()  
[USRP配置教程]()

## 使用介绍
--------------------
[如何使用gr-lora采集一个信号]()  
[如何使用gr-lora解析信号]()  
[如何使用matlab解析信号]()  
[如何使用USRP发送信号]()

## 详细讲解
--------------------
[gr-lora详细讲解]()  
[matlab详细讲解]()  
[USRP详细讲解]()  

