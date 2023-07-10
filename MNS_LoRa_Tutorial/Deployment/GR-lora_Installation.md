# 安装说明
&ensp;&ensp;系统： Ubuntu 20.04 TLS  
&ensp;&ensp;依赖项：GNU Radio 3.8.2.0   
&ensp;&ensp;Version： rpp0/gr-lora 0.6.2

# 安装liquid-dsp：
&ensp;&ensp;依赖项安装，一个dsp库，用于进行一些快速浮点运算，复数计算，gr-lora中使用到了    
&ensp;&ensp;https://github.com/jgaeddert/liquid-dsp

    git clone git://github.com/jgaeddert/liquid-dsp.git
    ./bootstrap.sh
    ./configure
    make
    sudo make install
    sudo ldconfig

# 安装rtl-sdr驱动
&ensp;&ensp;[RTL-SDR_Driver](./RTL-SDR_Driver_Installation.md)


# 安装gr-osmocom
&ensp;&ensp;依赖项安装，提供对rtl-sdr等一些硬件设备的开发包，被gr-lora用于控制rtl-sdr

    git clone -b v0.2.2 git://git.osmocom.org/gr-osmosdr
    cd gr-osmosdr/
    mkdir build
    cd build/
    cmake ../ -DENABLE_RTL=ON  # 很重要
    make
    sudo make install
    sudo ldconfig


# 安装gr-lora

    git clone https://github.com/rpp0/gr-lora.git  # 注意版本
    mkdir build   
    cd build
    cmake ../  # Note to Arch Linux users: add "-DCMAKE_INSTALL_PREFIX=/usr"
    make
    make test
    sudo make install
    sudo ldconfig
    //由于之前安装GNU Radio的时候已经添加好了 PYTHONPATH 和 LD_LIBRARY_PATH 路径，所以不需要再添加路径了
