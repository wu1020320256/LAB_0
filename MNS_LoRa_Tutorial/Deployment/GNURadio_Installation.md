# 安装说明
--------------------
&ensp;&ensp;系统： _Ubuntu 20.04 TLS_  
&ensp;&ensp;Version： _GNURadio 3.8.2.0_  
&ensp;&ensp;Ps<font color='red'>：当时安装 gr-lora 时并不支持在 GNURadio3.9 平台上运行，所以安装3.8版本，但是gr-lora后续已经更新支持3.9了，同学们可以自己尝试部署最新版本的 gr-lora 和 GNURadio </font>

# 1.安装依赖项
&ensp;&ensp;可见：https://wiki.gnuradio.org/index.php/UbuntuInstall#Focal_Fossa_.2820.04.29

    sudo apt install git cmake g++ libboost-all-dev libgmp-dev swig python3-numpy \
    python3-mako python3-sphinx python3-lxml doxygen libfftw3-dev \
    libsdl1.2-dev libgsl-dev libqwt-qt5-dev libqt5opengl5-dev python3-pyqt5 \
    liblog4cpp5-dev libzmq3-dev python3-yaml python3-click python3-click-plugins \
    python3-zmq python3-scipy python3-gi python3-gi-cairo gobject-introspection gir1.2-gtk-3.0

# 2.安装Volk
&ensp;&ensp;VOLK是向量优化的内核库。它是一个包含用于不同数学运算的手写SIMD代码内核的库。由于每种SIMD体系结构都可能非常不同，并且还没有编译器能够正确或高效地处理矢量化，因此VOLK会以不同的方式解决此问题。  

    cd    //返回主目录
    git clone --recursive https://github.com/gnuradio/volk.git  //循环克隆git子项目
    （由于我的cpu用的是AMD的，有些CPU特性需要修改，需要参考README中的missing submodule，这里intel的cpu可以不加--recursive，因为可能会git不上）
    cd volk
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/usr/bin/python3 ../
    make
    make test
    sudo make install
    sudo ldconfig

# 3.安装GNU Radio 3.8版本

    cd
    git clone -b maint-3.8 https://github.com/gnuradio/gnuradio.git //选择3.8.2.0版本
    cd gnuradio
    git submodule update --init --recursive  
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/usr/bin/python3 ../
    make
    make test
    sudo make install
    sudo ldconfig

# 4.添加环境变量 PYTHONPATH 和 LD_LIBRARY_PATH
&ensp;&ensp;(1).找到安装目录{install-prefix}：gnuradio-config-info --prefix  
&ensp;&ensp;(2).找到python library{Py-version}：find {install-prefix} -name gnuradio | grep "packages"  
&ensp;&ensp;(3).添加路径到bash start-up file：  
sudo gedit ~/.bashrc 在文件的末尾添加：  

    export PYTHONPATH={install-prefix}/lib/{Py-version}/dist-packages:$PYTHONPATH  
    export LD_LIBRARY_PATH={install-prefix}/lib:$LD_LIBRARY_PATH
    Ps:其中{install-prefix}和{Py-version}是你分别找到的安装目录和python目录

# 5.检查
&ensp;&ensp;查看 GNU Radio 版本：gnuradio-config-info --version  
&ensp;&ensp;GNU Radio 运行：gnuradio-companion，有软件界面显示正常