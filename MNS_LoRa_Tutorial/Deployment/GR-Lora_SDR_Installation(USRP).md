# 1.安装GNURadio
&ensp;&ensp;[GNURadio_Installation](./GNURadio_Installation.md)

# 2.安装更新pip
&ensp;&ensp;通过pip+pybombs方便安装后续的依赖包
```
sudo apt install python3-pip
sudo pip install -U pip
```

# 3.安装pybombs
```
sudo pip install pybombs
```

# 4.获取安装库
```
pybombs recipes add gr-recipes git+https://github.com/gnuradio/gr-recipes.git
pybombs recipes add gr-etcetera git+https://github.com/gnuradio/gr-etcetera.git
```

# 5.安装相关库
```
sudo apt-get install gqrx-sdr
sudo pybombs install hackrf rtl-sdr gr-osmosdr osmo-sdr
注：其中安装到gnuradio部分能跳过（到了gnuradio部分已经够了），不用在这里安装完GNURadio，要不然会和之前安装gr-lora环境下的GNURadio版本冲突
```

# 6.下载uhd镜像文件
```
sudo uhd_images_downloader
```

# 7.安装gr-lora_sdr
```
git clone https://github.com/jkadbear/gr-lora.git
mkdir build
cd build
cmake ../
make
sudo make install
sudo ldconfig
```