# 安装rtl-sdr驱动
1. 可能的更新：
```
sudo apt-get update
```
2. 安装需要的工具:git,cmake,build-essential
```
sudo apt-get install git
sudo apt-get install cmake
sudo apt-get install build-essential
```
3. 安装用于接入USB设备的通用C库：libusb-1.0-0-dev
```
sudo apt-get install libusb-1.0-0-dev
```
4. 获取并编译构建 RTL2832U Osmocom 驱动
```
git clone git://git.osmocom.org/rtl-sdr.git
cd rtl-sdr/
mkdir build
cd build
cmake ../ -DINSTALL_UDEV_RULES=ON
(出现问题：--Could NOT find PkgConfig(missing:PKG_CONFIG_EXECUTABLE)
尝试：sudo apt-get install --reinstall pkg-config cmake-data)
make
sudo make install 
sudo ldconfig
sudo cp ../rtl-sdr.rules /etc/udev/rules.d/
```
5. 将可能导致冲突的驱动列入黑名单
```
cd /etc/modprobe.d/
sudo gedit blacklist-rtl.conf
(添加一行：blacklist dvb_usb_rtl28xxu并保存)
```
6. 测试驱动
```
rtl_test -t
```