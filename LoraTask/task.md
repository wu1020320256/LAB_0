1、读懂发送方的代码，熟悉Lora信号的形式，学习Lora，什么是FHSS？

Lora数据：

（1）一个symbol传输SF个bit的数据，一个symbol是由2^SF个chip组成的

（2）循环误差编码：取决于CR的值，例如CR=4/5，表示每4个bit需要多出1bit的冗余位来填充，

FHSS：

跳频扩频通信，在预定的跳频周期结束后，即该部分数据发送完成，则发射机和接收机切换到跳频预定义的下一个信道，以便继续发送和接收数据包的下一部分内容。跳频周期由setFHSSHoppingPeriod（n）决定，n为symbol的整数倍，即n个symbol进行一次FHSS

代码解析：

（1）subchirp_num:将一个symbol分为4个短的chirp

（2）C_fre：系统的工作频率

（3）78：true_bin的作用??

（4）bin_part表示的是一个短的chirp由多少个chip组成，在本次实验中，bin_part=256

（5）252-254三行的作用？？？

setInvertIQ_downchirp的作用：将射频模块设置为下行扫频中的反相IQ模式

2、熟悉RadioLib.h库

3、学习gr-lora和rtl-sdr

4、询问凯文接收方gr-lora的修改

5、接收到第一个lora信号

6、对信号使用matlab进行处理