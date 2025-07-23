运行cleanDynamic.dat





此文件夹包含纯GPS单点定位与INSAT代码，运行INSAT需要按照如下步骤进行：

1、将INSoff = 1，运行init得到星历、帧头等变量；

2、运行INStrj得到运动状态；

3、将INSoff = 0，运行init开始INSAT并得到欺骗抑制结果



算法原理图



捕获代码需要改善

运行



