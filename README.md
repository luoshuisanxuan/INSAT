This folder contains the code for pure GPS single-point positioning and INSAT. To run INSAT, follow these steps:

1. Set `settings.INSAT` to 0, then run `init` to obtain variables such as ephemeris and frame headers.
2. Set `settings.INSAT` to 1, then run `init` again to start INSAT and get the spoofing suppression results.

The schematic diagram of the algorithm is shown as follows:

![PixPin_2025-07-23_23-26-42.png](https://youke1.picui.cn/s1/2025/07/23/6881039f62ea3.png)

