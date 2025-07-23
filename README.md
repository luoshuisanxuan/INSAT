This folder contains the code for pure GPS single-point positioning and INSAT. To run INSAT, follow these steps:

1. Set `settings.INSAT` to 0, then run `init` to obtain variables such as ephemeris and frame headers.
2. Set `settings.INSAT` to 1, then run `init` again to start INSAT and get the spoofing suppression results.

The schematic diagram of the algorithm is shown as follows:

<img src="https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAJCjmiA_0ImdVJtlQMrqUnT3UOpL8QLAAKDFgAD1ghU7sgldn7s8rY2BA.png" alt="PixPin_2025-07-23_23-26-42.png" style="zoom:50%;" />