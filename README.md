### CMAP-dpeak

Peak deconvolution algorithm to extract gene expression levels. Outperforms original Harvard Medical School algorithm, +18% to accuracy and execution time.

> The goal of the Connectivity Map (CMap) is to build a resource of gene expression profiles to aid biomedical research and therapeutic discovery. The core CMap technology is an assay called L1000, which measures gene expression changes in cells. L1000 achieves a significant increase in throughput and a significant decrease in cost as compared to other gene expression technologies due in large part to the practice of measuring two genes using the same physical material (in this case, microscopic beads). After data collection, the expression values of the two genes are computationally extracted in a process called "peak deconvolution". The current implementation, an algorithm called D-Peak, is effective but imperfect, and hence the aim of this contest is to improve D-Peak, which we hypothesize will enable CMap to produce higher quality data more efficiently.

![1.png](https://lh3.googleusercontent.com/560LxnDXqbeYFCDt2OhLto8VHOCAfSIackvMa8ozpaXQltrJamOYI_2qYY8DwyJXDBkKO5ydt2fd-xjzI-XbLFbSmbXKitlER0FJFe8VIWsdc5fCkeOgyzWHlgCSvZt48Npw9pgM)

D-Peak algorithm schematics

![4.png](https://lh3.googleusercontent.com/-wtQ7pv4eLKwV3-AvHrkxmak3psQPbsHzCMEssbPxQYelMfqoPkNaLOGvTWjUGJmfbxq_DR_6Bau0MtfMNQhxk3VfAhAQ54R0erAP1s-XmFxRGIztRpPl7M-Wj89nGWngWWzY0xq)

D-Peak error modes (*)

![2.png](https://lh5.googleusercontent.com/IqTt4KC_thCmx27xGYcsXYRRzkLdcw3Vvx75HkheMt2myty3XbCB-CEgduqgCpY02s9gwQs2EsuByxoPyZys0QfV5Pf8U6C4Y1AHjjsL_pkyM3UQDDTOXYOexHFS4vWrGg5WzRY4)

Accuracy scoring schematic (*)

![3.png](https://lh6.googleusercontent.com/yZcyA6oOtnDa5D21z_dAfRxNHZx2pankRLsSyfmUt7RTIGuvSirCHwm7x1gJG5V1u45O9SlSaAEs5AaIMtvFJggMj5lxVAqghocRSYFAp-wTSP13xl8_AkNVLwqAJ2GulddgQC_B)

Luminex detection schematic (*)

-----

References:

* original algorithm: https://github.com/cmap/cmapM
* (*) contest (images and quotations source): https://www.topcoder.com/challenges/30076905
