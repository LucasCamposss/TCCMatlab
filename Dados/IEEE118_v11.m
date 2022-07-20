 % TITLE 
 % ****   A 118 Bus IEEE Test System  ****  
 
% Bus data 
% NB  T G   VT    Angle      PG     QG       QMIN   QMAX    PLOAD    QLOAD     QBAR  
 DBAR = [
    1  1 0 0.955  10.67     0.0     0.0      -5.      15.    51.0     27.0      0.0 
    2  0 0 0.971  11.22     0.0     0.0       0.       0.    20.0      9.0      0.0 
    3  0 0 0.968  11.56     0.0     0.0       0.       0.    39.0     10.0      0.0 
    4  1 1 0.998  15.28    -9.0     0.0    -300.     300.    30.0     12.0      0.0 
    5  0 0 1.002  15.73     0.0     0.0       0.       0.     0.0      0.0     -40.
    6  1 0 0.990  13.00     0.0     0.0     -13.      50.    52.0     22.0      0.0 
    7  0 0 0.989  12.56     0.0     0.0       0.       0.    19.0      2.0      0.0 
    8  1 1 1.015  20.77   -28.0     0.0    -300.     300.     0.0      0.0      0.0 
    9  0 0 1.043  28.02     0.0     0.0       0.       0.     0.0      0.0      0.0 
   10  1 1 1.050  35.61   450.0     0.0    -147.     200.     0.0      0.0      0.0 
   11  0 0 0.985  12.72     0.0     0.0       0.       0.    70.0     23.0      0.0 
   12  1 0 0.990  12.20    85.0     0.0     -35.     120.    47.0     10.0      0.0 
   13  0 0 0.968  11.35     0.0     0.0       0.       0.    34.0     16.0      0.0 
   14  0 0 0.984  11.50     0.0     0.0       0.       0.    14.0      1.0      0.0 
   15  1 0 0.970  11.23     0.0     0.0     -10.      30.    90.0     30.0      0.0 
   16  0 0 0.984  11.91     0.0     0.0       0.       0.    25.0     10.0      0.0 
   17  0 0 0.995  13.74     0.0     0.0       0.       0.    11.0      3.0      0.0 
   18  1 0 0.973  11.53     0.0     0.0     -16.      50.    60.0     34.0      0.0 
   19  1 0 0.963  11.05     0.0     0.0      -8.      24.    45.0     25.0      0.0 
   20  0 0 0.958  11.93     0.0     0.0       0.       0.    18.0      3.0      0.0 
   21  0 0 0.959  13.52     0.0     0.0       0.       0.    14.0      8.0      0.0 
   22  0 0 0.970  16.08     0.0     0.0       0.       0.    10.0      5.0      0.0 
   23  0 0 1.000  21.00     0.0     0.0       0.       0.     7.0      3.0      0.0 
   24  1 1 0.992  20.89   -13.0     0.0    -300.     300.     0.0      0.0      0.0 
   25  1 0 1.050  27.93   220.0     0.0     -47.     140.     0.0      0.0      0.0 
   26  1 1 1.015  29.71   314.0     0.0   -1000.    1000.     0.0      0.0      0.0 
   27  1 1 0.968  15.35    -9.0     0.0    -300.     300.    62.0     13.0      0.0 
   28  0 0 0.962  13.62     0.0     0.0       0.       0.    17.0      7.0      0.0 
   29  0 0 0.963  12.63     0.0     0.0       0.       0.    24.0      4.0      0.0 
   30  0 0 0.968  18.79     0.0     0.0       0.       0.     0.0      0.0      0.0 
   31  1 1 0.967  12.75     7.0     0.0    -300.     300.    43.0     27.0      0.0 
   32  1 0 0.964  14.80     0.0     0.0     -14.      42.    59.0     23.0      0.0 
   33  0 0 0.972  10.63     0.0     0.0       0.       0.    23.0      9.0      0.0 
   34  1 0 0.986  11.30     0.0     0.0      -8.      24.    59.0     26.0      14.
   35  0 0 0.981  10.87     0.0     0.0       0.       0.    33.0      9.0      0.0 
   36  1 0 0.980  10.87     0.0     0.0      -8.      24.    31.0     17.0      0.0 
   37  0 0 0.992  11.77     0.0     0.0       0.       0.     0.0      0.0     -25.
   38  0 0 0.962  16.91     0.0     0.0       0.       0.     0.0      0.0      0.0 
   39  0 0 0.970   8.41     0.0     0.0       0.       0.    27.0     11.0      0.0 
   40  1 1 0.970   7.35   -46.0     0.0    -300.     300.    20.0     23.0      0.0 
   41  0 0 0.967   6.92     0.0     0.0       0.       0.    37.0     10.0      0.0 
   42  1 1 0.985   8.53   -59.0     0.0    -300.     300.    37.0     23.0      0.0 
   43  0 0 0.978  11.28     0.0     0.0       0.       0.    18.0      7.0      0.0 
   44  0 0 0.985  13.82     0.0     0.0       0.       0.    16.0      8.0      10.
   45  0 0 0.987  15.67     0.0     0.0       0.       0.    53.0     22.0      10.
   46  1 0 1.005  18.49    19.0     0.0    -100.     100.    28.0     10.0      10.
   47  0 0 1.017  20.73     0.0     0.0       0.       0.    34.0      0.0      0.0 
   48  0 0 1.021  19.93     0.0     0.0       0.       0.    20.0     11.0      15.
   49  1 0 1.025  20.94   204.0     0.0     -85.     210.    87.0     30.0      0.0 
   50  0 0 1.001  18.90     0.0     0.0       0.       0.    17.0      4.0      0.0 
   51  0 0 0.967  16.28     0.0     0.0       0.       0.    17.0      8.0      0.0 
   52  0 0 0.957  15.32     0.0     0.0       0.       0.    18.0      5.0      0.0 
   53  0 0 0.946  14.35     0.0     0.0       0.       0.    23.0     11.0      0.0 
   54  1 1 0.955  15.26    48.0     0.0    -300.     300.   113.0     32.0      0.0 
   55  1 0 0.952  14.97     0.0     0.0      -8.      23.    63.0     22.0      0.0 
   56  1 0 0.954  15.16     0.0     0.0      -8.      15.    84.0     18.0      0.0 
   57  0 0 0.971  16.36     0.0     0.0       0.       0.    12.0      3.0      0.0 
   58  0 0 0.959  15.51     0.0     0.0       0.       0.    12.0      3.0      0.0 
   59  1 0 0.985  19.37   155.0     0.0     -60.     180.   277.0    113.0      0.0 
   60  0 0 0.993  23.15     0.0     0.0       0.       0.    78.0      3.0      0.0 
   61  1 1 0.995  24.04   160.0     0.0    -100.     300.     0.0      0.0      0.0 
   62  1 0 0.998  23.43     0.0     0.0     -20.      20.    77.0     14.0      0.0 
   63  0 0 0.969  22.75     0.0     0.0       0.       0.     0.0      0.0      0.0 
   64  0 0 0.984  24.52     0.0     0.0       0.       0.     0.0      0.0      0.0 
   65  1 0 1.005  27.65   391.0     0.0     -67.     200.     0.0      0.0      0.0 
   66  1 0 1.050  27.48   392.0     0.0     -67.     200.    39.0     18.0      0.0 
   67  0 0 1.020  24.84     0.0     0.0       0.       0.    28.0      7.0      0.0 
   68  0 0 1.003  27.55     0.0     0.0       0.       0.     0.0      0.0      0.0 
   69  2 1 1.035  30.00   516.4     0.0    -300.     300.     0.0      0.0      0.0 
   70  1 0 0.984  22.58     0.0     0.0     -10.      32.    66.0     20.0      0.0 
   71  0 0 0.987  22.15     0.0     0.0       0.       0.     0.0      0.0      0.0 
   72  1 0 0.980  20.98   -12.0     0.0    -100.     100.     0.0      0.0      0.0 
   73  1 0 0.991  21.94    -6.0     0.0    -100.     100.     0.0      0.0      0.0 
   74  1 0 0.958  21.64     0.0     0.0      -6.       9.    68.0     27.0      12.
   75  0 0 0.967  22.91     0.0     0.0       0.       0.    47.0     11.0      0.0 
   76  1 0 0.943  21.77     0.0     0.0      -8.      23.    68.0     36.0      0.0 
   77  1 0 1.006  26.72     0.0     0.0     -20.      70.    61.0     28.0      0.0 
   78  0 0 1.003  26.42     0.0     0.0       0.       0.    71.0     26.0      0.0 
   79  0 0 1.009  26.72     0.0     0.0       0.       0.    39.0     32.0      20.
   80  1 1 1.040  28.96   477.0     0.0    -165.     280.   130.0     26.0      0.0 
   81  0 0 0.997  28.10     0.0     0.0       0.       0.     0.0      0.0      0.0 
   82  0 0 0.989  27.24     0.0     0.0       0.       0.    54.0     27.0      20.
   83  0 0 0.985  28.42     0.0     0.0       0.       0.    20.0     10.0      10.
   84  0 0 0.980  30.95     0.0     0.0       0.       0.    11.0      7.0      0.0 
   85  1 0 0.985  32.51     0.0     0.0      -8.      23.    24.0     15.0      0.0 
   86  0 0 0.987  31.14     0.0     0.0       0.       0.    21.0     10.0      0.0 
   87  1 0 1.015  31.40     4.0     0.0    -100.    1000.     0.0      0.0      0.0 
   88  0 0 0.987  35.64     0.0     0.0       0.       0.    48.0     10.0      0.0 
   89  1 1 1.005  39.69   607.0     0.0    -210.     300.     0.0      0.0      0.0 
   90  1 1 0.985  33.29   -85.0     0.0    -300.     300.    78.0     42.0      0.0 
   91  1 0 0.980  33.31   -10.0     0.0    -100.     100.     0.0      0.0      0.0 
   92  1 0 0.993  33.80     0.0     0.0      -3.       9.    65.0     10.0      0.0 
   93  0 0 0.987  30.79     0.0     0.0       0.       0.    12.0      7.0      0.0 
   94  0 0 0.991  28.64     0.0     0.0       0.       0.    30.0     16.0      0.0 
   95  0 0 0.981  27.67     0.0     0.0       0.       0.    42.0     31.0      0.0 
   96  0 0 0.993  27.51     0.0     0.0       0.       0.    38.0     15.0      0.0 
   97  0 0 1.011  27.88     0.0     0.0       0.       0.    15.0      9.0      0.0 
   98  0 0 1.024  27.40     0.0     0.0       0.       0.    34.0      8.0      0.0 
   99  1 1 1.010  27.04   -42.0     0.0    -100.     100.     0.0      0.0      0.0 
  100  1 0 1.017  28.03   252.0     0.0     -50.     155.    37.0     18.0      0.0 
  101  0 0 0.993  29.61     0.0     0.0       0.       0.    22.0     15.0      0.0 
  102  0 0 0.991  32.30     0.0     0.0       0.       0.     5.0      3.0      0.0 
  103  1 0 1.001  24.44    40.0     0.0     -15.      40.    23.0     16.0      0.0 
  104  1 0 0.971  21.69     0.0     0.0      -8.      23.    38.0     25.0      0.0 
  105  1 0 0.965  20.57     0.0     0.0      -8.      23.    31.0     26.0      20.
  106  0 0 0.962  20.32     0.0     0.0       0.       0.    43.0     16.0      0.0 
  107  1 1 0.952  17.53   -22.0     0.0    -200.     200.    28.0     12.0      6.0
  108  0 0 0.967  19.38     0.0     0.0       0.       0.     2.0      1.0      0.0 
  109  0 0 0.967  18.93     0.0     0.0       0.       0.     8.0      3.0      0.0 
  110  1 0 0.973  18.09     0.0     0.0      -8.      23.    39.0     30.0      6.0
  111  1 1 0.980  19.74    36.0     0.0    -100.    1000.     0.0      0.0      0.0 
  112  1 1 0.975  14.99   -43.0     0.0    -100.    1000.    25.0     13.0      0.0 
  113  1 0 0.993  13.74    -6.0     0.0    -100.     200.     0.0      0.0      0.0 
  114  0 0 0.960  14.46     0.0     0.0       0.       0.     8.0      3.0      0.0 
  115  0 0 0.960  14.46     0.0     0.0       0.       0.    22.0      7.0      0.0 
  116  1 1 1.005  27.12  -184.0     0.0   -1000.    1000.     0.0      0.0      0.0 
  117  0 0 0.974  10.67     0.0     0.0       0.       0.    20.0      8.0      0.0 
  118  0 0 0.949  21.92     0.0     0.0       0.       0.    33.0     15.0      0.0 
                                                                              ];
% Line Data
% From        T0    r             x          Bsh     
 DLIN = [
    1         2      3.030      9.990       2.540        0.0     0.0      0.0     0.0     0.0     0.0
    1         3      1.290      4.240       1.082        0.0     0.0      0.0     0.0     0.0     0.0
    4         5      0.176      0.798       0.210        0.0     0.0      0.0     0.0     0.0     0.0
    3         5      2.410     10.800       2.840        0.0     0.0      0.0     0.0     0.0     0.0
    5         6      1.190      5.400       1.426        0.0     0.0      0.0     0.0     0.0     0.0
    6         7      0.459      2.080       0.550        0.0     0.0      0.0     0.0     0.0     0.0
    8         9      0.244      3.050      116.20        0.0     0.0      0.0     0.0     0.0     0.0
    8         5      0.000      2.670         0.0      0.985     0.0      0.0     0.0     0.0     0.0
    9        10      0.258      3.220      123.00        0.0     0.0      0.0     0.0     0.0     0.0
    4        11      2.090      6.880       1.748        0.0     0.0      0.0     0.0     0.0     0.0
    5        11      2.030      6.820       1.738        0.0     0.0      0.0     0.0     0.0     0.0
   11        12      0.595      1.960       0.502        0.0     0.0      0.0     0.0     0.0     0.0
    2        12      1.870      6.160       1.572        0.0     0.0      0.0     0.0     0.0     0.0
    3        12      4.840     16.000       4.060        0.0     0.0      0.0     0.0     0.0     0.0
    7        12      0.862      3.400       0.874        0.0     0.0      0.0     0.0     0.0     0.0
   11        13      2.225      7.310       1.876        0.0     0.0      0.0     0.0     0.0     0.0
   12        14      2.150      7.070       1.816        0.0     0.0      0.0     0.0     0.0     0.0
   13        15      7.440     24.440       6.268        0.0     0.0      0.0     0.0     0.0     0.0
   14        15      5.950     19.500       5.020        0.0     0.0      0.0     0.0     0.0     0.0
   12        16      2.120      8.340       2.140        0.0     0.0      0.0     0.0     0.0     0.0
   15        17      1.320      4.370       4.440        0.0     0.0      0.0     0.0     0.0     0.0
   16        17      4.540     18.010       4.660        0.0     0.0      0.0     0.0     0.0     0.0
   17        18      1.230      5.050       1.298        0.0     0.0      0.0     0.0     0.0     0.0
   18        19      1.119      4.930       1.142        0.0     0.0      0.0     0.0     0.0     0.0
   19        20      2.520     11.700       2.980        0.0     0.0      0.0     0.0     0.0     0.0
   15        19      1.200      3.940       1.010        0.0     0.0      0.0     0.0     0.0     0.0
   20        21      1.830      8.490       2.160        0.0     0.0      0.0     0.0     0.0     0.0
   21        22      2.090      9.700       2.460        0.0     0.0      0.0     0.0     0.0     0.0
   22        23      3.420     15.900       4.040        0.0     0.0      0.0     0.0     0.0     0.0
   23        24      1.350      4.920       4.980        0.0     0.0      0.0     0.0     0.0     0.0
   23        25      1.560      8.000       8.640        0.0     0.0      0.0     0.0     0.0     0.0
   26        25      0.000      3.820         0.0      0.960     0.0      0.0     0.0     0.0     0.0
   25        27      3.180     16.300      17.640        0.0     0.0      0.0     0.0     0.0     0.0
   27        28      1.913      8.550       2.160        0.0     0.0      0.0     0.0     0.0     0.0
   28        29      2.370      9.430       2.380        0.0     0.0      0.0     0.0     0.0     0.0
   30        17      0.000      3.880         0.0      0.960     0.0      0.0     0.0     0.0     0.0
    8        30      0.431      5.040      51.400        0.0     0.0      0.0     0.0     0.0     0.0
   26        30      0.799      8.600      90.800        0.0     0.0      0.0     0.0     0.0     0.0
   17        31      4.740     15.630       3.990        0.0     0.0      0.0     0.0     0.0     0.0
   29        31      1.080      3.310       0.830        0.0     0.0      0.0     0.0     0.0     0.0
   23        32      3.170     11.530      11.730        0.0     0.0      0.0     0.0     0.0     0.0
   31        32      2.980      9.850       2.510        0.0     0.0      0.0     0.0     0.0     0.0
   27        32      2.290      7.550       1.926        0.0     0.0      0.0     0.0     0.0     0.0
   15        33      3.800     12.440       3.194        0.0     0.0      0.0     0.0     0.0     0.0
   19        34      7.520     24.700       6.320        0.0     0.0      0.0     0.0     0.0     0.0
   35        36      0.224      1.020       0.268        0.0     0.0      0.0     0.0     0.0     0.0
   35        37      1.100      4.970       1.318        0.0     0.0      0.0     0.0     0.0     0.0
   33        37      4.150     14.200       3.660        0.0     0.0      0.0     0.0     0.0     0.0
   34        36      0.871      2.680       0.568        0.0     0.0      0.0     0.0     0.0     0.0
   34        37      0.256      0.940       0.984        0.0     0.0      0.0     0.0     0.0     0.0
   38        37      0.000      3.750         0.0      0.935     0.0      0.0     0.0     0.0     0.0
   37        39      3.210     10.600       2.700        0.0     0.0      0.0     0.0     0.0     0.0
   37        40      5.930     16.800       4.200        0.0     0.0      0.0     0.0     0.0     0.0
   30        38      0.464      5.400      42.200        0.0     0.0      0.0     0.0     0.0     0.0
   39        40      1.840      6.050       1.552        0.0     0.0      0.0     0.0     0.0     0.0
   40        41      1.450      4.870       1.222        0.0     0.0      0.0     0.0     0.0     0.0
   40        42      5.550     18.300       4.660        0.0     0.0      0.0     0.0     0.0     0.0
   41        42      4.100     13.500       3.440        0.0     0.0      0.0     0.0     0.0     0.0
   43        44      6.080     24.540       6.068        0.0     0.0      0.0     0.0     0.0     0.0
   34        43      4.130     16.810       4.226        0.0     0.0      0.0     0.0     0.0     0.0
   44        45      2.240      9.010       2.240        0.0     0.0      0.0     0.0     0.0     0.0
   45        46      4.000     13.560       3.320        0.0     0.0      0.0     0.0     0.0     0.0
   46        47      3.800     12.700       3.160        0.0     0.0      0.0     0.0     0.0     0.0
   46        48      6.010     18.900       4.720        0.0     0.0      0.0     0.0     0.0     0.0
   47        49      1.910      6.250       1.604        0.0     0.0      0.0     0.0     0.0     0.0
   42        49      7.150     32.300       8.600        0.0     0.0      0.0     0.0     0.0     0.0
   42        49      7.150     32.300       8.600        0.0     0.0      0.0     0.0     0.0     0.0
   45        49      6.840     18.600       4.440        0.0     0.0      0.0     0.0     0.0     0.0
   48        49      1.790      5.050       1.258        0.0     0.0      0.0     0.0     0.0     0.0
   49        50      2.670      7.520       1.874        0.0     0.0      0.0     0.0     0.0     0.0
   49        51      4.860     13.700       3.420        0.0     0.0      0.0     0.0     0.0     0.0
   51        52      2.030      5.880       1.396        0.0     0.0      0.0     0.0     0.0     0.0
   52        53      4.050     16.350       4.058        0.0     0.0      0.0     0.0     0.0     0.0
   53        54      2.630     12.200       3.100        0.0     0.0      0.0     0.0     0.0     0.0
   49        54      7.300     28.900       7.380        0.0     0.0      0.0     0.0     0.0     0.0
   49        54      8.690     29.100       7.300        0.0     0.0      0.0     0.0     0.0     0.0
   54        55      1.690      7.070       2.020        0.0     0.0      0.0     0.0     0.0     0.0
   54        56      0.275      0.955       0.732        0.0     0.0      0.0     0.0     0.0     0.0
   55        56      0.488      1.510       0.374        0.0     0.0      0.0     0.0     0.0     0.0
   56        57      3.430      9.660       2.420        0.0     0.0      0.0     0.0     0.0     0.0
   50        57      4.740     13.400       3.320        0.0     0.0      0.0     0.0     0.0     0.0
   56        58      3.430      9.660       2.420        0.0     0.0      0.0     0.0     0.0     0.0
   51        58      2.550      7.190       1.788        0.0     0.0      0.0     0.0     0.0     0.0
   54        59      5.030     22.930       5.980        0.0     0.0      0.0     0.0     0.0     0.0
   56        59      8.250     25.100       5.690        0.0     0.0      0.0     0.0     0.0     0.0
   56        59      8.030     23.900       5.360        0.0     0.0      0.0     0.0     0.0     0.0
   55        59      4.739     21.580       5.646        0.0     0.0      0.0     0.0     0.0     0.0
   59        60      3.170     14.500       3.760        0.0     0.0      0.0     0.0     0.0     0.0
   59        61      3.280     15.000       3.880        0.0     0.0      0.0     0.0     0.0     0.0
   60        61      0.264      1.350       1.456        0.0     0.0      0.0     0.0     0.0     0.0
   60        62      1.230      5.610       1.468        0.0     0.0      0.0     0.0     0.0     0.0
   61        62      0.824      3.760       0.980        0.0     0.0      0.0     0.0     0.0     0.0
   63        59      0.000      3.860         0.0      0.960     0.0      0.0     0.0     0.0     0.0
   63        64      0.172      2.000      21.600        0.0     0.0      0.0     0.0     0.0     0.0
   64        61      0.000      2.680         0.0      0.985     0.0      0.0     0.0     0.0     0.0
   38        65      0.901      9.860      104.60        0.0     0.0      0.0     0.0     0.0     0.0
   64        65      0.269      3.020      38.000        0.0     0.0      0.0     0.0     0.0     0.0
   49        66      1.800      9.190       2.480        0.0     0.0      0.0     0.0     0.0     0.0
   49        66      1.800      9.190       2.480        0.0     0.0      0.0     0.0     0.0     0.0
   62        66      4.820     21.800       5.780        0.0     0.0      0.0     0.0     0.0     0.0
   62        67      2.580     11.700       3.100        0.0     0.0      0.0     0.0     0.0     0.0
   65        66      0.000      3.700         0.0      0.935     0.0      0.0     0.0     0.0     0.0
   66        67      2.240     10.150       2.682        0.0     0.0      0.0     0.0     0.0     0.0
   65        68      0.138      1.600      63.800        0.0     0.0      0.0     0.0     0.0     0.0
   47        69      8.440     27.780       7.092        0.0     0.0      0.0     0.0     0.0     0.0
   49        69      9.850     32.400       8.280        0.0     0.0      0.0     0.0     0.0     0.0
   68        69      0.000      3.700         0.0      0.935     0.0      0.0     0.0     0.0     0.0
   69        70      3.000     12.700      12.200        0.0     0.0      0.0     0.0     0.0     0.0
   24        70      0.221     41.150      10.198        0.0     0.0      0.0     0.0     0.0     0.0
   70        71      0.882      3.550       0.878        0.0     0.0      0.0     0.0     0.0     0.0
   24        72      4.880     19.600       4.880        0.0     0.0      0.0     0.0     0.0     0.0
   71        72      4.460     18.000       4.444        0.0     0.0      0.0     0.0     0.0     0.0
   71        73      0.866      4.540       1.178        0.0     0.0      0.0     0.0     0.0     0.0
   70        74      4.010     13.230       3.368        0.0     0.0      0.0     0.0     0.0     0.0
   70        75      4.280     14.100       3.600        0.0     0.0      0.0     0.0     0.0     0.0
   69        75      4.050     12.200      12.400        0.0     0.0      0.0     0.0     0.0     0.0
   74        75      1.230      4.060       1.034        0.0     0.0      0.0     0.0     0.0     0.0
   76        77      4.440     14.800       3.680        0.0     0.0      0.0     0.0     0.0     0.0
   69        77      3.090     10.100      10.380        0.0     0.0      0.0     0.0     0.0     0.0
   75        77      6.010     19.990       4.978        0.0     0.0      0.0     0.0     0.0     0.0
   77        78      0.376      1.240       1.264        0.0     0.0      0.0     0.0     0.0     0.0
   78        79      0.546      2.440       0.648        0.0     0.0      0.0     0.0     0.0     0.0
   77        80      1.700      4.850       4.720        0.0     0.0      0.0     0.0     0.0     0.0
   77        80      2.940     10.500       2.280        0.0     0.0      0.0     0.0     0.0     0.0
   79        80      1.560      7.040       1.870        0.0     0.0      0.0     0.0     0.0     0.0
   68        81      0.175      2.020      80.800        0.0     0.0      0.0     0.0     0.0     0.0
   81        80      0.000      3.700         0.0      0.935     0.0      0.0     0.0     0.0     0.0
   77        82      2.980      8.530       8.174        0.0     0.0      0.0     0.0     0.0     0.0
   82        83      1.120      3.665       3.796        0.0     0.0      0.0     0.0     0.0     0.0
   83        84      6.250     13.200       2.580        0.0     0.0      0.0     0.0     0.0     0.0
   83        85      4.300     14.800       3.480        0.0     0.0      0.0     0.0     0.0     0.0
   84        85      3.020      6.410       1.234        0.0     0.0      0.0     0.0     0.0     0.0
   85        86      3.500     12.300       2.760        0.0     0.0      0.0     0.0     0.0     0.0
   86        87      2.828     20.740       4.450        0.0     0.0      0.0     0.0     0.0     0.0
   85        88      2.000     10.200       2.760        0.0     0.0      0.0     0.0     0.0     0.0
   85        89      2.390     17.300       4.700        0.0     0.0      0.0     0.0     0.0     0.0
   88        89      1.390      7.120       1.934        0.0     0.0      0.0     0.0     0.0     0.0
   89        90      5.180     18.800       5.280        0.0     0.0      0.0     0.0     0.0     0.0
   89        90      2.380      9.970      10.600        0.0     0.0      0.0     0.0     0.0     0.0
   90        91      2.540      8.360       2.140        0.0     0.0      0.0     0.0     0.0     0.0
   89        92      0.990      5.050       5.480        0.0     0.0      0.0     0.0     0.0     0.0
   89        92      3.930     15.810       4.140        0.0     0.0      0.0     0.0     0.0     0.0
   91        92      3.870     12.720       3.268        0.0     0.0      0.0     0.0     0.0     0.0
   92        93      2.580      8.480       2.180        0.0     0.0      0.0     0.0     0.0     0.0
   92        94      4.810     15.800       4.060        0.0     0.0      0.0     0.0     0.0     0.0
   93        94      2.230      7.320       1.876        0.0     0.0      0.0     0.0     0.0     0.0
   94        95      1.320      4.340       1.110        0.0     0.0      0.0     0.0     0.0     0.0
   80        96      3.560     18.200       4.940        0.0     0.0      0.0     0.0     0.0     0.0
   82        96      1.620      5.300       5.440        0.0     0.0      0.0     0.0     0.0     0.0
   94        96      2.690      8.690       2.300        0.0     0.0      0.0     0.0     0.0     0.0
   80        97      1.830      9.340       2.540        0.0     0.0      0.0     0.0     0.0     0.0
   80        98      2.380     10.800       2.860        0.0     0.0      0.0     0.0     0.0     0.0
   80        99      4.540     20.600       5.460        0.0     0.0      0.0     0.0     0.0     0.0
   92       100      6.480     29.500       4.720        0.0     0.0      0.0     0.0     0.0     0.0
   94       100      1.780      5.800       6.040        0.0     0.0      0.0     0.0     0.0     0.0
   95        96      1.710      5.470       1.474        0.0     0.0      0.0     0.0     0.0     0.0
   96        97      1.730      8.850       2.400        0.0     0.0      0.0     0.0     0.0     0.0
   98       100      3.970     17.900       4.760        0.0     0.0      0.0     0.0     0.0     0.0
   99       100      1.800      8.130       2.160        0.0     0.0      0.0     0.0     0.0     0.0
  100       101      2.770     12.620       3.280        0.0     0.0      0.0     0.0     0.0     0.0
   92       102      1.230      5.590       1.464        0.0     0.0      0.0     0.0     0.0     0.0
  101       102      2.460     11.200       2.940        0.0     0.0      0.0     0.0     0.0     0.0
  100       103      1.600      5.250       5.360        0.0     0.0      0.0     0.0     0.0     0.0
  100       104      4.510     20.400       5.410        0.0     0.0      0.0     0.0     0.0     0.0
  103       104      4.660     15.840       4.070        0.0     0.0      0.0     0.0     0.0     0.0
  103       105      5.350     16.250       4.080        0.0     0.0      0.0     0.0     0.0     0.0
  100       106      6.050     22.900       6.200        0.0     0.0      0.0     0.0     0.0     0.0
  104       105      0.994      3.780       0.986        0.0     0.0      0.0     0.0     0.0     0.0
  105       106      1.400      5.470       1.434        0.0     0.0      0.0     0.0     0.0     0.0
  105       107      5.300     18.300       4.720        0.0     0.0      0.0     0.0     0.0     0.0
  105       108      2.610      7.030       1.844        0.0     0.0      0.0     0.0     0.0     0.0
  106       107      5.300     18.300       4.720        0.0     0.0      0.0     0.0     0.0     0.0
  108       109      1.050      2.880       0.760        0.0     0.0      0.0     0.0     0.0     0.0
  103       110      3.906     18.130       4.610        0.0     0.0      0.0     0.0     0.0     0.0
  109       110      2.780      7.620       2.020        0.0     0.0      0.0     0.0     0.0     0.0
  110       111      2.200      7.550       2.000        0.0     0.0      0.0     0.0     0.0     0.0
  110       112      2.470      6.400       6.200        0.0     0.0      0.0     0.0     0.0     0.0
   17       113      0.913      3.010       0.768        0.0     0.0      0.0     0.0     0.0     0.0
   32       113      6.150     20.300       5.180        0.0     0.0      0.0     0.0     0.0     0.0
   32       114      1.350      6.120       1.628        0.0     0.0      0.0     0.0     0.0     0.0
   27       115      1.640      7.410       1.972        0.0     0.0      0.0     0.0     0.0     0.0
  114       115      0.230      1.040       0.276        0.0     0.0      0.0     0.0     0.0     0.0
   68       116      0.034      0.405      16.400        0.0     0.0      0.0     0.0     0.0     0.0
   12       117      3.290     14.000       3.580        0.0     0.0      0.0     0.0     0.0     0.0
   75       118      1.450      4.810       1.198        0.0     0.0      0.0     0.0     0.0     0.0
   76       118      1.640      5.440       1.356        0.0     0.0      0.0     0.0     0.0     0.0
                                                                                  ];
% [NL,NCOL]    = size(DLIN);
% DLIN(:,10)   = 100*ones(NL,1);
% DLIN(121,10) = 85;

% Dados para Limites de Circuitos - DVIO

DVIO = [ 3    50
         4   100
         5    60
         8   110
         11   60
         15   50
         16   50
         21   50
         22   50
         23   80
         26   80
         30  250
         33  100
         36  180
         37  330
         41  120
         50  100
         51   80
         54  330
         62   90
         78  290
         96  500
         97  160
         104 630
         107 380
         108 230
         111 110
         116 180
         119 320
         121 100
         123 120
         128 160
         148 100
         152 120
         176 80];


% Choose the lower and upper limits of voltage
 DTEN = [
 0  0.80  1.20
 1  0.90  1.05
 ];

% Active Power Generation Data
% Bus  PGmin  PGmax  Cost($/Mwh)
DGER = [ 
    1    0    130     1
    2    0    1000    400
    3    0    1000    400
    4    0    1000    400
    5   -100  0      -200
    6    0    150     6
    7    0    1000    400
    8    0    470     8
    9    0    1000    400
   10    0    1       10
   11    0    1000    400
   12    0    1000    400
   13    0    1000    400
   14  -100   0      -200
   15    0    450     15
   16    0    1000    400
   17    0    1000    400
   18    0    1000    400
   19    0    1000    400
   20    0    1000    400
   21    0    1000    400
   22    0    1000    400
   23    0    1000    400
   24    0    250     24
   25    0    1       25
   26  -100   0      -200
   27    0    1000    400  
   28    0    1000    400
   29    0    1000    400
   31    0    1000    400 
   32    0    1000    400  
   33    0    1000    400  
   34    0    1000    400
   35    0    1000    400  
   36    0    1000    400  
   37  -100   0      -200
   39    0    1000    400  
   40    0    1000    400
   41    0    1000    400  
   42    0    140     42
   43    0    1000    400  
   44  -100   0      -200
   45    0    1000    400  
   46    0    120     46
   47    0    1000    400  
   48    0    1000    400  
   49    0    450     49
   50    0    1000    400  
   51    0    1000    400  
   52  -100   0      -200
   53    0    1000    400
   54    0    570     54
   55    0    1000    400  
   56    0    1000    400  
   59    0    1000    400
   60    0    1000    400  
   61  -100   0      -200
   62    0    1000    400  
   65    0    1       65
   66    0    1000    400  
   67    0    1000    400  
   69    0    1000    69
   70    0    1000    400  
   72    0    1       72
   73    0    1       73
   74    0    1000    400  
   75    0    1000    400  
   76    0    1000    400  
   77    0    1000    400  
   78    0    1000    400  
   79    0    1000    400  
   80    0    1000    80
   82    0    1000    400
   83    0    1000    400
   84    0    1000    400
   85    0    300     85
   86    0    1000    400
   87  -100   0      -200
   88    0    1000    400
   89    0    1       89
   90    0    250     90
   91  -100   0      -200
   92    0    1000    400
   93    0    1000    400
   94    0    1000    400
   95    0    1000    400
   96    0    1000    400
   97    0    1000    400
   98    0    1000    400
   99  -100   0      -200
  100    0    1000    400
  101    0    1000    400
  102    0    1000    400
  103    0    1000    400
  104    0    1000    400
  105    0    1000    400
  106    0    1000    400
  107    0    1000    400
  108    0    1000    400
  109    0    1000    400
  110    0    1000    400
  111  -100   0      -200
  112    0    1000    400
  113    0    1       113
  114    0    1000    400
  115    0    1000    400
  116    0    1       116
  117    0    1000    400
  118    0    1000    400
                             ]; 
 
DGW = [
         5  0.0 100
        14  0.0 100
        26  0.0 100
        37  0.0 100
        44  0.0 100
        52  0.0 100
        61  0.0 100
        87  0.0 100
        91  0.0 100
        99  0.0 100
       111  0.0 100];
    
% % Choose the Limit of powerflow 
% FLG_LIM = 1; % Considerar limite de LT
% FLG_LIM = 0; % Não Considerar limite de LT

FLG_LIM = 1;

% % Choose the Objective Function
% 
 OBJF = 1;   % Choose OBJF = 1 for minimal operational cost 
%            % Choose OBJF = 2 for minimal transmission losses
%            
% % Choose number 1 to PrntScr 
% 
 R_BAR = 0;   % R_BAR = 1; Imprime Relatório de Barra 
%             % R_BAR = 0; Não Imprime Relatório de Barra 
% 
 R_GER = 0;   % R_GER = 1; Imprime Relatório de Geração 
%             % R_GER = 0; Não Imprime Relatório de Geração 
% 
 R_LIN = 0;   % R_BAR = 1; Imprime Relatório de Linha 
%             % R_BAR = 0; Não Imprime Relatório de Linha 
         
