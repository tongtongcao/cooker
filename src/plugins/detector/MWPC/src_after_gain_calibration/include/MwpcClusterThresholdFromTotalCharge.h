const double totChargeCutC2XL[56]={37.2726, 82.9899, 128.707, 174.425, 220.142, 265.859, 311.576, 357.294, 403.011, 448.728, 494.446, 540.163, 585.88, 631.598, 677.315, 723.032, 768.75, 814.467, 860.184, 905.901, 951.619, 997.336, 1043.05, 1088.77, 1134.49, 1180.21, 1225.92, 1271.64, 1317.36, 1363.07, 1408.79, 1454.51, 1500.23, 1545.94, 1591.66, 1637.38, 1683.1, 1728.81, 1774.53, 1820.25, 1865.96, 1911.68, 1957.4, 2003.12, 2048.83, 2094.55, 2140.27, 2185.99, 2231.7, 2277.42, 2323.14, 2368.86, 2414.57, 2460.29, 2506.01, 2551.72};
const double totChargeCutErrC2XL[56]={23.838, 24.7495, 26.1983, 28.1015, 30.3736, 32.9385, 35.7332, 38.708, 41.8243, 45.0529, 48.3713, 51.7623, 55.2123, 58.7111, 62.2505, 65.8237, 69.4258, 73.0523, 76.6998, 80.3655, 84.0469, 87.7422, 91.4495, 95.1676, 98.8951, 102.631, 106.375, 110.125, 113.882, 117.644, 121.411, 125.182, 128.958, 132.738, 136.521, 140.307, 144.096, 147.888, 151.683, 155.48, 159.279, 163.08, 166.883, 170.688, 174.494, 178.302, 182.111, 185.922, 189.734, 193.547, 197.361, 201.177, 204.993, 208.811, 212.629, 216.448};
const double totChargeCutC2XR[56]={48.0734, 92.3153, 136.557, 180.799, 225.041, 269.283, 313.525, 357.767, 402.009, 446.251, 490.493, 534.735, 578.977, 623.219, 667.461, 711.703, 755.945, 800.187, 844.429, 888.671, 932.913, 977.155, 1021.4, 1065.64, 1109.88, 1154.12, 1198.36, 1242.61, 1286.85, 1331.09, 1375.33, 1419.57, 1463.82, 1508.06, 1552.3, 1596.54, 1640.78, 1685.03, 1729.27, 1773.51, 1817.75, 1861.99, 1906.24, 1950.48, 1994.72, 2038.96, 2083.2, 2127.45, 2171.69, 2215.93, 2260.17, 2304.41, 2348.66, 2392.9, 2437.14, 2481.38};
const double totChargeCutErrC2XR[56]={28.0862, 29.3504, 31.3442, 33.9394, 37.0096, 40.4468, 44.1654, 48.1002, 52.2023, 56.4352, 60.7717, 65.191, 69.6774, 74.2188, 78.8055, 83.4303, 88.0869, 92.7708, 97.4778, 102.205, 106.949, 111.709, 116.482, 121.267, 126.062, 130.866, 135.678, 140.498, 145.325, 150.158, 154.996, 159.839, 164.687, 169.539, 174.394, 179.253, 184.116, 188.981, 193.849, 198.72, 203.593, 208.468, 213.345, 218.224, 223.105, 227.987, 232.871, 237.757, 242.644, 247.532, 252.422, 257.312, 262.204, 267.097, 271.99, 276.885};
const double totChargeCutC2YL[16]={78.4454, 102.253, 126.062, 149.87, 173.678, 197.486, 221.294, 245.102, 268.91, 292.718, 316.526, 340.334, 364.142, 387.95, 411.759, 435.567};
const double totChargeCutErrC2YL[16]={40.0048, 41.6002, 44.1311, 47.4481, 51.3993, 55.8502, 60.6909, 65.8356, 71.2184, 76.7892, 82.5099, 88.3515, 94.2915, 100.312, 106.4, 112.545};
const double totChargeCutC2YR[16]={40.0797, 88.5549, 137.03, 185.505, 233.98, 282.455, 330.931, 379.406, 427.881, 476.356, 524.831, 573.306, 621.781, 670.256, 718.732, 767.207};
const double totChargeCutErrC2YR[16]={27.9018, 30.656, 34.7648, 39.8109, 45.4835, 51.5762, 57.9568, 64.5398, 71.2693, 78.1073, 85.0277, 92.0119, 99.0464, 106.121, 113.229, 120.363};
const double totChargeCutC3XL[64]={36.6442, 66.8923, 97.1404, 127.388, 157.637, 187.885, 218.133, 248.381, 278.629, 308.877, 339.125, 369.373, 399.621, 429.869, 460.117, 490.365, 520.613, 550.861, 581.11, 611.358, 641.606, 671.854, 702.102, 732.35, 762.598, 792.846, 823.094, 853.342, 883.59, 913.838, 944.086, 974.334, 1004.58, 1034.83, 1065.08, 1095.33, 1125.57, 1155.82, 1186.07, 1216.32, 1246.57, 1276.82, 1307.06, 1337.31, 1367.56, 1397.81, 1428.06, 1458.3, 1488.55, 1518.8, 1549.05, 1579.3, 1609.54, 1639.79, 1670.04, 1700.29, 1730.54, 1760.78, 1791.03, 1821.28, 1851.53, 1881.78, 1912.02, 1942.27};
const double totChargeCutErrC3XL[64]={24.2296, 25.9965, 28.7007, 32.1062, 36.0146, 40.2797, 44.7998, 49.5051, 54.3475, 59.2935, 64.3191, 69.407, 74.5446, 79.7221, 84.9323, 90.1696, 95.4294, 100.708, 106.003, 111.312, 116.633, 121.964, 127.305, 132.653, 138.009, 143.37, 148.738, 154.11, 159.487, 164.868, 170.253, 175.641, 181.032, 186.426, 191.822, 197.221, 202.622, 208.025, 213.429, 218.836, 224.244, 229.653, 235.064, 240.476, 245.889, 251.304, 256.719, 262.135, 267.553, 272.971, 278.39, 283.81, 289.23, 294.651, 300.073, 305.495, 310.918, 316.342, 321.766, 327.19, 332.615, 338.041, 343.467, 348.893};
const double totChargeCutC3XR[64]={38.8024, 69.0901, 99.3778, 129.666, 159.953, 190.241, 220.529, 250.816, 281.104, 311.392, 341.68, 371.967, 402.255, 432.543, 462.83, 493.118, 523.406, 553.694, 583.981, 614.269, 644.557, 674.844, 705.132, 735.42, 765.708, 795.995, 826.283, 856.571, 886.859, 917.146, 947.434, 977.722, 1008.01, 1038.3, 1068.58, 1098.87, 1129.16, 1159.45, 1189.74, 1220.02, 1250.31, 1280.6, 1310.89, 1341.17, 1371.46, 1401.75, 1432.04, 1462.33, 1492.61, 1522.9, 1553.19, 1583.48, 1613.76, 1644.05, 1674.34, 1704.63, 1734.91, 1765.2, 1795.49, 1825.78, 1856.07, 1886.35, 1916.64, 1946.93};
const double totChargeCutErrC3XR[64]={26.9271, 29.5116, 33.3774, 38.1368, 43.4974, 49.2633, 55.308, 61.5494, 67.9332, 74.4228, 80.9928, 87.6252, 94.3067, 101.028, 107.781, 114.56, 121.361, 128.181, 135.016, 141.865, 148.725, 155.596, 162.475, 169.362, 176.255, 183.155, 190.06, 196.969, 203.884, 210.801, 217.723, 224.648, 231.575, 238.505, 245.438, 252.373, 259.31, 266.248, 273.189, 280.131, 287.075, 294.02, 300.966, 307.914, 314.862, 321.812, 328.763, 335.714, 342.667, 349.62, 356.574, 363.529, 370.485, 377.441, 384.397, 391.355, 398.313, 405.271, 412.23, 419.189, 426.149, 433.109, 440.07, 447.031};
const double totChargeCutC3YL[16]={40.7279, 79.6215, 118.515, 157.409, 196.302, 235.196, 274.09, 312.983, 351.877, 390.77, 429.664, 468.558, 507.451, 546.345, 585.238, 624.132};
const double totChargeCutErrC3YL[16]={25.3884, 27.8778, 31.5938, 36.1602, 41.2959, 46.8139, 52.594, 58.5586, 64.6568, 70.8539, 77.1263, 83.4568, 89.8333, 96.2465, 102.69, 109.157};
const double totChargeCutC3YR[16]={46.6103, 75.7098, 104.809, 133.909, 163.008, 192.108, 221.207, 250.307, 279.406, 308.505, 337.605, 366.704, 395.804, 424.903, 454.003, 483.102};
const double totChargeCutErrC3YR[16]={25.9377, 26.7182, 27.9705, 29.6351, 31.6468, 33.9441, 36.4729, 39.1886, 42.0549, 45.0431, 48.1304, 51.2991, 54.5348, 57.8264, 61.1649, 64.5429};
const double totChargeCutC4XL[72]={60.8444, 76.4919, 92.1394, 107.787, 123.434, 139.082, 154.729, 170.377, 186.024, 201.672, 217.319, 232.967, 248.614, 264.262, 279.909, 295.557, 311.204, 326.852, 342.499, 358.147, 373.794, 389.442, 405.089, 420.737, 436.384, 452.032, 467.679, 483.327, 498.974, 514.622, 530.269, 545.917, 561.564, 577.212, 592.859, 608.507, 624.154, 639.802, 655.449, 671.097, 686.744, 702.392, 718.039, 733.687, 749.334, 764.982, 780.629, 796.277, 811.924, 827.572, 843.219, 858.867, 874.514, 890.162, 905.809, 921.457, 937.104, 952.752, 968.399, 984.047, 999.694, 1015.34, 1030.99, 1046.64, 1062.28, 1077.93, 1093.58, 1109.23, 1124.87, 1140.52, 1156.17, 1171.82};
const double totChargeCutErrC4XL[72]={52.5909, 54.8576, 58.4404, 63.1155, 68.6602, 74.8815, 81.625, 88.7716, 96.2316, 103.938, 111.839, 119.896, 128.081, 136.37, 144.745, 153.192, 161.7, 170.259, 178.863, 187.505, 196.18, 204.885, 213.614, 222.366, 231.138, 239.927, 248.732, 257.552, 266.384, 275.228, 284.082, 292.946, 301.818, 310.698, 319.586, 328.48, 337.38, 346.287, 355.198, 364.114, 373.035, 381.96, 390.889, 399.821, 408.757, 417.696, 426.638, 435.583, 444.531, 453.481, 462.433, 471.388, 480.344, 489.303, 498.263, 507.226, 516.19, 525.155, 534.122, 543.09, 552.06, 561.031, 570.004, 578.977, 587.952, 596.927, 605.904, 614.882, 623.86, 632.84, 641.82, 650.802};
const double totChargeCutC4XR[72]={31.4642, 59.6361, 87.8081, 115.98, 144.152, 172.324, 200.496, 228.668, 256.84, 285.012, 313.184, 341.356, 369.528, 397.7, 425.872, 454.044, 482.216, 510.388, 538.56, 566.732, 594.903, 623.075, 651.247, 679.419, 707.591, 735.763, 763.935, 792.107, 820.279, 848.451, 876.623, 904.795, 932.967, 961.139, 989.311, 1017.48, 1045.65, 1073.83, 1102, 1130.17, 1158.34, 1186.51, 1214.69, 1242.86, 1271.03, 1299.2, 1327.37, 1355.55, 1383.72, 1411.89, 1440.06, 1468.23, 1496.41, 1524.58, 1552.75, 1580.92, 1609.09, 1637.27, 1665.44, 1693.61, 1721.78, 1749.95, 1778.13, 1806.3, 1834.47, 1862.64, 1890.81, 1918.99, 1947.16, 1975.33, 2003.5, 2031.67};
const double totChargeCutErrC4XR[72]={27.6021, 29.6852, 32.8649, 36.8585, 41.4313, 46.4125, 51.684, 57.1657, 62.8024, 68.5561, 74.3995, 80.313, 86.2823, 92.2966, 98.3475, 104.429, 110.535, 116.663, 122.809, 128.971, 135.146, 141.332, 147.529, 153.735, 159.949, 166.169, 172.396, 178.629, 184.867, 191.109, 197.355, 203.605, 209.859, 216.116, 222.375, 228.637, 234.901, 241.168, 247.437, 253.708, 259.98, 266.254, 272.53, 278.806, 285.085, 291.364, 297.645, 303.927, 310.21, 316.493, 322.778, 329.064, 335.35, 341.637, 347.925, 354.213, 360.503, 366.792, 373.083, 379.374, 385.665, 391.957, 398.249, 404.542, 410.835, 417.129, 423.423, 429.718, 436.013, 442.308, 448.603, 454.899};
const double totChargeCutC4YL[16]={28.5017, 66.5246, 104.547, 142.57, 180.593, 218.616, 256.639, 294.662, 332.685, 370.708, 408.731, 446.754, 484.776, 522.799, 560.822, 598.845};
const double totChargeCutErrC4YL[16]={21.3222, 22.2382, 23.6863, 25.5763, 27.8182, 30.3342, 33.0617, 35.9526, 38.9705, 42.0882, 45.285, 48.5454, 51.8572, 55.2114, 58.6005, 62.0189};
const double totChargeCutC4YR[16]={32.6953, 62.7766, 92.8579, 122.939, 153.021, 183.102, 213.183, 243.264, 273.346, 303.427, 333.508, 363.59, 393.671, 423.752, 453.834, 483.915};
const double totChargeCutErrC4YR[16]={29.9724, 31.7784, 34.5794, 38.157, 42.3145, 46.8981, 51.7948, 56.9237, 62.2276, 67.6652, 73.2068, 78.8304, 84.5197, 90.2623, 96.0486, 101.871};