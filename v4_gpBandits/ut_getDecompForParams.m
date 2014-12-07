% Unit test for getDecompForParams

boParams.decompStrategy = 'known';
[numGroups, decomp, boAddParams] = getDecompForParams(10, 3,boParams),decomp{1},
[numGroups, decomp, boAddParams] = getDecompForParams(24, 8, boParams),decomp{1},
[numGroups, decomp, boAddParams] = getDecompForParams(100, 8, boParams),decomp{1},


boParams.decompStrategy = 'partialLearn';
[numGroups, decomp, boAddParams] = getDecompForParams(10, 3, boParams),
[numGroups, decomp, boAddParams] = getDecompForParams(24, 8, boParams),
[numGroups, decomp, boAddParams] = getDecompForParams(100, 8, boParams),

