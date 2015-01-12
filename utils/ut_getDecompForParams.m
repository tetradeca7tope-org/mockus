% Unit test for getDecompForParams

boParams.decompStrategy = 'known';
[decomp, numGroups, boAddParams] = getDecompForParams(10, 3,boParams),decomp{1},
[decomp, numGroups, boAddParams] = getDecompForParams(24, 8, boParams),decomp{1},
[decomp, numGroups, boAddParams] = getDecompForParams(100, 8, boParams),decomp{1},


boParams.decompStrategy = 'partialLearn';
[decomp, numGroups, boAddParams] = getDecompForParams(10, 3, boParams),
[decomp, numGroups, boAddParams] = getDecompForParams(24, 8, boParams),
[decomp, numGroups, boAddParams] = getDecompForParams(100, 8, boParams),

