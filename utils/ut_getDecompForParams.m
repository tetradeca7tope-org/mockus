% Unit test for getDecompForParams

boParams.decompStrategy = 'known';
[decomp, numGroups, boAddParams] = getDecompForParams(10, 3,boParams  ,true),decomp{1},
[decomp, numGroups, boAddParams] = getDecompForParams(24, 8, boParams ,true),decomp{1},
[decomp, numGroups, boAddParams] = getDecompForParams(100, 8, boParams,true),decomp{1},


boParams.decompStrategy = 'partialLearn';
[decomp, numGroups, boAddParams] = getDecompForParams(10, 3, boParams ,true),
[decomp, numGroups, boAddParams] = getDecompForParams(24, 8, boParams ,true),
[decomp, numGroups, boAddParams] = getDecompForParams(100, 8, boParams,true),

