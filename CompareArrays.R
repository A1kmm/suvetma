full <- read.table('microarray_test_results_fulltfbs.txt')
perturbed <- read.table('microarray_test_results_perturbedtfbs.txt')
perturbedRatios <- (perturbed[,2] / perturbed[,3])[full[,3] != 0]
fullRatios <- (full[,2] / full[,3])[full[,3] != 0]
