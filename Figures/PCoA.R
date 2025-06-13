library(ConQuR)
library(doParallel)

Plot_PCoA(TAX=taxa, factor=batchid, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, main="ConQuR (Penalized), Bray-Curtis")

Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")