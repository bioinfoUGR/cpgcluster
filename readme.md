CpGcluster
----------

**CpGcluster website:** http://bioinfo2.ugr.es/CpGcluster/

Despite their involvement in the regulation of gene expression and their importance as genomic markers for promoter prediction, no objective standard exists for defining CpG islands (CGIs), since all current approaches rely on a large parameter space formed by the thresholds of length, CpG fraction and G+C content.

Given the higher frequency of CpG dinucleotides at CGIs, as compared to bulk DNA, the distance distributions between neighboring CpGs should differ for bulk and island CpGs. The algorithm CpGcluster [1, 2], based on the physical distance between neighboring CpGs on the chromosome is able to predict directly clusters of CpGs, while not depending on the subjective criteria mentioned above. By assigning a p-value to each of these clusters, the most statistically significant ones can be predicted as CGIs. CpGcluster was benchmarked against five other CGI finders by using a test sequence set assembled from an experimental CGI library. CpGcluster reached the highest overall accuracy values, while showing the lowest rate of false-positive predictions. Since a minimum-length threshold is not required, CpGcluster can find short but fully functional CGIs usually missed by other algorithms. The CGIs predicted by CpGcluster present the lowest degree of overlap with Alu retrotransposons and, simultaneously, the highest overlap with vertebrate Phylogenetic Conserved Elements (PhastCons). CpGcluster’s CGIs overlapping with the Transcription Start Site (TSS) show the highest statistical significance, as compared to the islands in other genome locations, thus qualifying CpGcluster as a valuable tool in discriminating functional CGIs from the remaining islands in the bulk genome.

CpGcluster uses only integer arithmetic, thus being a fast and computationally efficient algorithm able to predict statistically significant clusters of CpG dinucleotides. Another outstanding feature is that all predicted CGIs start and end with a CpG dinucleotide, which should be appropriate for a genomic feature whose functionality is based precisely on CpG dinucleotides. The only search parameter in CpGcluster is the distance between two consecutive CpGs, in contrast to previous algorithms. Therefore, none of the main statistical properties of CpG islands (neither G+C content, CpG fraction nor length threshold) are needed as search parameters, which may lead to the high specificity and low overlap with spurious Alu elements observed for CpGcluster predictions.

[1] Hackenberg M, Previti C, Luque-Escamilla PL, Carpena P, Martínez-Aroza J, Oliver JL. 2006.
CpGcluster: A distance-based algorithm for CpG-island detection.
BMC Bioinformatics 7: 446
http://dx.doi.org/10.1186/1471-2105-7-446

[2] http://bioinfo2.ugr.es/CpGcluster/