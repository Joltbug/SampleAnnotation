bedtools intersect -a sample0.bed -b sample1.bed -f 0.5 -wa> intersect_0to1_0.5.bed
bedtools intersect -b sample0.bed -a sample1.bed -f 0.5 -wa> intersect_1to0_0.5.bed
bedtools intersect -a sample0.bed -b sample1.bed -f 0.3 -wa> intersect_0to1_0.3.bed
bedtools intersect -b sample0.bed -a sample1.bed -f 0.3 -wa> intersect_1to0_0.3.bed
bedtools intersect -a sample0.bed -b sample1.bed -f 0.7 -wa> intersect_0to1_0.7.bed
bedtools intersect -b sample0.bed -a sample1.bed -f 0.7 -wa> intersect_1to0_0.7.bed