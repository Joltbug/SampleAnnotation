sort -k1,1 -k2,2n GSE97887_occ1_MAINSPLITS.Ast.diffPeaks.bed > GSE97887_Ast_sorted.bed
sort -k1,1 -k2,2n GSE97887_occ1_MAINSPLITS.End.diffPeaks.bed > GSE97887_End_sorted.bed
sort -k1,1 -k2,2n GSE97887_occ1_MAINSPLITS.Ex.diffPeaks.bed > GSE97887_Ex_sorted.bed
bedtools intersect -a sample0.bed -b GSE97887_Ast_sorted.bed -f 0.5 -wa> intersect_0toAst_0.5.bed
bedtools intersect -a sample0.bed -b GSE97887_End_sorted.bed -f 0.5 -wa> intersect_0toEnd_0.5.bed
bedtools intersect -a sample0.bed -b GSE97887_Ex_sorted.bed -f 0.5 -wa> intersect_0toEx_0.5.bed
bedtools intersect -b sample0.bed -a GSE97887_Ast_sorted.bed -f 0.5 -wa> intersect_Astto0_0.5.bed
bedtools intersect -b sample0.bed -a GSE97887_End_sorted.bed -f 0.5 -wa> intersect_Endto0_0.5.bed
bedtools intersect -b sample0.bed -a GSE97887_Ex_sorted.bed -f 0.5 -wa> intersect_Exto0_0.5.bed