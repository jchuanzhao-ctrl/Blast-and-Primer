这段代码是一个用于自动化设计PCR引物的Python脚本。它从Excel文件中读取短DNA序列，通过BLAST比对定位其在参考基因组（如GRCh38）中的位置，提取上下游各300bp的序列片段，使用primer3设计引物，并通过BLAST验证引物特异性。最终将结果（包括染色体位置、引物序列、Tm值、特异性等）输出到Excel文件中。

This code is an automated PCR primer design tool. It reads short DNA sequences from an Excel file, locates their positions in a reference genome (e.g., GRCh38) using BLAST, extracts 600bp flanking sequences (300bp upstream/downstream), designs primers with primer3, and validates specificity via BLAST. Results (including chromosomal position, primer sequences, Tm values, specificity, etc.) are exported to an Excel file.
