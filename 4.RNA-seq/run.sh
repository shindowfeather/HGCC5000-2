for i in `ls *.fastq.gz |sed 's/_R1_001.fastq.gz//g'|sed 's/_R2_001.fastq.gz//g'|uniq`;do bash trim.sh $i & done 
for i in `ls *R1_001.trim.fastq.gz|sed 's/_R1_001.trim.fastq.gz//g'`;do echo "bash mapping.sh $i" >> runmapping.sh 
bash runmapping.sh 
for i in `ls *R1_001.trim.fastq.gz|sed 's/_R1_001.trim.fastq.gz//g'`;do echo "bash samtools2.sh $i" >> run.samtools.sh ;done
nohup cat run.samtools.sh|xargs -n 2 -P 6 bash &
ls *gtf > gtf_list

stringtie --merge -p 40 -G /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/database/database/hg38/genenome/hg38/gene/hg38.refGene.gtf -o merged.gtf gtf_list

for i in `ls *gtf |grep -v merged |sed 's/.gtf//g'`;do bash stringtie2.sh $i & done

for i in `ls */*.gtf|sed 's/.gtf//g'`;do cat ${i}.gtf|grep -v chrM|grep -v chrY > ${i}.gtf2 & done

cat gtf_list |sed 's/.gtf//g' > li
ls */*.gtf2|grep -v merge > li2
paste li li2 > list2
python /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/pipeline/prepDE3.py -i list2

/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/gffcompare/gffcompare  -r /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/database/database/hg38/genenome/hg38/gene/hg38.refGene.gtf -G merged.gtf
#cat gene_count_matrix.csv|sed '1d'|cat header4 - > gene_count_matrix2.csv
perl convert_stringtie_gene_to_reference_gene.pl gffcmp.merged.gtf.refmap gene_count_matrix.csv |grep -v replicates |grep -v unknown > gene_count_matrix.annotation.txt
Rscript sum.R 

#DEG analysis 
Rscript edger.without_replicates.P13A.R
le P13Acomparison.txt|sed '1d'|sort -k5g |awk '{OFS="\t"} { if ($2>0) printf "%s\t%4.3e\n", $1, 1/$5 ;else printf "%s\t%4.3e\n", $1, -1/$5 }'|sort -k2gr > P13Acomparison.rnk
#GSEA hallmark analysis 
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh  GSEAPreranked -rnk P13Acomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh  GSEAPreranked -rnk P13Acomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c7.immunesigdb.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh  GSEAPreranked -rnk P13Acomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -rnk P13Acomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.symbols.gmt -collapse false 

Rscript edger.without_replicates.P13E.R
le P13Ecomparison.txt|sed '1d'|sort -k5g |awk '{OFS="\t"} { if ($2>0) printf "%s\t%4.3e\n", $1, 1/$5 ;else printf "%s\t%4.3e\n", $1, -1/$5 }'|sort -k2gr > P13Ecomparison.rnk
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh  GSEAPreranked -rnk P13Ecomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -rnk P13Ecomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c7.immunesigdb.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -rnk P13Ecomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.symbols.gmt -collapse false
/proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -rnk P13Ecomparison.rnk -gmx /proj/sens2023512/nobackup/wharf/lucycy/lucycy-sens2023512/software/GSEA/msigdb_v2023.2.Hs_files_to_download_locally/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt -collapse false
