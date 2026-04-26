while IFS=" " read -r name
do
	rm -rf /home/wt/GEO_data/SRA/gse_out/${name}/jcast_out
	mkdir /home/wt/GEO_data/SRA/gse_out/${name}/jcast_out
	jcast -o /home/wt/GEO_data/SRA/gse_out/${name}/jcast_out/ -q 0 0.1 /home/wt/GEO_data/SRA/gse_out/${name}/rmats_filter_out/ /home/wt/Genome_data/Homo_sapiens.GRCh38.115.gtf /home/wt/Genome_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	#cat /home/wt/GEO_data/SRA/gse_out/${name}/jcast_out/psq_T* > /home/wt/GEO_data/SRA/gse_out/${name}/jcast_out/all_seq.fasta
	echo "${name} done"
done < "test"
