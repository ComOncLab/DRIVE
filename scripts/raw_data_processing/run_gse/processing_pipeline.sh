#!/bin/bash
###第一个参数为 SRR IDs，如果有多个 SRR 就用逗号分割
###第二个参数是 输入SRA文件路径，所有的 SRR 都在这个里面
###第三个参数是 输出路径，该样本的输出路径
###第四个参数是 这个样本测序是单端还是双端
###第五个参数是 调用多少核来计算
###第六个参数是 参考基因组的index所在folder
###第七个参数是 salmon index
sra_id=$1
in_sra=$2
out_path=$3
is_pair="${4:-false}"  # 默认值为false
workers=$5
genome_idx=$6
salmon_idx=$7
sample=$(basename ${out_path})
# 验证输入
if [[ "$is_pair" != "true" && "$is_pair" != "false" ]]; then
    echo "错误：is_pair 参数必须是 'true' 或 'false'"
    exit 1
fi
###判断目标文件夹是否有 bam 文件，如果有就不需要运行
if find ${out_path} -maxdepth 1 -type f -name "*_Aligned.sortedByCoord.out.bam" -print -quit | grep -q .
then
    echo "存在 BAM 比对文件。"
    has_bam="Yes"
    exit 1 ##退出
fi
###判断是否有 Salmon 定量文件
if find ${out_path} -type f -name "gene_counts.rds" -print -quit | grep -q .
then
    echo "存在 Salmon 定量文件。"
    has_salmon="Yes"
fi
###判断是否有 HLA typing 文件
if find ${out_path} -type f -name "*genotype.json" -print -quit | grep -q .
then
    echo "存在 HLA typing 定量文件。"
    has_hla="Yes"
fi

##分割srr
IFS=',' read -r -a parts <<< "$sra_id"
###extract fastq
for ((i=0; i<${#parts[@]}; i++))
do
	echo "开始从 ${parts[i]} 中提取 Fastq..."
	mkdir ${out_path}/${parts[i]}
	cp ${in_sra}/${parts[i]}.sra ${out_path}/${parts[i]}/
	fasterq-dump ${out_path}/${parts[i]}/ -O ${out_path}/${parts[i]}/ -e ${workers} --split-files
	###将提取出的fastq移到上级样本文件夹下，并删除sra文件
	mv ${out_path}/${parts[i]}/*fastq ${out_path}/
	rm -rf ${out_path}/${parts[i]}
	echo "${parts[i]} 提取完成"
done
###fastp 质控
for ((i=0; i<${#parts[@]}; i++))
do
	echo "开始对 ${parts[i]} 进行质控..."
	#双端还是单端
	if [[ "$is_pair" == "true" ]]; then
	    fastp -w ${workers} -i ${out_path}/${parts[i]}_1.fastq -I ${out_path}/${parts[i]}_2.fastq -o ${out_path}/${parts[i]}_qc_1.fastq -O ${out_path}/${parts[i]}_qc_2.fastq -j ${out_path}/${parts[i]}.json -h ${out_path}/${parts[i]}.html
	elif [[ "$is_pair" == "false" ]]; then
	    fastp -w ${workers} -i ${out_path}/${parts[i]}.fastq -o ${out_path}/${parts[i]}_qc.fastq -j ${out_path}/${parts[i]}.json -h ${out_path}/${parts[i]}.html
	fi
	echo "${parts[i]} 质控完成"
done
###STAR 比对
echo "开始比对样本 ${sample} 到参考基因组 ${genome_idx}"
if [[ "$is_pair" == "true" ]]; then
	echo "双端测序"
	temp=${sra_id//,/_qc_1.fastq,${out_path}/}
	r1="${out_path}/${temp}_qc_1.fastq"
	temp=${sra_id//,/_qc_2.fastq,${out_path}/}
	r2="${out_path}/${temp}_qc_2.fastq"
	STAR --runThreadN ${workers} --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignIntronMax 299999 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 80000000000 --outSAMstrandField intronMotif --outSAMattributes All --limitSjdbInsertNsj 2000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --twopassMode Basic --outSAMmultNmax 1 --outFileNamePrefix ${out_path}/${sample}_ --genomeDir ${genome_idx} --readFilesIn ${r1} ${r2}
elif [[ "$is_pair" == "false" ]]; then
	echo "单端测序"
	temp=${sra_id//,/_qc.fastq,${out_path}/}
	r1="${out_path}/${temp}_qc.fastq"
	STAR --runThreadN ${workers} --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignIntronMax 299999 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 80000000000 --outSAMstrandField intronMotif --outSAMattributes All --limitSjdbInsertNsj 2000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --twopassMode Basic --outSAMmultNmax 1 --outFileNamePrefix ${out_path}/${sample}_ --genomeDir ${genome_idx} --readFilesIn ${r1}
fi
echo "比对完成"
###Salmon 定量表达
##最新版 1.10.3
if [ "$has_salmon" != "Yes" ]; then
    echo "未发现 Salmon 定量文件"
    mkdir ${out_path}/salmon_res
    echo "开始对样本 ${sample} 进行 Salmon 转录本定量"
    if [[ "$is_pair" == "true" ]]; then
	    r1=$(echo "${out_path}/"${sra_id//,/_qc_1.fastq ${out_path}/}"_qc_1.fastq")
	    r2=$(echo "${out_path}/"${sra_id//,/_qc_2.fastq ${out_path}/}"_qc_2.fastq")
	    conda run -n salmon salmon quant -i ${salmon_idx} -l A -1 ${r1} -2 ${r2} --validateMappings -o ${out_path}/salmon_res/
    elif [[ "$is_pair" == "false" ]]; then
	    r1=$(echo "${out_path}/"${sra_id//,/_qc.fastq ${out_path}/}"_qc.fastq")
	    conda run -n salmon salmon quant -i ${salmon_idx} -l A -r ${r1} --validateMappings -o ${out_path}/salmon_res/
    fi	
    echo "转录本定量完成"
    ###将转录本定量转化为基因定量（tximport）
    echo "将转录本表达转化为基因表达..."
    /home/data/R_442/bin/Rscript /home/wt/Drug_splicing/scripts/trans2gene.R -i ${out_path}/salmon_res/quant.sf -o ${out_path}/salmon_res/gene_counts.rds
    echo "基因表达转化完成"
fi
##删除所有的fastq文件
echo "删除所有 Fastq 文件..."
rm ${out_path}/*fastq
echo "删除完成"

####MHC typing
if [ "$has_hla" != "Yes" ]; then
	echo "进行 HLA Typing..."
	#conda run -n arcas-hla samtools reheader -c 'grep -v ^@HD' ${out_path}/${sample}_Aligned.sortedByCoord.out.bam > ${out_path}/${sample}.bam
	if [[ "$is_pair" == "true" ]]; then
		conda run -n arcas-hla /home/wt/software/arcasHLA/arcasHLA extract ${out_path}/${sample}_Aligned.sortedByCoord.out.bam -t 10 -v -o /home/wt/GEO_data/SRA/typing_out/ --temp /home/wt/GEO_data/SRA/typing_tmp/
		conda run -n arcas-hla /home/wt/software/arcasHLA/arcasHLA genotype /home/wt/GEO_data/SRA/typing_out/${sample}_Aligned.sortedByCoord.out.extracted.1.fq.gz /home/wt/GEO_data/SRA/typing_out/${sample}_Aligned.sortedByCoord.out.extracted.2.fq.gz -o /home/wt/GEO_data/SRA/typing_out/ --temp /home/wt/GEO_data/SRA/typing_tmp/ -v -t 10
	elif [[ "$is_pair" == "false" ]]; then
		conda run -n arcas-hla /home/wt/software/arcasHLA/arcasHLA extract ${out_path}/${sample}_Aligned.sortedByCoord.out.bam -t 10 -o /home/wt/GEO_data/SRA/typing_out/ --temp /home/wt/GEO_data/SRA/typing_tmp/ --single
		conda run -n arcas-hla /home/wt/software/arcasHLA/arcasHLA genotype /home/wt/GEO_data/SRA/typing_out/${sample}_Aligned.sortedByCoord.out.extracted.fq.gz -o /home/wt/GEO_data/SRA/typing_out/ --temp /home/wt/GEO_data/SRA/typing_tmp/ -t 10 --single
	fi
	mv /home/wt/GEO_data/SRA/typing_out/${sample}_Aligned.genotype.json ${out_path}/${sample}.json
	rm -rf /home/wt/GEO_data/SRA/typing_out/${sample}.*
	echo "HLA Typing 完成。"
fi
