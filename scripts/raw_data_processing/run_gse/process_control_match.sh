#!/bin/bash
# 参数
#第一个参数为对照组的样本，用逗号分割
#第二个参数为包含对照组，处理组的meta文件，如gse_1000.txt，三列：IDs,.treat_sample. control samples
#第三个参数为储存样本meta信息的文本文件，四列：GSM3093315 SRR6985601 Homo_sapiens PAIRED
con_samples=$1
con_treat_meta=$2
meta_file=$3
echo "开始处理对照组..."
IFS=',' read -r -a parts <<< "$con_samples"
####extract fastq
for ((i=0; i<${#parts[@]}; i++))
do
	# 查找匹配行（精确匹配第一列）
	line=$(grep "^${parts[i]}\s" ${meta_file} || awk -v key="${parts[i]}" '$1 == key' ${meta_file})
	# 如果找到匹配行
	if [ -n "$line" ]; then
		# 读取后三列
    		srrs=$(echo "$line" | awk '{print $2}')
    		organ=$(echo "$line" | awk '{print $3}')
    		pair=$(echo "$line" | awk '{print $4}')
		# 检查文件夹是否存在
		gsm=/home/wt/GEO_data/SRA/out/${parts[i]}/
		if [ ! -d "$gsm" ]; then
			echo "创建样本输出文件夹: $gsm"
   			mkdir $gsm
		else
   			echo "文件夹已存在: $gsm"
		fi
		if [[ "$pair" == "PAIRED" ]]; then is_pair="true"; else is_pair="false"; fi
		star_index="/home/wt/Genome_data/STAR_index/"
		salmon_index="/home/wt/Genome_data/salmon_index_21/human_transcripts_index/"
		bash processing_pipeline.sh ${srrs} /home/wt/GEO_data/SRA/sra/ ${gsm} ${is_pair} 15 ${star_index} ${salmon_index} 
	else
   		echo "样本 ${parts[i]} 未在 Metainfo 文件中找到。"
	fi
done
echo "对照组样本处理完成。"

#依据control样本，匹配处理组
while IFS=' ' read -r ids treat control; do
    # 检查第三列是否包含搜索字符
    if [[ "$control" == *"$con_samples"* ]]; then
        bash processing_gse.sh ${treat} ${control} ${meta_file} ${ids}
    fi
done < "$con_treat_meta"

echo "删除对照组 BAM 文件..."
IFS=',' read -r -a parts <<< "$con_samples"
###extract fastq
for ((i=0; i<${#parts[@]}; i++))
do
	rm /home/wt/GEO_data/SRA/out/${parts[i]}/*bam
	echo "已删除样本 ${parts[i]} 的 BAM 文件"
done
echo "运行完成"

