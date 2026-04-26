#!/bin/bash
# 参数
#第一个参数为处理组的样本，用逗号分割
#第二个参数为对照组的样本，用逗号分割
#第三个参数为储存样本meta信息的文本文件，四列：GSM3093315 SRR6985601 Homo_sapiens PAIRED
#第四个参数为输出的文件夹 GSE233609_HDMBO3_A485
treat_group=$1
con_group=$2
meta_file=$3
out_dir=$4
##分割 sample
##
echo "开始处理药物处理组样本..."
IFS=',' read -r -a parts <<< "$treat_group"
###extract fastq
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
		echo "样本 ${parts[i]} 处理完成"
	else
		echo "样本 ${parts[i]} 未在 Metainfo 文件中找到。"
	fi
done
##用最后一个样本得到 read length
read_len=$(cat ${gsm}/*html | grep -oP "mean length before filtering:.*?\K\d+bp" | head -1 | grep -oP '\d+')

echo "所有样本处理完成。"
###运行 rmats
mkdir /home/wt/GEO_data/SRA/gse_out/${out_dir}
mkdir /home/wt/GEO_data/SRA/gse_out/${out_dir}/rmats_out
mkdir /home/wt/GEO_data/SRA/gse_out/${out_dir}/tmp
##制作样本信息 txt 文件
IFS=',' read -ra samples <<< "$treat_group"
temp="/home/wt/GEO_data/SRA/out/sample/sample_Aligned.sortedByCoord.out.bam"
declare -a treat_samples
# 遍历原数组，替换并添加到新数组
for name in "${samples[@]}"; do
    new_str="${temp//sample/$name}"
    treat_samples+=("$new_str")
done
printf -v treat_samples '%s,' "${treat_samples[@]}"
treat_samples="${treat_samples%,}"  # 移除最后一个逗号和空格

IFS=',' read -ra samples <<< "$con_group"
declare -a control_samples
# 遍历原数组，替换并添加到新数组
for name in "${samples[@]}"; do
    new_str="${temp//sample/$name}"
    control_samples+=("$new_str")
done
printf -v control_samples '%s,' "${control_samples[@]}"
control_samples="${control_samples%,}"  # 移除最后一个逗号和空格

echo "${treat_samples}" > /home/wt/GEO_data/SRA/gse_out/${out_dir}/treat_samples.txt
echo "${control_samples}" > /home/wt/GEO_data/SRA/gse_out/${out_dir}/control_samples.txt
###
echo "运行 rmats..."
echo "读长为 ${read_len}"
if [[ "$pair" == "PAIRED" ]]; then is_pair="paired"; else is_pair="single"; fi
conda run -n rmats rmats.py --b1 /home/wt/GEO_data/SRA/gse_out/${out_dir}/treat_samples.txt --b2 /home/wt/GEO_data/SRA/gse_out/${out_dir}/control_samples.txt --gtf /home/wt/Genome_data/Homo_sapiens.GRCh38.115.gtf -t ${is_pair} --readLength ${read_len} --nthread 15 --od /home/wt/GEO_data/SRA/gse_out/${out_dir}/rmats_out/ --tmp /home/wt/GEO_data/SRA/gse_out/${out_dir}/tmp/ --variable-read-length
echo "rmats 运行完成"
echo "删除处理组 BAM 文件..."
IFS=',' read -r -a parts <<< "$treat_group"
###extract fastq
for ((i=0; i<${#parts[@]}; i++))
do
	rm /home/wt/GEO_data/SRA/out/${parts[i]}/*bam
	echo "已删除样本 ${parts[i]} 的 BAM 文件"
done
echo "运行完成"
