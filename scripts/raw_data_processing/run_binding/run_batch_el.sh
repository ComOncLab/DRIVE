cat all_splits | while read i
do
	split_name=$(basename $i .fasta)
	tmux new-session -d -s el_${split_name} "
		conda activate iedb
		cd /home/wt/cy/cell_pipeline/MC38_invitro/A_C/jcast_out/binding_out
		sed "s/test/${split_name}/g" el.json > el_${split_name}.json
		python3 /home/wt/software/ng_tc1-0.1.2-beta/src/tcell_mhci.py -j el_${split_name}.json -o res/el_${split_name} -f tsv
	"
	echo "el_${split_name} 已启动"
done

