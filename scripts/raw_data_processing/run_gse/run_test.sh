while IFS=" " read -r control
do
	bash process_control_match.sh ${control} ./rem_gse ./meta_info.txt	
done < "test"
