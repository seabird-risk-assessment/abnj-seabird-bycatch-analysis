find .. -type f -name "*.r"  | xargs grep input-data > inputs.txt
sed -i 's|.*/input-data/||' inputs.txt
sed -i 's|,.*$||' inputs.txt
sed -i "s|')$||" inputs.txt
sed -i "s|'$||" inputs.txt
sed -i 's|"$||' inputs.txt
sort inputs.txt | uniq > file-names.txt
rm inputs.txt

