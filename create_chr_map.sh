output_file="chr_map.txt"

for i in {1..22}; do
    echo -e "chr$i\t$i" >> "$output_file"
done
