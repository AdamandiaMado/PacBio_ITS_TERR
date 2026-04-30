echo -e "file\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > seqkit_summary.tsv
for f in *hifi_reads; do
    tail -n 1 "$f" | awk -v file="$(basename "$f")" '{
        gsub(",", "");
        print file"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
    }' >> seqkit_summary.tsv
done
