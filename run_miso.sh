## [MISO](https://miso.readthedocs.io/en/fastmiso/)

source activate miso

source ./config.txt

caselist_string=`less $case_bam_list_file` ; caselist=(${caselist_string// / })
ctrllist_string=`less $ctrl_bam_list_file` ; ctrllist=(${ctrllist_string// / })
total_bam_list= (${caselist[@]} ${caselist[@]})



index_db=$output_dir/`basename $GFF3_anno|awk -F'.gff3' '{print$1}"_index_db"'`
exon_utils_dir=$output_dir/`basename $GFF3_anno|awk -F'.gff3' '{print$1"_exon_utils"}'`
insert_dist_dir=$output_dir/`basename $GFF3_anno|awk -F'.gff3' '{print$1"_insert_dist"}'`

mkdir -p $index_db
mkdir -p $insert_dist_dir
mkdir -p $exon_utils_dir

index_gff --index $GFF3_anno $index_db/ || echo index_gff failed !
exon_utils --get-const-exons $GFF3_anno --min-exon-size 1000 --output-dir $exon_utils_dir/ || echo exon_utils failed !

each_bam_miso(){
  bam_file=$1
  exon_utils_file=$2
  insert_dist_dir=$3
  each_outdir=$4
  read_length=$5

  settingfile=`miso|tail -n 1|awk -F': ' '{print$2}'`

  pe_utils --compute-insert-len $bam_file $exon_utils_file --output-dir $insert_dist_dir/
  mn=`head -n 1 $insert_dist_dir/$(basename $bam_file).insert_len| cut -f2 -d '='|cut -f1 -d ','|xargs  printf "%.*f" 0`
  sd=`head -n 1 $insert_dist_dir/$(basename $bam_file).insert_len| cut -f3 -d '='|cut -f1 -d ','|xargs  printf "%.*f" 0`
  miso --run $index_db  $bam_file \
    --output-dir $each_outdir/ \
    --read-len $readlength \
    --paired-end $mn $sd \
    --settings-filename $settingfile \
    -p 12 \
  summarize_miso --summarize-samples $each_outdir/ $each_outdir/
}


for bamfile in $total_bam_list
do 
  exon_utils_file=$exon_utils_dir/`ls $exon_utils_dir`
  each_outdir=$outdir/$(basename $bamfile|awk -F'.sort.bam' '{print$1}')
  each_bam_miso $bam_file $exon_utils_file $insert_dist_dir $each_outdir $read_length
done



bam_list=`find $bamdir -name "*.sort.bam"` ;\
bam_prefix_list=`for i in $bam_list;do echo $i|awk -F'/' '{print$NF}'|awk -F'.' '{print$1}' ;done`

for i in $bam_prefix_list
do
  pe_utils --compute-insert-len $bamdir/*/$i.sort.bam $(ls $exon_utils_dir) --output-dir $insert_dist_dir/
  mn=`head -n 1 $insert_dist_dir/$i.sort.bam.insert_len| cut -f2 -d '='|cut -f1 -d ','|xargs  printf "%.*f" 0`
  sd=`head -n 1 $insert_dist_dir/$i.sort.bam.insert_len| cut -f3 -d '='|cut -f1 -d ','|xargs  printf "%.*f" 0`
  outdir=`dirname $bamdir | xargs -i echo {}/result/$i `
  mkdir -p $outdir
  echo miso --run $index_db  $bamdir/*/$i.sort.bam --output-dir $outdir/ --read-len $readlength --paired-end $mn $sd --settings-filename $settingfile -p 12 --prefilter
  miso --run $index_db  $bamdir/*/$i.sort.bam \
    --output-dir $outdir/ \
    --read-len $readlength \
    --paired-end $mn $sd \
    --settings-filename $settingfile \
    -p 12 \
  summarize_miso --summarize-samples $outdir/ $outdir/

done

compare_miso --compare-samples control/ knockdown/ comparisons/
