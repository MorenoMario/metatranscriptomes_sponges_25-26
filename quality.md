1. 
```
for fa in Trinity_Dendrilla_out.Trinity.fasta Trinity_Mycale_out.Trinity.fasta Trinity_Myxilla_out.Trinity.fasta
do
  $TRINITY_HOME/util/TrinityStats.pl $fa > ${fa}.basic_stats.txt
done
```
