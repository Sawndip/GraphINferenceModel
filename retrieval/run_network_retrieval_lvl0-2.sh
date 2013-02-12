#if [ $# -lt 1 ]; then
#  echo "Usage: run_network_retrieval_lvl0-2.sh <param-file> [<qrel>]"
#  exit
#fi

params=params/cons-medical-topics101-185-active.params
if [ $# -gt 0 ]; then
  params=$1
fi
qrel=params/medtrack-all.qrel
if [ $# -gt 1 ]; then
  qrel=$2
fi

echo params: $params qrel: $qrel

for i in `seq 0 1 2`; do
	./network_retrieval $params $i 
	base=$(echo $params | sed 's/.params//g')
	cp ${base}.results $params.$i.results
	trec_eval -q $qrel ${base}.results > $params.$i.eval
done
paste $params.*.eval > ${base}.eval
cat ${base}.eval