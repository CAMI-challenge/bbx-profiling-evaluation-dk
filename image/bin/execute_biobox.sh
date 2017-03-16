#!/bin/bash  

set -o errexit
set -o nounset
set -o xtrace

METADATA=/bbx/metadata

GROUND_TRUTH=$(biobox_args.sh 'select(has("ground_truth")) | .ground_truth | .value ')

PREDICTION=$(biobox_args.sh 'select(has("prediction")) | .prediction | .value ')

CMD=$(fetch_task_from_taskfile.sh ${TASKFILE} $1)

OUTPUT_FILE="${OUTPUT}/metrics.txt"

eval $CMD

cat << EOF > ${OUTPUT}/biobox.yaml
version: 0.1.1
results:
  - name: David Koslicki
    type: tsv
    inline: false
    description: This tool produces the metrics True Positives, False Positives, False Negatives, Precision, Sensitivity, L1-Norm, Unifrac
    value: metrics.txt
EOF
