#!/bin/bash
#BSUB -J analyze_data
#BSUB -q cpuqueue
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/analyze_data_%J.stdout
#BSUB -eo logs/analyze_data_%J.stderr

# Initialize conda environment
source ~/.bashrc
conda activate rnacompete
cd $LS_SUBCWD

# HybID00025_00103
# HybID00110_00146
# HybID00149_00243
# HybID00256_00295
# HybID00298_00386
# HybID00389_00410
# HybID00428_00526
# HybID00537_00754
# HybID00763_00836
# HybID00841_00866
# HybID00899_01195
# HybID01196_01496
# HybID01497_01618
# HybID01619_01828
# HybID01829_02042
# HybID02043_02210
# HybID02211_02318
# HybID02323_02479

# Run pipeline
python3 analyze_data.py HybID00025_00035
