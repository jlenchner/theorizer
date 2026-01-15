#!/bin/bash

# SynPAT Dataset Generation Script
# This script runs the generate_dataset.py with commonly used parameters.
# Modify the parameters below as needed for your use case.

# Exit on error
set -e

# Default parameters - modify these as needed
VARS="['m1', 'm2', 'd1', 'd2', 'W', 'theta', 'sinTheta', 'p', 'Fc', 'Fg']"
DERIVS="['dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2']"
NUM_VARS="[6,7,8,9]"
NUM_DERIVS="[2,3,4]"
NUM_EQUATIONS="[4,5]"
NUM_REPLACEMENTS=5
NUM_SYSTEMS=3
NUM_CONST_CONSEQ=1
NUM_DERIV_CONSEQ=1
CONSEQ_DATA_RANGE="[1, 10]"
SYS_DATA_RANGE="[1, 10]"
MAX_CONSEQ_TERMS=8
TIMEOUT=3600
SEED=42

# Generation flags - set to True or False
GEN_REPLACEMENTS=True
GEN_SYS_DATA=True
GEN_CONSEQUENCE=True
GEN_CONSEQUENCE_DATA=True

echo "=========================================="
echo "SynPAT Dataset Generation"
echo "=========================================="
echo "Variables: $VARS"
echo "Derivatives: $DERIVS"
echo "Number of variables options: $NUM_VARS"
echo "Number of derivatives options: $NUM_DERIVS"
echo "Number of equations options: $NUM_EQUATIONS"
echo "Number of systems per config: $NUM_SYSTEMS"
echo "Max constants in consequence: $NUM_CONST_CONSEQ"
echo "Max derivatives in consequence: $NUM_DERIV_CONSEQ"
echo "=========================================="

python generate_dataset.py \
    --vars "$VARS" \
    --derivs "$DERIVS" \
    --numVars "$NUM_VARS" \
    --numDerivs "$NUM_DERIVS" \
    --numEquations "$NUM_EQUATIONS" \
    --numReplacements $NUM_REPLACEMENTS \
    --numSystems $NUM_SYSTEMS \
    --numConstConseq $NUM_CONST_CONSEQ \
    --numDerivConseq $NUM_DERIV_CONSEQ \
    --conseqDataRange "$CONSEQ_DATA_RANGE" \
    --sysDataRange "$SYS_DATA_RANGE" \
    --maxConseqTerms $MAX_CONSEQ_TERMS \
    --timeout $TIMEOUT \
    --seed $SEED \
    --genReplacements $GEN_REPLACEMENTS \
    --genSysData $GEN_SYS_DATA \
    --genConsequence $GEN_CONSEQUENCE \
    --genConsequenceData $GEN_CONSEQUENCE_DATA

echo "=========================================="
echo "Generation complete!"
echo "=========================================="
