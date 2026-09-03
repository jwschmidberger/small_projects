#!/usr/bin/env bash
# Detect tandem CMD cores, validate their internal similarity, and compare the
# intact duplication products with the two DehI halves.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
results_dir="${script_dir}/results"
project_root="$(cd "${script_dir}/../.." && pwd)"
python_cmd="${PYTHON:-${project_root}/.venv/bin/python}"

for program in mafft reformat.pl hhmake hhalign; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: ${program} is not on PATH." >&2
        exit 1
    fi
done

"${python_cmd}" "${script_dir}/identify_cmd_ahpd_duplications.py" prepare
mafft --thread -1 --retree 1 --maxiterate 0 \
    "${results_dir}/cmd_single_core_seed.fasta" \
    > "${results_dir}/cmd_single_core_seed_msa.fasta"
reformat.pl fas a3m "${results_dir}/cmd_single_core_seed_msa.fasta" \
    "${results_dir}/cmd_single_core_seed.a3m" -M 50
reformat.pl a3m sto "${results_dir}/cmd_single_core_seed.a3m" \
    "${results_dir}/cmd_single_core_seed.sto"
"${python_cmd}" "${script_dir}/identify_cmd_ahpd_duplications.py" scan

build_profile() {
    local stem="$1"
    local name="$2"
    mafft --thread -1 --retree 1 --maxiterate 0 \
        "${results_dir}/${stem}.fasta" > "${results_dir}/${stem}_msa.fasta"
    reformat.pl fas a3m "${results_dir}/${stem}_msa.fasta" \
        "${results_dir}/${stem}.a3m" -M 50
    hhmake -i "${results_dir}/${stem}.a3m" \
        -o "${results_dir}/${stem}.hhm" -name "${name}"
}

compare() {
    local query_stem="$1"
    local target_stem="$2"
    local output_stem="$3"
    hhalign -i "${results_dir}/${query_stem}.hhm" \
        -t "${results_dir}/${target_stem}.hhm" \
        -o "${results_dir}/${output_stem}.hhr" \
        -Oa3m "${results_dir}/${output_stem}.a3m" \
        -atab "${results_dir}/${output_stem}.tsv"
}

build_profile cmd_double_core CMD_double_core_intact
build_profile cmd_double_core_repeat1 CMD_double_core_repeat1
build_profile cmd_double_core_repeat2 CMD_double_core_repeat2
build_profile ahpd_double_core AhpD_double_core_intact
build_profile ahpd_double_core_repeat1 AhpD_double_core_repeat1
build_profile ahpd_double_core_repeat2 AhpD_double_core_repeat2
build_profile cmd_compact_120_150 CMD_compact_120_150_intact
build_profile cmd_compact_120_150_repeat1 CMD_compact_120_150_repeat1
build_profile cmd_compact_120_150_repeat2 CMD_compact_120_150_repeat2
build_profile ahpd_compact_120_150 AhpD_compact_120_150_intact
build_profile ahpd_compact_120_150_repeat1 AhpD_compact_120_150_repeat1
build_profile ahpd_compact_120_150_repeat2 AhpD_compact_120_150_repeat2

compare cmd_double_core_repeat1 cmd_double_core_repeat2 CMD_repeat1_vs_repeat2
compare ahpd_double_core_repeat1 ahpd_double_core_repeat2 AhpD_repeat1_vs_repeat2
compare cmd_compact_120_150_repeat1 cmd_compact_120_150_repeat2 CMD_compact_repeat1_vs_repeat2
compare ahpd_compact_120_150_repeat1 ahpd_compact_120_150_repeat2 AhpD_compact_repeat1_vs_repeat2
compare dehi_N_terminal cmd_double_core dehi_N_vs_CMD_double_core
compare dehi_C_terminal cmd_double_core dehi_C_vs_CMD_double_core
compare dehi_N_terminal ahpd_double_core dehi_N_vs_AhpD_double_core
compare dehi_C_terminal ahpd_double_core dehi_C_vs_AhpD_double_core
compare dehi_N_terminal cmd_compact_120_150 dehi_N_vs_CMD_compact
compare dehi_C_terminal cmd_compact_120_150 dehi_C_vs_CMD_compact
compare dehi_N_terminal ahpd_compact_120_150 dehi_N_vs_AhpD_compact
compare dehi_C_terminal ahpd_compact_120_150 dehi_C_vs_AhpD_compact

echo "Built duplication-aware CMD/AhpD profiles in ${results_dir}"
