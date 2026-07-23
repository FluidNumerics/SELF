#!/bin/bash
# Sync the local working tree to the MI300X cluster login node so the in-container
# build job (build_gfx942_sqsh.slurm) can copy it over /opt/self/src.
#
# Usage:  bench/mi300x/sync_source.sh [CLUSTER] [REMOTE_SRC]
#   CLUSTER     ssh target (default: jschoonover@iad-mj-login)
#   REMOTE_SRC  destination directory on the cluster (default: ~/self-src)
#
# Run from the repository root. Excludes local build trees, git internals, and
# large container/output artifacts that must not be shipped.
set -euo pipefail

CLUSTER="${1:-jschoonover@iad-mj-login}"
REMOTE_SRC="${2:-self-src}"
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

echo "Syncing ${REPO_ROOT}/ -> ${CLUSTER}:~/${REMOTE_SRC}/"
rsync -a --delete \
  --exclude '.git/' \
  --exclude 'build*/' \
  --exclude 'install*/' \
  --exclude '*.sqsh' \
  --exclude 'output/' \
  --exclude 'test_output/' \
  --exclude '*-results/' \
  "${REPO_ROOT}/" "${CLUSTER}:${REMOTE_SRC}/"
echo "Done. On the cluster: sbatch bench/mi300x/build_gfx942_sqsh.slurm"
