### Sync from deNBI cloud cluster to MMCI/MPI

Executed as daily cronjob at 6am on d3compute09

```bash
rsync --recursive --delete-before --prune-empty-dirs \
    --exclude-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/exclude.txt \
    --include-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/include.txt \
    centos@valet:/mnt/vol/beeond_backup/projects/diploid-assembly \
    /scratch/bioinf/projects/diploid-genome-assembly/sync/denbi
```
