### Sync from deNBI cloud cluster to MMCI/MPI

Executed as daily cronjob at 6am on d3compute09

```bash
rsync --recursive --delete-before --prune-empty-dirs \
    --exclude-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/exclude.txt \
    --include-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/include.txt \
    centos@valet:/mnt/vol/beeond_backup/projects/diploid-assembly \
    /scratch/bioinf/projects/diploid-genome-assembly/sync/denbi
```

Executed as daily cronjob at 6pm on lap-13-72

```bash
rsync --recursive --delete-before --prune-empty-dirs \
    -e "ssh contact.mpi-inf.mpg.de ssh" \
    --exclude-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/exclude.txt \
    --include-from=/home/pebert/work/code/github/project-diploid-assembly/environment/sync/include.txt \
    /mnt/sshfs/hhu/project/ebertp/projects/rfdga \
    /scratch/bioinf/projects/diploid-genome-assembly/sync/hhu
```
