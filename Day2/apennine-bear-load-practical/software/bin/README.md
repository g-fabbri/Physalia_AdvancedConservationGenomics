# Day 2 command links

The Conda environment installs `minimap2`, `transanno`, `liftOver`, and `bigWigToBedGraph`. After activating it, run:

```bash
bash software/link_conda_tools.sh
```

The script creates local symbolic links in this directory. Consequently, the lessons can consistently use `software/bin/` while the package versions remain managed by Conda. Do not commit machine-specific generated links.
