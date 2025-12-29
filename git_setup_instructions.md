## Git Setup Instructions

### 1. Local Setup (MacBook)


# Make folder easier to get to—call it "docs"
nano ~/.zshrc
# add this line below to this zsh file
ln -s "$HOME/Library/Mobile Documents/com~apple~CloudDocs/Documents - Dylan's MacBook Air" ~/docs
# exit file and run below line
source ~/.zshrc


```bash
cd ~/docs
cd data_local
```

# Pull files from Github
git init
git remote set-url origin <your-github-repo-url>.git
git pull origin main

# might need to remove a file for some reason
rm .DS_Store



### Cluster Setup (Pod/Bridges-2)
```bash
# SSH to cluster
ssh pod-login1.cnsi.ucsb.edu

cd ~/Documents
mkdir -p lammps_data/input_data

# Clone repository
git clone <your-github-repo-url>
cd lammps_work
chmod +x scripts/run_lammps.sh

# Upload specific .data files you need
# Use scp or rsync from your Mac:
# scp ~/Documents/lammps_data/input_data/slab_*.data user@pod:~/Documents/lammps_data/input_data/
```

# For VS Code modifications: commit and push
git add .
git commit -m "Initial commit"
git push -u origin main
```


