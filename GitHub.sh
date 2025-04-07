git pull
git add --all
git commit -m "updating package"
git push
find . -name .DS_Store -print0 | xargs -0 git rm -f --ignore-unmatch
find . -name .Rapp.history -print0 | xargs -0 git rm -f --ignore-unmatch
git add --all
git commit -m "updating package"
git push
