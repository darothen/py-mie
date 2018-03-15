#!/usr/bin/env bash

git checkout master

echo "Building and deploying docs to gh-pages"
cd docs
make clean
make html
cd ..

# commit and push
git add -A
git commit -m "building and pushing new docs"
git push origin master

# switch branches
git checkout gh-pages
rm -rf *
touch .nojekyll
git chekcout master docs/_build/html
mv ./docs/_build/html/* ./
rm -rf ./docs

git add -A
git commit -m "publish new docs"
git push origin gh-pages

# switch back to the master branch
git checkout master
