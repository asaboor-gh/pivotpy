Write-Host "Running Notebooks..." -ForegroundColor Green

Resolve-Path "*.ipynb" | ForEach-Object { jupyter nbconvert --execute --to notebook --inplace $_}

Write-Host "Building Lib/Docs..." -ForegroundColor Green
nbdev_build_lib.exe 
nbdev_build_docs.exe

Write-Host "Running Notebooks and Stripping Output..." -ForegroundColor Green
Resolve-Path "*[O,s].ipynb" | ForEach-Object { jupyter nbconvert --execute --clear-output $_}

Write-Host "Simply Cleaning Now..." -ForegroundColor Green
nbdev_install_git_hooks.exe
nbdev_clean_nbs.exe #All clean