#!/bin/bash
# Push to GitHub - Replace YOUR_USERNAME and YOUR_REPO_NAME with your actual values

# Add your GitHub repository as origin
# Replace YOUR_USERNAME with your GitHub username
# Replace YOUR_REPO_NAME with the name you gave your repository
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git

# Push to GitHub
git push -u origin main

echo "Repository pushed to GitHub!"
echo "Don't forget to replace YOUR_USERNAME and YOUR_REPO_NAME in this script!"