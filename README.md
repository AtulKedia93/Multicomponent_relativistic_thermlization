# Multicomponent_relativistic_thermlization

# README #

This is a repostitory by University of Notre Dame PhD students Atul Kedia and Nishanth Sasankan

Here you will find a lot of codes and data that look very similar to each other. As a general rule, do not use
codes/data from here that are older than Dec-2019. Those are legacy codes used to make a more current version
which has a similar name. There might be some repetition in newer codes as well, but for all, in general go 
for the most recent version.



#creating a nice readme edit and push request.

This repository "should" contain all the codes we write in the cosmology group for solving the Cosmic Lithium problem.
Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)


How to use Git (first-time)
===========================

1. Use git clone to create a local copy of the repo on your machine. The cmd
should look like: 
git clone https://username@bitbucket.org/Atul_Kedia/bbn_pmf.git

2. Make a new branch using
git branch new_branch_name

3. Switch current branch to new_branch_name using
git checkout new_branch_name

4. Make necessary edits to files. Check compilation and output PDF.

5. IMP: Ensure that only the source code files have been edited/tracked using
git status (otherwise contact Atul/read .gitignore section from Pro Git)

6. Save all your changes with a useful message
git commit -a -m "helpful message here"

7. Push to the web repo using
git push --set-upstream origin new_branch_name

8. Use the web interface for creating a pull request (since this is easier)
* Branches (left bar) 
* Check that new_branch_name is ahead by one or two commits  (from Filter: Active) 
* Select/Click on new_branch_name
* Create pull request
* Check that new_branch_name is being merged into master
* Write appropriate description and title.
* Enter other person's username in Reviewers.

9. Other person will do merging using web interface after checking compilation.
* Run "git fetch". Similar to "git clone". Clone => install, fetch => update.
* git checkout origin/new_branch_name
* Check compilation without errors. Also check whether output looks okay.
* Suggest changes using whatever method, if needed.
* If everything is okay, go to Pull requests (left bar on website)
* Click Merge option
* Double-check that "Close source branch" option is selected.
* Merge

How to use Git (second-time onwards)
====================================

1. Use "git pull" instead of "git clone ...".

2. Repeat the steps as earlier.

Hints
=====
* You can always see what is the current branch using "git branch --list"# Multicomponent_relativistic_thermlization