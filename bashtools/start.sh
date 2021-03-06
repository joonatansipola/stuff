#
# BASH CONFIG
#

# source bashrc
#if [ -f ~/.bashrc ]; then
#	. ~/.bashrc
#fi

# misc
#PS1='\[\033[01;34m\][\W] \[\033[01;32m\]€\[\033[00m\] '
alias lsa="ls -Alh"
shopt -s globstar

# history
#HISTSIZE=1000
#HISTFILESIZE=2000
HISTCONTROL=ignoredups:erasedups
shopt -s histappend
PROMPT_COMMAND="history -n; history -w; history -c; history -r; $PROMPT_COMMAND"

#
# GITHUB CONFIG
#

# git graph
git config --global alias.graph "log --all --graph --decorate --oneline"

# quick wget from joonatansipola github
jsgget() {
        f=https://raw.githubusercontent.com/joonatansipola/stuff/master/$1
        wget $f
}

# quick add-commit-push
gitup () {
        git reset --
        git add $1
        git commit -m'gitup push'
        git push
}

# quick github upload
#upstuff () {
#
#        wd=$(pwd)
#        cd <LOCAL_REPO_PATH>
#
#        fname=$(basename $1)
#        fpath=<FILE_BASE_PATH>/$1
#
#	cp -i $fpath $2/$fname
#
#	git reset --
#	git add $2/$fname
#	git commit -m'upstuff upload'
#	git push
#	
#	echo
#        ls
#        cd $wd
#}
